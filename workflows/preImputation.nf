// Imputation workflow
nextflow.enable.dsl = 2

// import modules
include {set_chrom_code} from '../modules/imputation.nf' // F1
include {run_snpflip} from '../modules/popStrat.nf' // F2, reuse D4
include {flip_snps} from '../modules/popStrat.nf' // F2, reuse D4
include {remove_ambiguous_snps} from '../modules/popStrat.nf' // F2b
include {fix_duplicates} from '../modules/imputation.nf' // F3
include {create_afreq_file} from '../modules/imputation.nf' // F3b
include {run_mccarthy_tools} from '../modules/imputation.nf' // F3c
include {prepare_VCFs_for_TOPMed_imputation} from '../modules/imputation.nf' // F3d
include {to_vcf} from '../modules/imputation.nf' // F4
include {to_bcf} from '../modules/imputation.nf' // F5
include {check_ref_allele} from '../modules/imputation.nf' // F6
include {bcf_to_vcf} from '../modules/imputation.nf' // F7
include {parse_logs} from '../modules/qc.nf' // E13

// workflow component for snpqt pipeline
workflow preImputation {
  take:
    ch_bed
    ch_bim
    ch_fam

  main:
    set_chrom_code(ch_bed, ch_bim, ch_fam)
    Channel
      .fromPath("${params.db}/h37_squeezed.fasta", checkIfExists: true)
      .set { g37 }
    Channel
      .fromPath("${params.db}/hg38_modified.fa", checkIfExists: true)
      .set { g38 }

    // run_snpflip(set_chrom_code.out.bed, set_chrom_code.out.bim, set_chrom_code.out.fam, g37, g38)
    // flip_snps(ch_bed, ch_bim, ch_fam, run_snpflip.out.rev) // , run_snpflip.out.ambig
    // remove_ambiguous_snps(flip_snps.out.bed, flip_snps.out.bim, flip_snps.out.fam, run_snpflip.out.ambig)
    // fix_duplicates(remove_ambiguous_snps.out.bed, remove_ambiguous_snps.out.bim, remove_ambiguous_snps.out.fam)
    // fix_duplicates(flip_snps.out.bed, flip_snps.out.bim, flip_snps.out.fam)

    // Prepare files for local imputation
    if (params.pre_impute == true && params.pre_impute_TOPMed == false) {
      // Run SNPflip
      run_snpflip(set_chrom_code.out.bed, set_chrom_code.out.bim, set_chrom_code.out.fam, g37, g38)
      flip_snps(ch_bed, ch_bim, ch_fam, run_snpflip.out.rev, run_snpflip.out.ambig)
      // Remove duplicated SNPs
      fix_duplicates(flip_snps.out.bed, flip_snps.out.bim, flip_snps.out.fam)
      to_vcf(fix_duplicates.out.bed, fix_duplicates.out.bim, fix_duplicates.out.fam)
      to_bcf(to_vcf.out.vcf)
      Channel
        .fromPath("${params.db}/All_20180418.vcf.gz", checkIfExists: true) // file previously used is called SNPs_g38
        .set{ dbsnp38 }
      Channel
        .fromPath("${params.db}/All_20180423.vcf.gz", checkIfExists: true) 
        .set{ dbsnp37 }
      Channel
        .fromPath("${params.db}/All_20180418.vcf.gz.tbi", checkIfExists: true) 
        .set{ dbsnp_idx38 }
      Channel
        .fromPath("${params.db}/All_20180423.vcf.gz.tbi", checkIfExists: true)
        .set{ dbsnp_idx37 }
      check_ref_allele(to_bcf.out.bcf, g37, g38, dbsnp37, dbsnp38, dbsnp_idx38, dbsnp_idx37) // dbsnp, dbsnp_idx,
      bcf_to_vcf(check_ref_allele.out.bcf)
      // logs = fix_duplicates.out.log.concat(to_vcf.out.log).collect()
      // parse_logs("pre_imputation", logs, "pre_imputation_log.txt")
      emit:
      vcf = bcf_to_vcf.out.vcf
    }

    // Prepare files for TOPMed imputation server 
    if (params.pre_impute == true && params.pre_impute_TOPMed == true) {
      // Remove duplicated SNPs
      fix_duplicates(set_chrom_code.out.bed, set_chrom_code.out.bim, set_chrom_code.out.fam)
      // Not necessary to run SNPflip because McCarthy Group Tools takes care of it
      // ...
      // Create allele frequency file
      create_afreq_file(fix_duplicates.out.bed, fix_duplicates.out.bim, fix_duplicates.out.fam)
      Channel
        .fromPath("$baseDir/bin/mccarthy_tools_preimputation/HRC-1000G-check-bim.pl", checkIfExists: true)
        .set{ mccarthy_tools_script }
      Channel
        .fromPath("${params.db}/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz", checkIfExists: true)
        .set{ topmed_ref_hg38 }
      // run_mccarthy_tools(set_chrom_code.out.bed, set_chrom_code.out.bim, set_chrom_code.out.fam, create_afreq_file.out.afreq, mccarthy_tools_script, topmed_ref_hg38)
      run_mccarthy_tools(fix_duplicates.out.bed, fix_duplicates.out.bim, fix_duplicates.out.fam, create_afreq_file.out.afreq, mccarthy_tools_script, topmed_ref_hg38)
      prepare_VCFs_for_TOPMed_imputation(run_mccarthy_tools.out.chr_vcf)

      // to_vcf(fix_duplicates.out.bed, fix_duplicates.out.bim, fix_duplicates.out.fam)
      // to_bcf(to_vcf.out.vcf)
      // Channel
      //   .fromPath("${params.db}/All_20180418.vcf.gz", checkIfExists: true) // file previously used is called SNPs_g38
      //   .set{ dbsnp38 }
      // Channel
      //   .fromPath("${params.db}/All_20180423.vcf.gz", checkIfExists: true) 
      //   .set{ dbsnp37 }
      // Channel
      //   .fromPath("${params.db}/All_20180418.vcf.gz.tbi", checkIfExists: true) 
      //   .set{ dbsnp_idx38 }
      // Channel
      //   .fromPath("${params.db}/All_20180423.vcf.gz.tbi", checkIfExists: true)
      //   .set{ dbsnp_idx37 }
      // check_ref_allele(to_bcf.out.bcf, g37, g38, dbsnp37, dbsnp38, dbsnp_idx38, dbsnp_idx37) // dbsnp, dbsnp_idx,
      // bcf_to_vcf(check_ref_allele.out.bcf)
      // logs = fix_duplicates.out.log.concat(to_vcf.out.log).collect()
      // parse_logs("pre_imputation", logs, "pre_imputation_log.txt")
      emit:
      vcf_gz = prepare_VCFs_for_TOPMed_imputation.out.vcf_gz
      vcf_gz_tbi = prepare_VCFs_for_TOPMed_imputation.out.vcf_gz_tbi
    }
}  

