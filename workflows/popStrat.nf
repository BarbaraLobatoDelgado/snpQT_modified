// Population stratification workflow
nextflow.enable.dsl = 2

// import modules
include {filter_maf} from '../modules/popStrat.nf' // D3
include {run_snpflip} from '../modules/popStrat.nf' // D4 
include {flip_snps} from '../modules/popStrat.nf' // D4
include {remove_ambiguous_snps} from '../modules/popStrat.nf' // D4b
include {align} from '../modules/popStrat.nf' // D5
include {merge} from '../modules/popStrat.nf' // D6
include {pca_prep} from '../modules/popStrat.nf' // D6
include {popfile} from '../modules/popStrat.nf' // D7
include {eigensoft} from '../modules/popStrat.nf' // D8
include {pca_plink} from '../modules/popStrat.nf' // D8
include {plot_plink_pca} from '../modules/popStrat.nf' // D8
include {extract_homogenous} from '../modules/popStrat.nf' // D9
include {parse_logs} from '../modules/qc.nf' // E13
include {report} from '../modules/qc.nf' // E14

// workflow component for snpqt pipeline
workflow pop_strat {
  take:
    ch_bed
    ch_bim
    ch_fam

  main:
    Channel
      .fromPath("${params.db}/all_phase3_1.bed", checkIfExists: true)
      .set{ ref_bed }
    Channel
      .fromPath("${params.db}/all_phase3_1.bim", checkIfExists: true)
      .set{ ref_bim }
    Channel
      .fromPath("${params.db}/all_phase3_1.fam", checkIfExists: true)
      .set{ ref_fam } 
    Channel
      .fromPath("${params.db}/PCA.exclude.regions.b37.txt", checkIfExists: true)
      .set{ exclude }
    filter_maf(ch_bed, ch_bim, ch_fam, exclude)
    Channel
      .fromPath("${params.db}/h37_squeezed.fasta", checkIfExists: true)
      .set { g37 }
    Channel
      .fromPath("${params.db}/hg38.fa", checkIfExists: true)
      .set { g38 }
    run_snpflip(filter_maf.out.bed, filter_maf.out.bim, filter_maf.out.fam, g37, g38)
    flip_snps(filter_maf.out.bed, filter_maf.out.bim, filter_maf.out.fam, run_snpflip.out.rev) 
    remove_ambiguous_snps(flip_snps.out.bed, flip_snps.out.bim, flip_snps.out.fam, run_snpflip.out.ambig)
    align(remove_ambiguous_snps.out.bed, remove_ambiguous_snps.out.bim, remove_ambiguous_snps.out.fam, ref_bed, ref_bim, ref_fam)
    merge(align.out.bed, align.out.bim, align.out.fam, ref_bed, ref_bim, ref_fam)
    pca_prep(merge.out.bed, merge.out.bim, merge.out.fam, exclude)
    Channel
      .fromPath("${params.db}/integrated_call_samples_v3.20130502.ALL.panel", checkIfExists: true)
      .set{ panel }
    popfile(panel)
    if (params.popfile == "super") {
      rf = popfile.out.super
    } else if (params.popfile == "sub") {
      rf = popfile.out.sub
    }
    if (params.parfile == false ) {
      Channel
        .fromPath("$baseDir/bootstrap/default_parfile", checkIfExists: true)
        .set{ parfile_ch }
    } else {
      Channel
        .fromPath(params.parfile, checkIfExists: true)
        .set{ parfile_ch }
    }    
    eigensoft(pca_prep.out.bed, pca_prep.out.bim, pca_prep.out.fam, rf, filter_maf.out.fam, parfile_ch)
    pca_plink(pca_prep.out.bed, pca_prep.out.bim, pca_prep.out.fam, eigensoft.out.eigenvec)
    // prepare plink pca files for plotting
    // want list of tuples: id, eigenvec, popfile
    after_ch = Channel.from("after").concat(pca_plink.out.after, rf).toList()
    plot_plink_pca(Channel.from("before").concat(pca_plink.out.before, rf).toList().concat(after_ch))
    extract_homogenous(ch_bed, ch_bim, ch_fam, eigensoft.out.keep_samples)
    logs = filter_maf.out.log.concat(flip_snps.out.log, align.out.log, merge.out.log, pca_prep.out.log, extract_homogenous.out.log).collect()
    parse_logs("pop_strat", logs, "pop_strat_log.txt")

    figures = plot_plink_pca.out.figure
      .concat(plot_plink_pca.out.rds, parse_logs.out.figure)
      .collect()
    Channel
      .fromPath("$baseDir/bootstrap/popstrat_report.Rmd", checkIfExists: true)
      .set{ rmd }
    // // Pass yaml file containig parameters for nextflow run
    // yaml_file_path = paramsFile
    // Channel
    //   .path(params.yaml_file_path)
    //   .set { yaml_file }
    report("pop_strat", figures, rmd, "${params}")   

  emit:
    bed = extract_homogenous.out.bed
    bim = extract_homogenous.out.bim
    fam = extract_homogenous.out.fam
} 
