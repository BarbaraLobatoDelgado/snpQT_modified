// GWAS workflow
nextflow.enable.dsl = 2

// import modules
include {select_eligible_patients} from '../modules/gwas_lmm.nf'
include {snp_pruning} from '../modules/gwas_lmm.nf'
include {run_gwas_lmm} from '../modules/gwas_lmm.nf'

workflow gwas_lmm {
  take:
  clean_bed
  clean_bim
  clean_fam
  patients_list

  main:
    select_eligible_patients(clean_bed, clean_bim, clean_fam, patients_list)
    snp_pruning(select_eligible_patients.out.bed, select_eligible_patients.out.bim, select_eligible_patients.out.fam)
    run_gwas_lmm(snp_pruning.out.bed, snp_pruning.out.bim, snp_pruning.out.fam, patients_list)
}

