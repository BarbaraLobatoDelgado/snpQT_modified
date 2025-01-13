// Post-imputation QC workflow after imputation in TOPMed server
nextflow.enable.dsl = 2

// import modules
// include {filter_imp} from '../modules/postImputation.nf' // H1
include {filter_maf} from '../modules/postImputation.nf' // H2
include {duplicates_cat1} from '../modules/postImputation.nf' // H3
include {duplicates_cat2} from '../modules/postImputation.nf' // H4
include {duplicates_cat3} from '../modules/postImputation.nf' // H5
include {update_ids} from '../modules/postImputation.nf' // H6
include {update_phenotype} from '../modules/postImputation.nf' // H7
include {parse_logs} from '../modules/qc.nf' // E13

// specific for TOPMed imputation
include {unzip_vcfs} from '../modules/postImputation_topmed.nf'
include {create_index} from '../modules/postImputation_topmed.nf'
include {concat_vcfs} from '../modules/postImputation_topmed.nf'
include {sort_vcf} from '../modules/postImputation_topmed.nf'
include {trim_ids} from '../modules/postImputation_topmed.nf'
include {filter_imp_quality_topmed} from '../modules/postImputation_topmed.nf'


workflow postImputation_topmed {
  take:
  qc_fam
  topmed_results_dir

  main:
  unzip_vcfs(topmed_results_dir)
  create_index(unzip_vcfs.out.vcf)
  concat_vcfs(create_index.out.csi, unzip_vcfs.out.vcf)
  sort_vcf(concat_vcfs.out.concat_vcf)
  trim_ids(sort_vcf.out.sorted_vcf)
  filter_imp_quality_topmed(trim_ids.out.correct_ids_vcf)
  filter_maf(filter_imp_quality_topmed.out.bed, filter_imp_quality_topmed.out.bim, filter_imp_quality_topmed.out.fam)
  duplicates_cat1(filter_maf.out.bed, filter_maf.out.bim, filter_maf.out.fam)
  duplicates_cat2(duplicates_cat1.out.bed, duplicates_cat1.out.bim, duplicates_cat1.out.fam)
  duplicates_cat3(duplicates_cat2.out.bed, duplicates_cat2.out.bim, duplicates_cat2.out.fam)
  update_phenotype(duplicates_cat3.out.bed, duplicates_cat3.out.bim, duplicates_cat3.out.fam, qc_fam)
  // Logs
  logs = filter_imp_quality_topmed.out.log.concat(filter_maf.out.log, duplicates_cat1.out.log, duplicates_cat2.out.log, duplicates_cat3.out.log, update_phenotype.out.log).collect()
  parse_logs("post_imputation_topmed", logs, "post_impute_topmed_log.txt")

  emit:
  bed = update_phenotype.out.bed
  bim = update_phenotype.out.bim
  fam = update_phenotype.out.fam
}