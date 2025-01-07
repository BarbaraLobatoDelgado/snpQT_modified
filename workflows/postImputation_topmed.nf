// Post-imputation QC workflow
nextflow.enable.dsl = 2

// import modules
// include {filter_imp} from '../modules/postImputation.nf' // H1
// include {filter_maf} from '../modules/postImputation.nf' // H2
// include {duplicates_cat1} from '../modules/postImputation.nf' // H3
// include {duplicates_cat2} from '../modules/postImputation.nf' // H4
// include {duplicates_cat3} from '../modules/postImputation.nf' // H5
// include {update_ids} from '../modules/postImputation.nf' // H6
// include {update_phenotype} from '../modules/postImputation.nf' // H7
// include {parse_logs} from '../modules/qc.nf' // E13

include {unzip_chr_files} from '../modules/postImputation_topmed.nf'

workflow postImputation_topmed {
  take:
  topmed_results_dir

  main:
  unzip_chr_files(topmed_results_dir)

  emit:
  extracted_files_dose = unzip_chr_files.out.dose
  extracted_files_info = unzip_chr_files.out.info
}