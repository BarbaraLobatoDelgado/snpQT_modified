// Post-imputation QC workflow after imputation in TOPMed server
nextflow.enable.dsl = 2

// import modules
include {filter_imp} from '../modules/postImputation.nf' // H1
include {filter_maf} from '../modules/postImputation.nf' // H2
include {duplicates_cat1} from '../modules/postImputation.nf' // H3
include {duplicates_cat2} from '../modules/postImputation.nf' // H4
include {duplicates_cat3} from '../modules/postImputation.nf' // H5
include {update_ids} from '../modules/postImputation.nf' // H6
include {update_phenotype} from '../modules/postImputation.nf' // H7
include {parse_logs} from '../modules/qc.nf' // E13

include {unzip_vcfs} from '../modules/postImputation_topmed.nf'
// include {gunzip_vcf_files} from '../modules/postImputation_topmed.nf'
include {create_tbi_files} from '../modules/postImputation_topmed.nf'
include {concat_vcfs} from '../modules/postImputation_topmed.nf'
include {convert_vcfs_to_plink} from '../modules/postImputation_topmed.nf'
include {merge_plink_files} from '../modules/postImputation_topmed.nf'

workflow postImputation_topmed {
  take:
  topmed_results_dir

  main:
  unzip_vcfs(topmed_results_dir)
  create_tbi_files(unzip_vcfs.out.vcf)
  concat_vcfs(create_tbi_files.out.csi, unzip_vcfs.out.vcf)

  // gunzip_vcf_files(unzip_vcfs.out.dose)
  // convert_vcfs_to_plink(unzip_vcfs.out.vcf)
  // merge_plink_files(convert_vcfs_to_plink.out.bim, convert_vcfs_to_plink.out.bed, convert_vcfs_to_plink.out.fam)
  // filter_imp(convert_vcfs_to_plink.out.bim, convert_vcfs_to_plink.out.bed, convert_vcfs_to_plink.out.fam)
  // filter_imp(merge_vcfs.out.vcf_gz)



  // emit:
  // extracted_files_dose = unzip_vcfs.out.dose
  // extracted_files_info = unzip_vcfs.out.info
}