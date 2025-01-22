// STEP 1: exclude patients that are non-eligible, not signed IC, etc.
process select_eligible_patients {
  label 'plink1'

  input:
  path(bed)
  path(bim)
  path(fam)
  path(patients_list)

  output:
  path "eligible_patients.bed", emit: bed
  path "eligible_patients.bim", emit: bim
  path "eligible_patients.fam", emit: fam

  script:
  '''
  plink1.9 --bfile !{bim.baseName} --keep !{patients_list} --make-bed --out eligible_patients
  '''
}

// STEP 2: do SNP pruning to remove correlated SNPs before GWAS
process snp_pruning {
  label 'plink1'

  input:
  path(bed)
  path(bim)
  path(fam)

  output:
  path "pruned_data.bed", emit: bed
  path "pruned_data.bim", emit: bim
  path "pruned_data.fam", emit: fam

  script:
  '''
  plink1.9 --bfile !{bim.baseName} --indep-pairwise !{params.indep_pairwise} --out pruning_output 
  plink1.9 --bfile !{bim.baseName} --extract pruning_output.prune.in --make-bed --out pruned_data
  '''
}


// STEP 3: 
process run_gwas_lmm {
  publishDir "${params.results}/gwas_lmm/", mode: 'copy'

  input:
  path(bed)
  path(bim)
  path(fam)

  // output:
  // path "gwas_report.html"

  script:
  '''
  GENESIS_GWAS.R !{bed} !{bim} !{fam}
  '''
}

