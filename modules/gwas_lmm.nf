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





// Step I1: Run logistic regression, adjusting for covariates
process run_gwas {    
    label 'plink2'
    publishDir "${params.results}/gwas/files/", mode: 'copy'

    input:
    path(bed)
    path(bim)
    path(fam)
    path(covar)

    output:
    path "*.glm.*", emit: gwas
    path "*.log", emit: log
   
    shell:
    '''
    plink2 --bfile !{bed.baseName} \
      -covar !{covar} \
      --ci 0.95 \
      --glm hide-covar firth-fallback\
      --output-chr 26 \
      --adjust \
      --out gwas
    plink2 --bfile !{bed.baseName} \
      --ci 0.95 \
      --glm hide-covar firth-fallback\
      --adjust \
      --output-chr 26 \
      --out gwas_nocovars
    '''
}

process plot {
    label 'small'
    publishDir "${params.results}/gwas/figures/", mode: 'copy'
    
    input:
    tuple val(id), path(gwas), path(log)
	
    output:
    path "*_qqplot.png", emit: qqplot
    path "*_manhattan.png", emit: manhattan
 
    shell:
    '''
    # Identify lambda value calculated by plink2
    awk '/lambda/  {print $11+0}' !{log} > lambda.txt
    # NA p values sometimes happen
    awk '!/'NA'/' !{gwas} > no_na.gwas
    # Remove hash from beginning of line
    sed 's/#//' no_na.gwas | sed 's/LOG(OR)_SE/SE/' | sed 's/CHROM/CHR/'|  sed 's/POS/BP/' | sed 's/ID/SNP/' > valid.gwas
    # Run plots 
    qqplot.R valid.gwas lambda.txt
    manhattan.R valid.gwas
    cp qqplot.png !{id}_qqplot.png 
    cp manhattan.png !{id}_manhattan.png
    '''
}
