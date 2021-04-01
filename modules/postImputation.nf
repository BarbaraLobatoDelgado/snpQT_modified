// STEP E1: Convert vcf to binary plink files and filter all poorly imputed variants based on info score

process filter_imp {
    input:
    path(imp)

    output:
    path "E1.bed", emit: bed
    path "E1.bim", emit: bim
    path "E1.fam", emit: fam
    path "E1.log", emit: log
    
    shell:
    '''
    plink2 --vcf !{imp} \
        --extract-if-info INFO '>'= !{params.info} \
        --make-bed \
        --out E1
    '''
}

// STEP E2: Filter based on MAF 

process filter_maf {
    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "E2.bed", emit: bed
    path "E2.bim", emit: bim
    path "E2.fam", emit: fam
    path "E2.log", emit: log
    
    shell:
    '''
    plink2 --bfile !{bed.baseName} \
        --maf !{params.impute_maf} \
        --make-bed \
        --out E2
    '''
}

// STEP E3: Identify and remove exact duplicated variants

process duplicates_cat1 {
    input:
    path(bed)
    path(bim)
    path(fam)

     output:
    path "E3.bed", emit: bed
    path "E3.bim", emit: bim 
    path "E3.fam", emit: fam
    path "E3.log", emit: log
    
    shell:
    '''
    # Annotate all variants to this format chr:pos:ref:alt and remove exact duplicates
    plink2 --bfile !{bed.baseName} \
        --rm-dup force-first list \
        --make-bed \
        --out E3
    '''
}

// STEP E4: Identify and remove multi-allelics

process duplicates_cat2 {
    input:
    path(bed)
    path(bim)
    path(fam)
    
    output:
    path "E4.bed", emit: bed
    path "E4.bim", emit: bim 
    path "E4.fam", emit: fam
    path "E4.log", emit: log
  
	shell:
    '''
    # Identify the multi-allelics based on position and reference allele
    cut -f 1,4,6 !{bim} | sort | uniq -d | cut -f 2 | grep -w -F -f - !{bim} | cut -f 2 > multi_allelics.txt
    if [[ $(wc -l < multi_allelics.txt) -gt 0 ]]
    then
	plink2 --bfile !{bed.baseName} \
        --exclude multi_allelics.txt \
        --make-bed \
        --out E4
	else
      plink -bfile !{bim.baseName} \
        --make-bed \
        --out E4
    fi
    '''
}

// STEP E5: Identify and remove merged variants

process duplicates_cat3 {
    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "E5.bed", emit: bed
    path "E5.bim", emit: bim
    path "E5.fam", emit: fam
    path "E5.log", emit: log 
    
    shell:
    '''
    cut -f 2 !{bim} | sort | uniq -d > merged_variants.txt

    if [[ $(wc -l < merged_variants.txt) -gt 0 ]]
    then
      plink2 --bfile !{bim.baseName} \
        -extract merged_variants.txt \
        --make-bed \
        --out merged_snps
      plink2 --bfile !{bim.baseName} \
        --exclude merged_variants.txt \
        --make-bed \
        --out excluded_snps
      plink2 --bfile merged_snps \
        --set-all-var-ids @:#:\\$r:\\$a \
        --new-id-max-allele-len 300 missing\
        --make-bed \
        --out annotated
      plink --bfile excluded_snps \
        --bmerge annotated \
        --make-bed \
        --out E5   
    else
      plink -bfile !{bim.baseName} \
        --make-bed \
        --out E5
    fi
    '''
}

// STEP E6: update ids information

process update_ids {
    publishDir "${params.results}/post_imputation/bfiles", mode: 'copy'
    
    input:
    path(bed)
    path(bim)
    path(fam)
    path(user_fam)

    output:
    path "E6.bed", emit: bed
    path "E6.bim", emit: bim
    path "E6.fam", emit: fam
    path "E6.log", emit: log 
    
    shell:
    '''
	awk '{print $1, $2}' !{fam} > old_pids.txt
	awk '{print $1, $2}' !{user_fam} > new_pids.txt
	paste old_pids.txt new_pids.txt > update_ids.txt
	
    plink2 --bfile !{bed.baseName} \
        --update-ids update_ids.txt \
        --make-bed \
        --out E6
    '''
}

// STEP E7: update phenotype information

process update_phenotype {
    publishDir "${params.results}/post_imputation/bfiles", mode: 'copy'
    
    input:
    path(bed)
    path(bim)
    path(fam)
    path(user_fam)

    output:
    path "E7.bed", emit: bed
    path "E7.bim", emit: bim
    path "E7.fam", emit: fam
    path "E7.log", emit: log 
    
    shell:
    '''
    plink2 --bfile !{bed.baseName} \
        --fam !{user_fam} \
        --make-bed \
        --out E7
    '''
}