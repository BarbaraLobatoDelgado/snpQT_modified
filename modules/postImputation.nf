// STEP H1: Convert vcf to binary plink files and filter all poorly imputed variants based on info score
process filter_imp {
    label 'plink2'
	
    input:
    path(imp)

    output:
    path "H1.bed", emit: bed
    path "H1.bim", emit: bim
    path "H1.fam", emit: fam
    path "H1.log", emit: log
    
    shell:
    '''
    plink2 --vcf !{imp} \
        --extract-if-info INFO '>'= !{params.info} \
        --make-bed \
        --out H1
    '''
}

// STEP H2: Filter based on MAF 
process filter_maf {
    label 'plink2'
	
    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "H2.bed", emit: bed
    path "H2.bim", emit: bim
    path "H2.fam", emit: fam
    path "H2.log", emit: log
    
    shell:
    '''
    plink2 --bfile !{bed.baseName} \
        --maf !{params.impute_maf} \
        --make-bed \
        --out H2
    '''
}

// STEP H3: Identify and remove duplicated variants
process duplicates_cat1 {
    label 'plink2' 
 
    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "H3.bed", emit: bed
    path "H3.bim", emit: bim 
    path "H3.fam", emit: fam
    path "H3.log", emit: log
    
    shell:
    '''
    # All variants are annotated to this format chr:pos:ref:alt and remove exact duplicates
    plink2 --bfile !{bed.baseName} \
        --rm-dup force-first list \
        --make-bed \
        --out H3
    '''
}

// STEP H4: Identify and remove multi-allelics
process duplicates_cat2 {
    label 'plink2'
	
    input:
    path(bed)
    path(bim)
    path(fam)
    
    output:
    path "H4.bed", emit: bed
    path "H4.bim", emit: bim 
    path "H4.fam", emit: fam
    path "H4.log", emit: log
  
    shell:
    '''
    # Identify the multi-allelics based on position and reference allele
    cut -f 1,4,6 !{bim} | sort | uniq -d | cut -f 2 | grep -w -F -f - !{bim} | cut -f 2 > multi_allelics.txt
    if [[ $(wc -l < multi_allelics.txt) -gt 0 ]]
    then
        plink2 --bfile !{bed.baseName} \
	    --exclude multi_allelics.txt \
	    --make-bed \
	    --out H4
        else
        plink2 --bfile !{bim.baseName} \
	    --make-bed \
	    --out H4
    fi
    '''
}

// STEP H5: Identify and remove merged variants
process duplicates_cat3 {
    label 'plink1'
 
    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "H5a.bed", emit: bed
    path "H5a.bim", emit: bim
    path "H5a.fam", emit: fam
    path "H5a.log", emit: log 
    // path "merge" ...
    
    shell:
    '''
    # Get rsIDs, sort them and save repeated SNPs
    cut -f 2 !{bim} | sort | uniq -d > merged_variants.txt

    # Count number of lines in merged_variants.txt
    # If there are any lines (i.e. duplicated variants), then run plink commands
    # Note: there are variants with no rsID, noted as ".". These are found as duplicates in many analyses. Not certain if they should be removed.
    if [[ $(wc -l < merged_variants.txt) -gt 0 ]]
    then
      # Create plink files with duplicated rsID variants (most have rsID of .)
      plink --bfile !{bim.baseName} \
          --extract merged_variants.txt \
          --make-bed \
          --out merged_snps
      # Create plink files with all SNPs after excluding those with same rsID
      plink --bfile !{bim.baseName} \
          --exclude merged_variants.txt \
          --make-bed \
          --out excluded_snps
      # Create plink files with duplicated rsID variants, changing their IDs to be ??
      plink --bfile merged_snps \
          --set-all-var-ids @:#:\\$r:\\$a \
          --new-id-max-allele-len 300 missing\
          --make-bed \
          --out annotated
      plink --bfile excluded_snps \
          --bmerge annotated \
          --make-bed \
          --out H5a
    # If there are no duplicated SNPs, generate plink files with different name
    else
      plink --bfile !{bim.baseName} \
          --make-bed \
          --out H5a
    fi
    '''
}

// STEP H5: Identify and remove merged variants
process duplicates_cat4 {
    label 'plink2'
 
    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "H5b.bed", emit: bed
    path "H5b.bim", emit: bim
    path "H5b.fam", emit: fam
    path "H5b.log", emit: log 
    
    shell:
    '''
    cut -f 2 !{bim} | sort | uniq -d > merged_variants.txt

    if [[ $(wc -l < merged_variants.txt) -gt 0 ]]
    then
      plink2 --bfile merged_snps \
          --set-all-var-ids @:#:\\$r:\\$a \
          --new-id-max-allele-len 300 \
          --make-bed \
          --out annotated
      plink2 --bfile excluded_snps \
          --bmerge annotated \
          --make-bed \
          --out H5   
    else
      plink2 -bfile !{bim.baseName} \
          --make-bed \
          --out H5b
    fi
    '''
}

// STEP H6: update ids information
process update_ids {    
    label 'plink2'
	
    input:
    path(bed)
    path(bim)
    path(fam)
    path(user_fam)

    output:
    path "H6.bed", emit: bed
    path "H6.bim", emit: bim
    path "H6.fam", emit: fam
    path "H6.log", emit: log 
    
    shell:
    '''
    awk '{print $1, $2}' !{fam} > old_pids.txt
    awk '{print $1, $2}' !{user_fam} > new_pids.txt
    paste old_pids.txt new_pids.txt > update_ids.txt
	
    plink2 --bfile !{bed.baseName} \
        --update-ids update_ids.txt \
        --make-bed \
        --out H6
    '''
}

// STEP H7: update phenotype information
process update_phenotype {
    label 'plink2'
    publishDir "${params.results}/post_imputation/bfiles", mode: 'copy'
    
    input:
    path(bed)
    path(bim)
    path(fam)
    path(user_fam)

    output:
    path "H7.bed", emit: bed
    path "H7.bim", emit: bim
    path "H7.fam", emit: fam
    path "H7.log", emit: log 
    
    shell:
    '''
    plink2 --bfile !{bed.baseName} \
        --fam !{user_fam} \
        --make-bed \
        --out H7
    '''
}