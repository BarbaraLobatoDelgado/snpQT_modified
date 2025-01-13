// STEP 1: unzip files using password
// There are two types of files: dose.vcf.gz and info.gz
process unzip_vcfs {
    
    input:
    path(topmed_results_dir)

    output:
    path "*.dose.vcf.gz", emit: vcf
    path "*.info.gz", emit: info

    script:
    """
    for zip_file in ${topmed_results_dir}/chr_*.zip; do
        echo "Unzipping: \$zip_file"
        unzip -P '${params.topmed_password}' \$zip_file
    done
    """
}

// STEP 2: create indices files
process create_index {
    label 'bcftools'

    input:
    path(vcf)

    output:
    path "chr*.dose.vcf.gz.csi", emit: csi

    script:
    """
    # Create index
    for file in ${vcf.join(' ')}; do
        if [[ ! -f "\${file}.tbi" ]]; then
            echo "Indexing \${file}..."
            bcftools index \${file}
        fi
    done
    """
}


// STEP 3: join all VCFs into one
process concat_vcfs {
    label 'bcftools1_5'

    input:
    path(csi)
    path(vcf)

    output:
    path "concat.vcf.gz", emit: concat_vcf

    script:
    """
    # Create list with VCFs in the order they should be concatenated
    echo ${vcf} | tr ' ' '\\n' | sort -V > list_sorted_vcfs.txt
    
    # Merge the sorted VCF files using bcftools merge
    bcftools concat --file-list list_sorted_vcfs.txt -o concat.vcf.gz -O z 
    """
}

// STEP 4: sort VCF files by chromosome and position
process sort_vcf {
    label 'bcftools'

    input:
    path(vcf)

    output:
    path "sorted.vcf.gz", emit: sorted_vcf 

    script:
    """
    echo bcftools --version
    # Sort VCF file by chromosome and position
    bcftools sort ${vcf} -O z -o sorted.vcf.gz
    """
}

// STEP 6: remove "0_" substring from IDs (added during TOPMed imputation)
process trim_ids {
    label 'bcftools1_5'

    input:
    path(vcf)

    output:
    path "correct_ids.vcf.gz", emit: correct_ids_vcf 

    shell:
    '''
    # Get old IDs
    bcftools query -l !{vcf} > old_ids.txt
    # Remove substring
    sed 's/^0_//g' old_ids.txt > new_ids.txt
    # Create file with old and new IDs, each in a column
    paste old_ids.txt new_ids.txt > ids_mapping.txt
    # Change IDs in VCF
    bcftools reheader -s ids_mapping.txt -o correct_ids.vcf.gz !{vcf}
    '''
}

// STEP 5: Convert vcf to binary plink files and filter all poorly imputed variants based on R2 score
process filter_imp_quality_topmed {
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
        --extract-if-info R2 '>'= !{params.r2} \
        --make-bed \
        --out H1
    '''
}

