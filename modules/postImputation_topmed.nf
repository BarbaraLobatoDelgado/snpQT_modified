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

// // STEP 2: gunzip VCF files
// process gunzip_vcf_files {

//     input:
//     path(compressed_vcf)

//     output:
//     path "*.dose.vcf", emit: vcf

//     script:
//     """
//     echo ${compressed_vcf}
//     for gz_file in ${compressed_vcf}; do
//         echo "Gunzipping: \$gz_file" 
//         gunzip -f \$gz_file # -f is necessary for the command to work with symlinks, otherwise it throws an error
//     done
//     """
// }

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

// STEP 6: remove "0_" substring from IDs (added during TOPMed imputation)
process trim_ids {
    label 'plink1'

    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    // path "H2.bed", emit: bed
    // path "H2.bim", emit: bim
    path "H2.fam", emit: fam
    // path "H2.log", emit: log

    shell:
    '''
    # plink --bfile !{bim.baseName} --freq --out H2
    awk '{print $1, $2, $3, $4, $5, $6}' !{fam} | sed 's/^0_//g' > H2.fam
    '''
}

// // STEP 2b: convert files to PLINK format and join them all
// process convert_vcfs_to_plink {
//     label 'plink1'

//     input:
//     path(vcf)

//     output:
//     path "chr*.bim", emit: bim 
//     path "chr*.bed", emit: bed
//     path "chr*.fam", emit: fam

//     script:
//     """
//     echo ${vcf}
//     for vcf_file in ${vcf}; do
//         # Get VCF base name
//         # base_name=\$(basename \$vcf_file .dose.vcf)
//         plink --vcf \$vcf_file --make-bed --out \$(basename \$vcf_file .dose.vcf.gz) --const-fid
//     done
//     """
// }

// // STEP 3: merge all chromosomes into one plink file
// process merge_plink_files {
//     label 'plink1'

//     input:
//     path(bim)
//     path(bed)
//     path(fam)

//     output:
//     path "merged_dataset.bim", emit: bim
//     path "merged_dataset.bed", emit: bed
//     path "merged_dataset.fam", emit: fam

//     script:
//     """
//     echo "Creating list of files to be merged..."
//     # Remove parentheses and blank spaces, and add new line after comma
//     sorted_files=\$(echo ${bim.baseName} | tr -d '[:space:][]' | tr ',' '\n' | sort -V) 
//     echo "\$sorted_files" >> list_sorted_files.txt

//     # Get first line from list (normally chr1)
//     echo "Identify first file to be merged..."
//     first_file=\$(head -n 1 list_sorted_files.txt)
//     echo "First file is: \$first_file"

//     # Create new file with the rest of files to be merged
//     echo "Remove first line from list"
//     tail -n +2 list_sorted_files.txt >> list_sorted_files_except_first.txt

//     # Unify all in PLINK format
//     echo "Merging all PLINK files into one..."
//     plink --bfile \$first_file --merge-list list_sorted_files_except_first.txt --make-bed --out merged_dataset


//     # Change IDs and add info from oncoth2.fam in ./
//     """

// }