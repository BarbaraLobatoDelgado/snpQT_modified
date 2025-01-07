// STEP 1: unzip files using password
process unzip_chr_files {
    
    input:
    path(topmed_results_dir)

    output:
    path "*.dose", emit: dose
    path "*.info", emit: info

    script:
    """
    echo "topmed_results_dir: ${topmed_results_dir}"
    echo "Checking files in ${topmed_results_dir}:"
    ls ${topmed_results_dir}
    unzip -P "${params.topmed_password}" ${topmed_results_dir}/chr_*.zip 
    """
}