// STEP 1: unzip files using password
process unzip_chr_files {
    
    input:
    path(topmed_results_dir)

    output:
    path "*.dose", emit: dose
    path "*.info", emit: info

    script:
    """
    unzip ${topmed_results_dir}/chr_*.zip -P ${params.topmed_password}
    """
}