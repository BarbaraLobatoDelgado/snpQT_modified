// Pre-imputation 
// =============================================================================

// STEP F1: Set chromosome codes ---------------------------------------------
process set_chrom_code {
    label 'plink2'
	
    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "F1.bed", emit: bed
    path "F1.bim", emit: bim
    path "F1.fam", emit: fam
    path "F1.log", emit: log
    
    shell:
    """
    plink2 --bfile !{bed.baseName} \
        --output-chr ${params.output_chr} \
	    --make-bed \
        --out F1  
    """
}

// STEP F2: D2: Remove ambiguous SNPs and flip reverse SNPs -------------------------
// note: taken care of by popstrat modules D4 now

// STEP F2a: 

// STEP F3: Remove one of each pair of duplicated SNPs 
process fix_duplicates {
    label 'plink2'
	
    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "F3.bed", emit: bed
    path "F3.bim", emit: bim
    path "F3.fam", emit: fam
    path "F3.log", emit: log
    
    shell:
    """
    # F3: Deduplicate variants
    plink2 --bfile !{bed.baseName} \
        --output-chr ${params.output_chr} \
        --rm-dup 'force-first' \
        --make-bed \
        --out F3
    """
}

// STEP F4: Convert Plink file to .bcf file ---------------------------------
process to_vcf {
    label 'plink2'
	
    input:
    path(bed)
    path(bim)
    path(fam)
    
    output:
    path "F4.vcf.gz", emit: vcf
    path "F4.log", emit: log

    shell:
    '''
    plink2 --bfile !{bed.baseName} \
        --export vcf bgz \
        --out F4
    '''
}

// STEP F5: Convert Plink file to .bcf file ---------------------------------
process to_bcf {
    label 'bcftools'
	
    input:
    path(vcf)
    
    output:
    path "F5.bcf", emit: bcf

    shell:
    '''
    bcftools index !{vcf}
    bcftools convert !{vcf} -Ou -o F5.bcf --threads !{task.cpus}
    '''
}

// STEP F6: Check and fix the REF allele --------------------------------------
// Checks and corrects reference alleles in the provided BCF/VCF file by comparing them to a reference genome
// When a mismatch is found, the tool attempts to correct it based on the reference provided
// Using a dbSNP file can help validate that the reference alleles are accurate, as it provides known allele information for each SNP position
process check_ref_allele {
    label 'bcftools'
	
    input:
    path(bcf)
    path(g37)
    path(g38)
    path(dbsnp37)
    path(dbsnp38)
    path(dbsnp_idx38)
    path(dbsnp_idx37)

    output:
    path "F6.bcf", emit: bcf

    script:
    // Define variables based on conditionals
    def genome_reference
    def dbsnp
    def dbsnp_idx

    if (params.input_build == 38) {
        genome_reference = g38
        dbsnp = dbsnp38
        dbsnp_idx = dbsnp_idx38
        // println "Using ${g38} as genome reference file and ${dbsnp38} as SNP database"
    } else {
        genome_reference = g37
        dbsnp = dbsnp37
        dbsnp_idx = dbsnp_idx37
        // println "Using ${g37} as genome reference file and ${dbsnp37} as SNP database"
    }   

    shell:
    """
    # Swap the alleles
    bcftools +fixref ${bcf} \
        -Ob -o F6.bcf -- \
        -d -f ${genome_reference} \
        -i ${dbsnp}
    """
}

// STEP F7: Sort BCF, convert .bcf file to .vcf.gz file and index the vcf.gz -------------------------------
process bcf_to_vcf {
    label 'bcftools'
    publishDir "${params.results}/preImputation/files", mode: 'copy'

    input:
    path(bcf)
    
    output:
    path "F7.vcf.gz", emit: vcf
    path "F7.vcf.gz.csi", emit: idx
	
    shell:
    '''
    bcftools sort !{bcf} | bcftools convert -Oz --threads !{task.cpus} > F7.vcf.gz
    bcftools index F7.vcf.gz 
    '''
}

// Imputation 
// =============================================================================

// STEP G1: Split vcf.gz file in chromosomes and index all chroms ---------------------------------
process split_user_chrom {
    label 'bcftools'
	
    input:
    path(vcf)
    path(idx)
    each chr
    
    output:
    tuple val(chr), file('G1.vcf.gz'), file('G1.vcf.gz.csi'), emit: chrom 

    shell:
    '''
    bcftools view -r !{chr} !{vcf} -Oz -o G1.vcf.gz --threads !{task.cpus}
    bcftools index G1.vcf.gz --threads !{task.cpus}
    '''
}

// STEP G2: Perform phasing using shapeit4 --------------------------
process phasing {
    label 'shapeit'
	
    input:
    tuple val(chr), file('G1.vcf.gz'), file('G1.vcf.gz.csi'), \
        file('genetic_maps.b37.tar.gz')  

    output:
    tuple val(chr), file('G2.vcf.gz'), emit: chrom, optional: true

    shell:
    '''
    gunzip genetic_maps.b37.tar.gz
    tar -xf genetic_maps.b37.tar
    gunzip chr!{chr}.b37.gmap.gz # decompress the chromosome we need 

    # || true allows optional output without an error
    # people might not always be imputing every chromosome
    shapeit4 --input G1.vcf.gz \
        --map chr!{chr}.b37.gmap \
        --region !{chr} \
        --thread !{task.cpus} \
        --output G2.vcf.gz \
        --log log_chr.txt || true     
    '''
}

// STEP G3: Index phased chromosomes --------------------------
process bcftools_index_chr {
    label 'bcftools'
	
    input:
    tuple val(chr), path('chr.vcf.gz')

    output:
    tuple val(chr), path('chr.vcf.gz'), path('chr.vcf.gz.csi'), emit: chrom_idx

    shell:
    '''
    bcftools index chr.vcf.gz 
    '''
}

// STEP G4: Tabix reference files ----------------------------
process tabix_chr {
    label 'small'
    label 'tabix'
	
	
    input:
    tuple val(chr), path('chr.vcf.gz')

    output:
    tuple val(chr), path('chr.vcf.gz'), path('chr.vcf.gz.tbi'), emit: chrom_idx

    shell:
    '''
    tabix -p vcf chr.vcf.gz
    '''
}

// STEP G5: Convert vcf reference genome into a .imp5 format for each chromosome
process convert_imp5 {
    label 'impute5'
	
    input:
    tuple val(chr), file('ref_chr.vcf.gz'), file('ref_chr.vcf.gz.tbi')

    output:
    tuple val(chr), file('1k_b37_reference_chr.imp5'), \
        file('1k_b37_reference_chr.imp5.idx'), emit: chrom

    shell:
    '''
    imp5Converter!{params.impute5_version} --h ref_chr.vcf.gz \
        --r !{chr} \
	--threads !{task.cpus} \
        --o 1k_b37_reference_chr.imp5
    '''
}

// STEP G6: Perform imputation using impute5 ---------------------------------
// join phased vcfs with imp5 based on chrom value 
// then combine so each tuple element has a shapeit4 map file 

process impute5 {
    label 'bigmem'
    label 'impute5'
	
    input:
    tuple val(chr), file('1k_b37_reference_chr.imp5'), \
        file('1k_b37_reference_chr.imp5.idx'), file('G2.vcf.gz'), \
        file('G2.vcf.gz.csi'), file('genetic_maps.b37.tar.gz') 
       
    output:
    path "${chr}.vcf.gz", emit: imputed

    shell:
    '''
    tar -xzf genetic_maps.b37.tar.gz
    gunzip chr!{chr}.b37.gmap.gz # decompress the chromosome we need

    impute5!{params.impute5_version} --h 1k_b37_reference_chr.imp5 \
        --m chr!{chr}.b37.gmap \
        --g G2.vcf.gz \
        --r !{chr} \
        --out-gp-field \
        --o !{chr}.vcf.gz
    '''
}

// STEP G7: Merge all imputed chromosomes with bcftools, so that multi-allelics can be merged, -n is used since files are already sorted after imputation
process merge_imp {
    label 'bcftools'
	
    input:
    path(imp)
    
    output:
    path 'merged_imputed.vcf.gz', emit: vcf
    
    shell:
    '''
    # file order is important so use command substition
    bcftools concat -n $(ls *.vcf.gz | sort -t . -k 1n) -Oz -o merged_imputed.vcf.gz --threads !{task.cpus}
    '''
}


// Steps for imputation in TOPMed server 
// =============================================================================


// STEP F3a: create allele frequency file and output chromosomes as 1-26
process create_afreq_file {
    label 'plink2'

    input:
    path(bed)
    path(bim)
    path(fam)

    output:
    path "F3a.afreq", emit: afreq
    path "F3a.log", emit : log

    shell:
    """
    # Create allele frequency file
    plink2 --bfile !{bed.baseName} \
        --freq \
        --out F3a
    """
}

// STEP F3b: Run McCarthy Group Tools Perl script (v4.3.0)
// Does a series of checks before imputation
// Checks: strand, alleles, position, REF/ALT assignments and frequency differences
// Removes: A/T & G/C SNPs if MAF > 0.4, SNPs with differing alleles, SNPs with > 0.2 allele frequency difference (can be removed/changed in V4.2.2), SNPs not in reference panel 
// Produces: a set of plink commands to update or remove SNPs based on the checks as well as a file (FreqPlot) of cohort allele frequency vs reference panel allele frequency
process run_mccarthy_tools {
    label 'perl'
    label 'plink1'

    // Number of CPUs
    cpus = 16

    input:
    path(bed)
    path(bim)
    path(fam)
    path(afreq)
    path(mccarthy_tools_script)
    path(topmed_ref_hg38)

    output:
    path "F3-updated-chr*.vcf", emit: chr_vcf 

    shell:
    '''
    # Run MacCarthy Group tools Perl script to generate directory with Bash commands and files # -output mccarthy_tools
    perl !{mccarthy_tools_script} \
        -b !{bim} \
        -f !{afreq} \
        -r !{topmed_ref_hg38} \
        -h \
        --output 

    # Go into folder
    # cd mccarthy_tools
    # Give permision for execution
    chmod +x Run-plink.sh
    # Copy the original script to a new file to prevent accumulative modifications
    # cp Run-plink.sh Run-plink_modified.sh
    # Change lines starting with plink to plink1.9
    sed -i 's/^plink/plink1.9/' Run-plink.sh 
    # Remove last line 
    # sed -i '\$d' Run-plink.sh
    # Run plink commands
    sh Run-plink.sh
    '''
}

// STEP F3c: Run .sh file containing PLINK commands and get VCF files separated for chromosome
process prepare_VCFs_for_TOPMed_imputation {
    label 'bcftools'

    // Directory in which output files is saved
    publishDir "${params.results}/VCFs_for_TOPMed_imputation/", mode: 'copy'

    // // Copy file in new directory instead of creating symlink
    // stageInMode 'copy'

    input:
    path(chr_vcf)
    
    output:
    // path "*.vcf.gz", "*vcf.gz.tbi"
    // tuple(path("*.vcf.gz"), path("*.vcf.gz.tbi")), emit: vcf_gz, vcf_gz_tbi
    path "*.vcf.gz", emit: vcf_gz
    path "*.vcf.gz.tbi", emit: vcf_gz_tbi

    shell:
    '''
    echo !{chr_vcf}
    
    # Loop over each VCF file and apply the sed command separately
    for file in !{chr_vcf}; do
        # Reformat chromosome field, such as to have "chr" before number
        sed '/^#/! s/^/chr/' "$file" > "${file%.vcf}_reformat.vcf"
    done

    # Loop over each VCF files
    for file in *_reformat.vcf; do
        # Compress VCF files using bgzip command and create index with tabix
        bgzip "$file" && tabix -p vcf "$file.gz";
    done
    '''
}

// Finished!

