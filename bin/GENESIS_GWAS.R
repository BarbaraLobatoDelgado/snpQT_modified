#!/usr/bin/env Rscript

# GWAS using GENESIS package (mixed models)

# Sources
# GWAS_Dosage.R (author: Gerard Temprano-Sagrera)
# Tutorial
# https://www.bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/assoc_test.html


# Load libraries
pacman::p_load(
  an9elproject,
  data.table,
  tidyverse, 
  magrittr,
  ggplot2,
  GWASTools,
  SNPRelate,
  GENESIS, 
  qqman,
  biomaRt,
  install = FALSE, update = FALSE
)

# When running this script from snpQT pipeline
# Args
# 1: ?
# 2: ?
# ...
# args <- commandArgs(trailingOnly = TRUE)

# Set to use half the cores
num_cores <- parallel::detectCores() - 4
# you can then use function mclapply() to run several operations in different cores


# ------------------------------------------------- # 
#         Load genotype data from PLINK files
# ------------------------------------------------- # 

use_pruned_snps <- FALSE

# Load PLINK files
plink_files_dir_path <- "/mnt/ir-bioinf02/home/blobato/oncothromb02/GWAS/snpQT_results/results/post_imputation/bfiles/"

if (use_pruned_snps) {
  # Use SNPs resulting from pruning
  plink_files_name <- "pruned_output"
} else {
  # Using all SNPs
  plink_files_name <- "oncoth2_eligible_patients" # "H7"
}
 

# ------------------------------------------------- # 
#         Manage genotype data as GDS
# ------------------------------------------------- # 

# GDS file path
gds_file <- paste0(plink_files_dir_path, "genotype.gds")

# Check if GDS file already exists and delete it
if (file.exists(gds_file)) {
  unlink(gds_file)
}

# Convert PLINK files to GDS to handle genotype data
snpgdsBED2GDS(
  bed.fn = paste0(plink_files_dir_path, plink_files_name, ".bed"), # "./GWAS/imputation_results/oncoth2.all_chr.TOPMED.hg38.bed",
  fam.fn = paste0(plink_files_dir_path, plink_files_name, ".fam"),
  bim.fn = paste0(plink_files_dir_path, plink_files_name, ".bim"),
  out.gdsfn = gds_file
)

# Load GDS (Genomic Data Structures) file
genotype_file <- GdsGenotypeReader(gds_file)
genotype_data <- GenotypeData(genotype_file)

# Get IDs 
genotype_ids <- getScanID(genotype_data)
# Get chromosomes
genotype_chrs <- table(getChromosome(genotype_data))


# Show some information from GDS file
print("IIDs in GDS file:")
print(genotype_ids)
print("Chromosomes loaded into GenotypeData:")
print(unique(getChromosome(genotype_data)))
print("Number of SNPs in GDS file:")
print(nsnp(genotype_data))
print("SNPs rsIDs in GDS file:")
print(getSnpID(genotype_data))

# tt = PCs_pcair_matrix %>% filter(rownames(.) %in% oncoth2_pheno$patient_code)
# 
# 
# genotype_data_subset <- subsetGenotypeData(genotype_data, sample.id = oncoth2_pheno$patient_code)

# ------------------------------------------------- # 
#              Load phenotype data
# ------------------------------------------------- # 

# Load an9elproject object 
oncoth2 <- get_project("oncoth2") # , version = "0.0.8008"
# oncoth2_ids <- oncoth2$data %>%
#   select(id, patient_code)

table(oncoth2$data$eligible, useNA = "ifany")
table(oncoth2$data$informed_consent, useNA = "ifany")
table(oncoth2$data$VTE, useNA = "ifany")


# Select dependent variable and covariates
# oncoth2_pheno <- oncoth2$data %>%
#   # Make sure to leave out those with 
#   filter(informed_consent == "Yes") %>%
#   filter(eligible == "Yes") %>%
#   # filter(patient_code %in% pc$IID) %>%
#   select(id, patient_code, VTE, age_cancer_dx, sex)


# --- Remove when patients are excluded beforehand --- #
oncoth2_pheno <- oncoth2$data %>%
  # Make sure to leave out those with
  filter(informed_consent == "Yes") %>%
  filter(eligible == "Yes") %>%
  filter(!is.na(VTE)) %>%
  select(id, patient_code, VTE, age_cancer_dx, sex) 

# ------------------------------------------------- # 
#                     Do PCA
# ------------------------------------------------- # 

# We think reasonable to use the first 4 PCs as covariates to correct for
# unknown effects caused by genetic variability, as well as ancestry (some individuals
# don't have European ancestry but we are going to try to keep them in the analysis) 

# Do PCA on genome-wide SNP data
# The principal components account for population structure in a sample
pca <- pcair(
  gdsobj = genotype_data, # geno_iterator, 
  num.cores = num_cores
)
# Get principal components
# Matrix
PCs_pcair_matrix <- pca$vectors

# Dataframe
PCs_pcair_df <- as.data.frame(pca$vectors)
# Rename columns
colnames(PCs_pcair_df) <- c(paste0("PC", 1:dim(PCs_pcair_df)[2]))

# Proportion of variance explained by each PC
eigenvalues <- pca$values


variance_explained_by_PC <- function(eigenvalues) {
  #'
  #'
  
  # Calculate proportion of variance explain
  prop_variance <- eigenvalues / sum(eigenvalues)
  
  # Calculate cumulative proportion of variance explained
  cumulative_variance <- cumsum(prop_variance)
  
  # Create dataframe for plotting
  variance_df <- data.frame(
    PC = 1:length(prop_variance),
    Proportion = prop_variance,
    Cumulative = cumulative_variance
  )
  
  # Create scree plot
  ggplot(variance_df, aes(x = PC, y = Proportion)) + 
    # Barplot for proportion of variance explained
    geom_bar(stat = "identity", fill = "skyblue") +
    geom_text(aes(label = sprintf("%.2f", Proportion)), vjust = -0.5, size = 1.5) +
    # Line for cumulative proportion of variance explained
    # geom_line(aes(y = Cumulative), color = "red", size = 1) + 
    # geom_point(aes(y = Cumulative), color = "red", size = 2) +
    labs(
      title = "Scree Plot", 
      x = "Principal Component", 
      y = "Proportion of Variance Explained"
    ) + 
    theme_minimal() + 
    scale_x_continuous(breaks = 1:length(prop_variance)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 14)
    )
 
}

# Show scree plot
variance_explained_by_PC(eigenvalues = eigenvalues)

# Select subset of patients
oncoth2_pheno %<>%
  filter(patient_code %in% getScanID(genotype_data))


# ------------------------------------------------- # 
#           Create ScanAnnotationDataFrame
# ------------------------------------------------- # 

# Create ScanAnnotationDataFrame 
# Contains outcome and covariate data, as well as unique identifier (scanID)
# It can be paired with genotype
scan_data <- data.frame(
  "scanID" = oncoth2_pheno$patient_code, 
  "age" = oncoth2_pheno$age_cancer_dx, 
  "sex" = as.factor(ifelse(oncoth2_pheno$sex == "Man", "M", "F")), # expected format for ScanAnnotationDataFrame
  "pc1" = PCs_pcair_df$PC1, 
  "pc2" = PCs_pcair_df$PC2,
  "pc3" = PCs_pcair_df$PC3,
  "pc4" = PCs_pcair_df$PC4, 
  "pc5" = PCs_pcair_df$PC5,
  "pc6" = PCs_pcair_df$PC6,
  "pc7" = PCs_pcair_df$PC7,
  "pc8" = PCs_pcair_df$PC8,
  "pc9" = PCs_pcair_df$PC9,
  "pc10" = PCs_pcair_df$PC10,
  "pheno" = ifelse(oncoth2_pheno$VTE == "Yes", 1, 0) # expected format for assocTestSingle
)

# Create ScanAnnotationDataFrame object
scanAnnot <- ScanAnnotationDataFrame(data = scan_data)


# ------------------------------------------------- # 
#           Estimate kinship coefficients
# ------------------------------------------------- # 

# Create GenotypeBlockIterator class to iterate over blocks of 10,000 SNPs
geno_iterator <- GenotypeBlockIterator(genotype_data, snpBlock = 10000)

# Model-free estimation of recent genetic relatedness
pcrel_result <- pcrelate(
  gdsobj = geno_iterator, 
  pcs = PCs_pcair_matrix
)

# ------------------------------------------------- # 
#         Create Genetic Relationship Matrix
# ------------------------------------------------- # 

# Genetic relationship matrix
grm <- pcrelateToMatrix(pcrel_result)
# grm[1:5, 1:5]


# ------------------------------------------------- # 
#                   Fit models
# ------------------------------------------------- # 

# Fit null model
# Mixed model that contains all covariates but not genotype data
# Null hypothesis is that all SNPs have no effect on outcome
nullmodel <- fitNullModel(
  scanAnnot, 
  outcome = "pheno",
  # Use multiple fixed effect covariates
  covars = c("age", "sex", "pc1", "pc2", "pc3", "pc4", "pc5", 
             "pc6", "pc7","pc8","pc9", "pc10"),
  # Use GRM to account for relatedness
  cov.mat = grm,
  family = "binomial"
  )

# Variance components estimates for the random effects
nullmodel$varComp

# Fixed effects covariates info
fixed_effects <- nullmodel$fixef

# 
nullmodel$fit

# Heritability of the trait
varCompCI(null.model = nullmodel, prop = TRUE)

# Run SNP-phenotype association tests
assoc_results_table <- assocTestSingle(
  gdsobj = geno_iterator, # establish how many SNPs are read at a time
  null.model = nullmodel, 
  test = "Score.SPA" # Saddle point approximation works best for imbalance case to control ratio
  # BPPARAM = BiocParallel::MulticoreParam(workers = num_cores) # Parallelise
)

# Close GenotypeBlockIterator connection to avoid assocTestSingle() to show
# only last chromosome
close(geno_iterator)

# ------------------------------------------------- # 
#               Results interpretation
# ------------------------------------------------- # 

# Annotate genes 
# Connect to the Emsembl database
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Annotate variants with gene information
gene_annot <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", 
                 "ensembl_gene_id", "external_gene_name"), 
  filters = c("chromosome_name", "start", "end"), 
  values = list(assoc_results_table$chr, 
                assoc_results_table$pos, 
                assoc_results_table$pos), 
  mart = ensembl
)

# Merge annotations with GWAS results
assoc_results_table_annot <- merge(
  assoc_results_table, gene_annot, 
  by.x = c("chr", "pos"), 
  by.y = c("chromosome_name", "start_position")
  )


# Convert chr from character to numeric
assoc_results_table_reformated <- assoc_results_table %>%
  mutate(chr = case_when(
    chr == "X" ~ "23",
    TRUE ~ chr
  )) %>%
  mutate(chr = as.numeric(chr)) %>%
  # Calculate -log10() of raw p-value
  mutate(raw_pval_log = -log10(Score.pval)) %>%
  # Apply Bonferroni correction to p-values
  mutate(bonferroni_adjusted_pvals = p.adjust(
    p = assoc_results_table$Score.pval, 
    method = "bonferroni"
  )) %>%
  # Apply FDR correction to p-values
  mutate(fdr_adjusted_pvals = p.adjust(
    p = assoc_results_table$Score.pval, 
    method = "fdr"
  ))

# Save top table results
write.csv(assoc_results_table_reformated, 
          "./GWAS/snpQT_results/results/post_imputation/GWAS_results.csv", 
          row.names = FALSE
          )

# 
# ggplot(assoc_results_table_reformated, aes(x = adjust_pvalues)) +
#   geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) +
#   labs(x = "Adjusted p-value", y = "Frequency", title = "Histogram of Adjusted p-values") +
#   theme_minimal()


calculate_gif <- function(pvals) {
  #' Calculate Genetic Inflation Factor (lambda in QQ plots)
  #'
  
  # Compute observed -log10(p-values)
  observed_logp <- -log10(pvals)
  
  # Compute expected -log10(p-values) under the null hypothesis
  # Null hypothesis: p-values follow a uniform distribution
  expected_logp <- -log10(ppoints(length(pvals)))
  
  # Calculate lambda
  lambda <- median(observed_logp) / median(expected_logp)
  
  return(lambda)
}



# QQ plots
# Raw p-values
lambda_raw_pvals <- calculate_gif(assoc_results_table_reformated$Score.pval)
qq(pvector = assoc_results_table_reformated$Score.pval, 
   main = "QQ plot for association with raw p-values", 
   sub = paste("Lambda = ", lambda_raw_pvals)
   )

# # Bonferroni adjusted p-values
# lambda_bonfe_pvals <- calculate_gif(assoc_results_table_reformated$bonferroni_adjusted_pvals)
# qq(pvector = assoc_results_table_reformated$bonferroni_adjusted_pvals,
#    main = "QQ plot for association with adjusted p-values using Bonferroni method",
#    sub = paste("Lambda = ", lambda_bonfe_pvals)
#    )

# FDR adjusted p-values
lambda_fdr_pvals <- calculate_gif(assoc_results_table_reformated$fdr_adjusted_pvals)
qq(pvector = assoc_results_table_reformated$fdr_adjusted_pvals,
   main = "QQ plot for association with adjusted p-values using FDR method", 
   sub = paste("Lambda = ", lambda_fdr_pvals)
)

# Manhattan plot
# Raw p-values
manhattan(
  assoc_results_table_reformated,
  chr = "chr",
  bp = "pos",
  p = "Score.pval",
  snp = "variant.id",
  main = "Manhattan Plot with raw p-values",
  ylim = c(0, -log10(min(assoc_results_table_reformated$P)) + 1), # Adjust for P-value range
  col = c("blue4", "orange3"), 
  suggestiveline = -log10(1e-05),
  genomewideline = -log10(5e-08), 
  annotateTop = TRUE
)


# # Bonferroni adjusted p-values
# manhattan(
#   assoc_results_table_reformated,
#   chr = "chr",
#   bp = "pos",
#   p = "bonferroni_adjusted_pvals",
#   snp = "variant.id",
#   main = "Manhattan Plot with Bonferroni adjusted p-values",
#   ylim = c(0, -log10(min(assoc_results_table_reformated$P)) + 1), # Adjust for P-value range
#   col = c("blue4", "orange3"), 
#   suggestiveline = -log10(1e-05),
#   genomewideline = -log10(5e-08), 
#   annotateTop = TRUE
# )

# FDR adjusted p-values
manhattan(
  assoc_results_table_reformated,
  chr = "chr",
  bp = "pos",
  p = "fdr_adjusted_pvals",
  snp = "variant.id",
  main = "Manhattan Plot with FDR adjusted p-values",
  ylim = c(0, -log10(min(assoc_results_table_reformated$P)) + 1), # Adjust for P-value range
  col = c("blue4", "orange3"), 
  suggestiveline = -log10(1e-05),
  genomewideline = -log10(5e-08)
)




