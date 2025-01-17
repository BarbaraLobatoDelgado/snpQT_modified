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

# Load PLINK files
plink_files_dir_path <- "/mnt/ir-bioinf02/home/blobato/oncothromb02/GWAS/snpQT_results/results/post_imputation/bfiles/"
# Using all SNPs
# plink_files_name <- "H7" 
# Use SNPs resulting from pruning
plink_files_name <- "pruned_output"


# ------------------------------------------------- # 
#         Convert genotype data to GDS
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

# Show some information from GDS file
print("IIDs in GDS file:")
print(getScanID(genotype_data))
print("Number of SNPs in GDS file:")
print(nsnp(genotype_data))
print("SNPs rsIDs in GDS file:")
print(getSnpID(genotype_data))


# ------------------------------------------------- # 
#              Load phenotype data
# ------------------------------------------------- # 

# Use ID from an9elproject object
oncoth2 <- get_project("oncoth2") # , version = "0.0.8008"
# oncoth2_ids <- oncoth2$data %>%
#   select(id, patient_code)

table(oncoth2$data$eligible, useNA = "ifany")
table(oncoth2$data$informed_consent, useNA = "ifany")
table(oncoth2$data$VTE, useNA = "ifany")


# We think reasonable to use the first 4 PCs as covariates to correct for
# unknown effects caused by genetic variability, as well as ancestry (some individuals
# don't have European ancestry but we are going to try to keep them in the analysis) 


# Create ScanAnnotationDataFrame 
# Contains outcome and covariate data, as well as unique identifier (scanID)
# It can be paired with genotype



# Select dependent variable and covariates
# oncoth2_pheno <- oncoth2$data %>%
#   # Make sure to leave out those with 
#   filter(informed_consent == "Yes") %>%
#   filter(eligible == "Yes") %>%
#   # filter(patient_code %in% pc$IID) %>%
#   select(id, patient_code, VTE, age_cancer_dx, sex)


# --- Remove when patients are excluded beforehand --- #
oncoth2_pheno <- oncoth2$data %>%
  # # Make sure to leave out those with
  filter(informed_consent == "Yes") %>%
  filter(eligible == "Yes") %>%
  select(id, eligible, VTE, age_cancer_dx, sex) 



# --- Remove when patients are excluded beforehand --- #

# Proportion of VTE
# print(table(oncoth2_pheno$VTE))

print(paste("There are", dim(oncoth2_pheno)[1], "patients in ONCOTHROMB2 with clinical information that can be analysed."))
# print(paste("There are only", dim(pc)[1], "individuals genotyped in ONCOTHROMB2."))

# # Get clinical info for genotyped patients
# oncoth2_pheno_genotyped <- oncoth2_pheno %>%
#   filter(patient_code %in% pc$IID)

# ------------------------------------------------- # 
#                     Do PCA
# ------------------------------------------------- # 

# Create GenotypeBlockIterator class to iterate over blocks of 10,000 SNPs
geno_iterator <- GenotypeBlockIterator(genotype_data)

# Do PCA on genome-wide SNP data
# The principal components account for population structure in a sample
pca <- pcair(
  gdsobj = geno_iterator, 
  num.cores = num_cores
)
# Get principal components
PCs_pcair_matrix <- as.data.frame(pca$vectors)
# Rename columns
colnames(PCs_pcair_matrix) <- c(paste0("PC", 1:dim(PCs_pcair_matrix)[2]))

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
oncoth2_pheno_test <- oncoth2_pheno %>%
  filter(patient_code %in% getScanID(genotype_data))


# ------------------------------------------------- # 
#           Create ScanAnnotationDataFrame
# ------------------------------------------------- # 

# Create dataframe for ScanAnnotationDataFrame
scan_data <- data.frame(
  "scanID" = oncoth2_pheno$patient_code, 
  "age" = oncoth2_pheno$age_cancer_dx, 
  "sex" = ifelse(oncoth2_pheno$sex == 0, "M", "F"), # valid format for ScanAnnotationDataFrame
  "pc1" = PCs_pcair_matrix$PC1, 
  "pc2" = PCs_pcair_matrix$PC2,
  "pc3" = PCs_pcair_matrix$PC3,
  "pc4" = PCs_pcair_matrix$PC4, 
  "pheno" = oncoth2_pheno$VTE
)
# Create ScanAnnotationDataFrame object
scanAnnot <- ScanAnnotationDataFrame(data = scan_data)


# ------------------------------------------------- # 
#           Estimate kinship coefficients
# ------------------------------------------------- # 

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
  covars = c("age", "sex", "pc1", "pc2", "pc3", "pc4"),
  cov.mat = grm,
  family = "binomial"
  )

# Run SNP-phenotype association tests
mixed_models <- assocTestSingle(
  gdsobj = geno_iterator, # establish how many SNPs are read at a time
  null.model = nullmodel, 
  BPPARAM = BiocParallel::SerialParam()
)

