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

# Load PLINK files
plink_files_dir_path <- "/mnt/ir-bioinf02/home/blobato/oncothromb02/GWAS/snpQT_results/results/post_imputation/bfiles/"
# Using all SNPs
# plink_files_name <- "H7" 
# Use SNPs resulting from pruning
plink_files_name <- "pruned_output"

# Load principal components from PCA done with PLINK1.9
pca_dir_path <- "/mnt/ir-bioinf02/home/blobato/oncothromb02/GWAS/snpQT_results/results/post_imputation/bfiles/"
pca_file_path <- "pca_results.eigenvec"
pca_eigenvec_path <- paste0(pca_dir_path, pca_file_path)

pc <- read.table(
  pca_eigenvec_path, 
  header = FALSE, 
  col.names = c("FID", "IID", "PC1", "PC2", "PC3", 
                "PC4", "PC5", "PC6", "PC7", "PC8",
                "PC9", "PC10", "PC11", "PC12", "PC113",
                "PC14", "PC15", "PC16", "PC17", "PC18",
                "PC19", "PC20"), 
  colClasses = c("character", "character", rep("numeric", times = 20))
)


# # Add patient_code to table
# pc %<>%
#   mutate(IID = paste(FID, IID, sep = "_")) %>%
#   mutate(FID = 0)


# Use ID from an9elproject object
oncoth2 <- get_project("oncoth2") # , version = "0.0.8008"
oncoth2_ids <- oncoth2$data %>%
  select(id, patient_code)





# We think reasonable to use the first 4 PCs as covariates to correct for
# unknown effects caused by genetic variability, as well as ancestry (some individuals
# don't have European ancestry but we are going to try to keep them in the analysis) 


# Create ScanAnnotationDataFrame 
# Contains outcome and covariate data, as well as unique identifier (scanID)
# It can be paired with genotype

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
  # # Make sure to leave out those with
  # filter(informed_consent == "Yes") %>%
  # filter(eligible == "Yes") %>%
  # filter(patient_code %in% pc$IID) %>%
  select(id, patient_code, VTE, age_cancer_dx, sex)

oncoth2_pheno <- right_join(oncoth2_pheno, pc, by = c("patient_code" = "IID"))

# --- Remove when patients are excluded beforehand --- #

# Proportion of VTE
# print(table(oncoth2_pheno$VTE))

print(paste("There are", dim(oncoth2_pheno)[1], "patients in ONCOTHROMB2 with clinical information that can be analysed."))
print(paste("There are only", dim(pc)[1], "individuals genotyped in ONCOTHROMB2."))

# # Get clinical info for genotyped patients
# oncoth2_pheno_genotyped <- oncoth2_pheno %>%
#   filter(patient_code %in% pc$IID)

# # Remove patients that are not eligible or haven't signed the IC
pc %<>%
  # Careful! Since IID and patient code are strings, to properly compare them
  # they have to be ordered
  filter(order(IID) %in% order(oncoth2_pheno$patient_code))

# Transform to matrix
pc_matrix <- pc
rownames(pc_matrix) <- pc_matrix$IID
pc_matrix %<>% select(-c(FID, IID))
pc_matrix <- as.matrix(pc_matrix)


# Create dataframe for ScanAnnotationDataFrame
scan_data <- data.frame(
  "scanID" = oncoth2_pheno$patient_code, 
  "age" = oncoth2_pheno$age_cancer_dx, 
  "sex" = ifelse(oncoth2_pheno$sex == 0, "M", "F"), # valid format for ScanAnnotationDataFrame
  "pc1" = pc$PC1, 
  "pc2" = pc$PC2,
  "pc3" = pc$PC3,
  "pc4" = pc$PC4, 
  "pheno" = oncoth2_pheno$VTE
)
# Create ScanAnnotationDataFrame object
scanAnnot = ScanAnnotationDataFrame(data = scan_data)



# Check if previous GDS file exists
gds_file <- paste0(plink_files_dir_path, "genotype.gds")

# if (!file.exists(gds_file)) {
#   # Convert PLINK files to GDS to handle genotype data
#   snpgdsBED2GDS(
#     bed.fn = paste0(plink_files_dir_path, plink_files_name, ".bed"), # "./GWAS/imputation_results/oncoth2.all_chr.TOPMED.hg38.bed",
#     fam.fn = paste0(plink_files_dir_path, plink_files_name, ".fam"),
#     bim.fn = paste0(plink_files_dir_path, plink_files_name, ".bim"),
#     out.gdsfn = gds_file
#   )
# } else {
#   print("GDS file already exists.")
# }

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


# Include Genetic Relationship Matrix (GRM) to account for genetic
# similarity among sample individuals

# Create GenotypeBlockIterator class to iterate over blocks of 10,000 SNPs
geno_iterator <- GenotypeBlockIterator(genotype_data)

# Do PCA on genome-wide SNP data
# The principal components account for population structure in a sample
pca <- pcair(
  gdsobj = geno_iterator, 
  num.cores = num_cores
)

# Get principal components
PCs_pcair_matrix <- pca$vectors

# Estimate kinship coefficients
# Model-free estimation of recent genetic relatedness
pcrel_result <- pcrelate(
  gdsobj = geno_iterator, 
  pcs = PCs_pcair_matrix
)

# Genetic relationship matrix
grm <- pcrelateToMatrix(pcrel_result)
# grm[1:5, 1:5]




# Fit null model
# Mixed model that contains all covariates but not genotype data
# Null hypothesis is that all SNPs have no effect on outcome
nullmodel = fitNullModel(
  scanAnnot, family = "binomial")



