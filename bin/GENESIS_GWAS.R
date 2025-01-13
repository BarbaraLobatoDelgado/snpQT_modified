# 2024-02-08
# GWAS in ONCOTHROMB2 cohort
# Trait CAT vs no CAT

# Sources
# GWAS_Dosage.R (author: Gerard Temprano-Sagrera)
# https://www.bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/assoc_test.html


# Load libraries
pacman::p_load(
  an9elproject,
  data.table,
  tidyverse, 
  magrittr,
  GWASTools,
  SNPRelate,
  GENESIS
  )

# Detect cores
# parallel::detectCores()



# Load principal components from PCA done with PLINK1.9
pca_dir_path = "/mnt/ir-bioinf02/home/blobato/oncothromb02/GWAS/imputation_results/"
pca_file_path = "oncoth2.all_chr.TOPMED.hg38.pruned.pca.eigenvec"
pca_eigenvec_path = paste0(pca_dir_path, pca_file_path)

pc = read.table(
  pca_eigenvec_path, 
  header = FALSE, 
  col.names = c("FID", "IID", "PC1", "PC2", "PC3", 
                "PC4", "PC5", "PC6", "PC7", "PC8",
                "PC9", "PC10", "PC11", "PC12", "PC113",
                "PC14", "PC15", "PC16", "PC17", "PC18",
                "PC19", "PC20"), 
  colClasses = c("character", "character", rep("numeric", times = 20))
  )


# Add patient_code to table
pc %<>%
  mutate(IID = paste(FID, IID, sep = "_")) %>%
  mutate(FID = 0)


# Use ID from an9elproject object
oncoth2 = get_project("oncoth2", version = "0.0.8008")
oncoth2_ids = oncoth2$data %>%
  select(id, patient_code)


test = right_join(oncoth2_ids, pc, by = c("patient_code" = "IID"))


# We think reasonable to use the first 4 PCs as covariates to correct for
# unknown effects caused by genetic variability, as well as ancestry (some individuals
# don't have European ancestry but we are going to try to keep them in the analysis) 


# Create ScanAnnotationDataFrame 
# Contains outcome and covariate data, as well as unique identifier (scanID)
# It can be paired with genotype


# Load ONCOTHROMB2 DB
oncoth2 = get_project("oncoth2", version = "0.0.8007")

# Select dependent variable and covariates
oncoth2_dat = oncoth2$data %>%
  # Make sure to leave out those with 
  filter(informed_consent == 1) %>%
  filter(eligible == 1) %>%
  # filter(patient_code %in% pc$IID) %>%
  select(id, patient_code, VTE, age_cancer_dx, sex) 

# Proportion of VTE
# print(table(oncoth2_dat$VTE))

print(paste("There are", dim(oncoth2_dat)[1], "patients in ONCOTHROMB2 with clinical information that can be analysed."))
print(paste("There are only", dim(pc)[1], "individuals genotyped in ONCOTHROMB2."))


# Remove patients that are not eligible or haven't signed the IC
pc %<>%
  # Careful! Since IID and patient code are strings, to properly compare them 
  # they have to be ordered
  filter(order(IID) %in% order(oncoth2_dat$patient_code))

# Transform to matrix
pc_matrix = pc
rownames(pc_matrix) = pc_matrix$IID
pc_matrix %<>% select(-c(FID, IID))
pc_matrix = as.matrix(pc_matrix)


# Create dataframe for ScanAnnotationDataFrame
scan_data = data.frame(
  "scanID" = oncoth2_dat$patient_code, 
  "age" = oncoth2_dat$age_cancer_dx, 
  "sex" = ifelse(oncoth2_dat$sex == 0, "M", "F"),
  "pc1" = pc$PC1, 
  "pc2" = pc$PC2,
  "pc3" = pc$PC3,
  "pc4" = pc$PC4, 
  "pheno" = oncoth2_dat$VTE
)

scanAnnot = ScanAnnotationDataFrame(data = scan_data)


# Read genotype data
# Convert PLINK files to GDS
snpgdsBED2GDS(
  bed.fn = "./GWAS/imputation_results/oncoth2.all_chr.TOPMED.hg38.bed",
  fam.fn = "./GWAS/imputation_results/oncoth2.all_chr.TOPMED.hg38.fam",
  bim.fn = "./GWAS/imputation_results/oncoth2.all_chr.TOPMED.hg38.bim",
  out.gdsfn = "./GWAS/association_studies/oncoth2_genotype.gds"
  )

# Load GDS file
genotype_file = GdsGenotypeReader("./GWAS/association_studies/oncoth2_genotype.gds")
genotype_data = GenotypeData(genotype_file)

print(getScanID(genotype_data))
print(nsnp(genotype_data))
print(getSnpID(genotype_data))


# Include Genetic Relationship Matrix (GRM) to account for genetic
# similarity among sample individuals

# Create GenotypeBlockIterator
geno_iterator = GenotypeBlockIterator(genotype_data)

# Do pcrelate()
pcrel_result = pcrelate(geno_iterator, pcs = pc_matrix)



# Fit null model
nullmodel = fitNullModel(, family = "binomial")





