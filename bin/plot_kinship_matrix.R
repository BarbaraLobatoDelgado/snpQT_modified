#!/usr/bin/env Rscript

library("data.table")
library("tidyverse")

# Args
# 1: kingship matrix file path

args <- commandArgs(trailingOnly = TRUE)
kinship_file <- args[[1]]

# Read the kinship data from the .kin0 file (ensure proper column headers)
kinship_data <- fread(kinship_file, header = TRUE)

# Debug: Print the first few rows to make sure the data is read correctly
print("Showing head of kinship data...")
print(head(kinship_data))

# Extract unique IDs from ID1 and ID2
all_ids <- unique(c(kinship_data$ID1, kinship_data$ID2))

# Create an empty kinship matrix with IDs as row and column names
kinship_matrix <- matrix(NA, nrow = length(all_ids), ncol = length(all_ids),
                         dimnames = list(all_ids, all_ids))

# Check if the matrix is initialized correctly
if (nrow(kinship_matrix) == 0 || ncol(kinship_matrix) == 0) {
  stop("Kinship matrix could not be created: matrix dimensions are zero.")
}

# Fill the matrix with kinship values
for (i in 1:nrow(kinship_data)) {
  id1 <- kinship_data$ID1[i]
  id2 <- kinship_data$ID2[i]
  kinship_matrix[id1, id2] <- kinship_data$KINSHIP[i]
  kinship_matrix[id2, id1] <- kinship_data$KINSHIP[i] # Fill symmetrically
}

# Debug: Ensure the matrix has the expected dimensions
print(paste("Kinship matrix dimensions:", dim(kinship_matrix)))

# Safely check a subset of the matrix
if (nrow(kinship_matrix) >= 5 && ncol(kinship_matrix) >= 5) {
  print(kinship_matrix[1:5, 1:5])  # Print the first 5 rows and columns
} else {
  print("Kinship matrix does not have enough rows/columns for subsetting")
}

# Replace NA values with 0 (assuming no kinship means 0)
kinship_matrix[is.na(kinship_matrix)] <- 0

# Convert the kinship matrix to a tibble and reshape to long format using tidyverse
kinship_long <- as.data.frame(kinship_matrix) %>%
  rownames_to_column(var = "ID1") %>%
  pivot_longer(cols = -ID1, names_to = "ID2", values_to = "Kinship")

# Debug: Print the first few rows of the long-format data to verify the transformation
print(head(kinship_long))

# Plot the heatmap using ggplot2
heatmap_plot <- ggplot(kinship_long, aes(x = ID1, y = ID2, fill = Kinship)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       name = "Kinship Coefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Kinship Matrix Heatmap", x = "Samples", y = "Samples")

# Save the heatmap to a file
ggsave("kinship_matrix_heatmap.png", plot = heatmap_plot, width = 10, height = 8)
