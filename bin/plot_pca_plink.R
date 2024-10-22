#!/usr/bin/env Rscript

library("tidyverse")
library("plotly")
library("htmlwidgets")

# Args
# 1: eigenvec
# 2: popfile 
# 3: before or after removal

args <- commandArgs(trailingOnly = TRUE)

eigenvec <- read_delim(args[[1]], delim = " ")
popfile <- read_delim(args[[2]], delim = " ")
# important to override colnames to make super / sub pop file input consistent
colnames(popfile) <- c("IID", "IID2", "pop")

eigenvec %>%
    left_join(popfile, by = "IID") %>%
    # missing pop means it's user data
    replace_na(list(pop = "OWN")) -> datafile


# Set colors to each population ancestry
pop_colors <- c("AFR" = "#F8766D", 
                "AMR" = "#EDAE49", 
                "EAS" = "#00BE67",       
                "EUR" = "lightblue",
                "OWN" = "steelblue",
                "SAS" = "#A58AFF")

# Set alpha values for each population ancestry
pop_alpha <- c("AFR" = 0.4, 
               "AMR" = 0.4, 
               "EAS" = 0.4,       
               "EUR" = 0.4,
               "OWN" = 0.8,
               "SAS" = 0.4)

# Set point shape for each population ancestry
pop_shape <- c("AFR" = 16, 
               "AMR" = 16, 
               "EAS" = 16,       
               "EUR" = 16,
               "OWN" = 18,
               "SAS" = 16)


datafile %>%
  ggplot() +
  # First layer: plot all points except for 'OWN'
  geom_point(data = subset(datafile, pop != "OWN"), 
             aes(x = PC1, y = PC2, color = pop, alpha = pop), 
             shape = 16, size = 2) +  # Use the same shape and size for all non-OWN populations
  # Second layer: plot only 'OWN' points with a different shape and size
  geom_point(data = subset(datafile, pop == "OWN"), 
             aes(x = PC1, y = PC2, color = pop, alpha = pop), 
             shape = 18, size = 3) +  # Different shape and size for 'OWN'
  # Manual scales for color and alpha
  scale_color_manual(values = pop_colors) +  
  scale_alpha_manual(values = pop_alpha) +   
  # Custom theme and labels
  theme_linedraw()
ggsave(paste0("PC1vsPC2_",args[[3]],".png"))


datafile %>%
  ggplot() +
  # First layer: plot all points except for 'OWN'
  geom_point(data = subset(datafile, pop != "OWN"), 
             aes(x = PC1, y = PC3, color = pop, alpha = pop), 
             shape = 16, size = 2) +  # Use the same shape and size for all non-OWN populations
  # Second layer: plot only 'OWN' points with a different shape and size
  geom_point(data = subset(datafile, pop == "OWN"), 
             aes(x = PC1, y = PC3, color = pop, alpha = pop), 
             shape = 18, size = 3) +  # Different shape and size for 'OWN'
  # Manual scales for color and alpha
  scale_color_manual(values = pop_colors) +  
  scale_alpha_manual(values = pop_alpha) +   
  # Custom theme and labels
  theme_linedraw()
ggsave(paste0("PC1vsPC3_", args[[3]], ".png"))


datafile %>%
  ggplot() +
  # First layer: plot all points except for 'OWN'
  geom_point(data = subset(datafile, pop != "OWN"), 
             aes(x = PC2, y = PC3, color = pop, alpha = pop), 
             shape = 16, size = 2) +  # Use the same shape and size for all non-OWN populations
  # Second layer: plot only 'OWN' points with a different shape and size
  geom_point(data = subset(datafile, pop == "OWN"), 
             aes(x = PC2, y = PC3, color = pop, alpha = pop), 
             shape = 18, size = 3) +  # Different shape and size for 'OWN'
  # Manual scales for color and alpha
  scale_color_manual(values = pop_colors) +  
  scale_alpha_manual(values = pop_alpha) +   
  # Custom theme and labels
  theme_linedraw()
ggsave(paste0("PC2vsPC3_", args[[3]], ".png"))


fancy_plot <- plot_ly(
  datafile,
  x =  ~ PC1,
  y =  ~ PC2,
  z =  ~ PC3,
  size = 3,
  type = "scatter3d",
  mode = "markers",
  color =  ~ pop
)
saveRDS(fancy_plot, paste0("plink_3D_pca", args[[3]], ".rds"))
