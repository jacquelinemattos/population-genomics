#PCA of environmental variables#
#J. Mattos, 2024

install.packages("raster")   # For handling BIOCLIM raster data
install.packages("sp")       # Spatial data handling
install.packages("stats")    # For PCA (this is base R)
install.packages("ggplot2")  # For plotting PCA results
install.packages("factoextra") # Optional for enhanced PCA visualization

library(raster)
library(sp)
library(ggplot2)
library(factoextra)

#This part here is just if you didn't have your environmental variables already extracted for each coordinate/individual
#Go to next section if variables are already extracted

#### Load the 19 BIOCLIM variables (Global dataset) from .tif or .asc files saved on a specific folder ####

# Set working directory
setwd("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/")

# List all BIOCLIM raster files in the folder, assuming they have the `.tif` or `.asc` extension
raster_files <- list.files(pattern = "\\.tif$|\\.asc$")
raster_files

#Now, use stack() or lapply() to load and combine all the rasters into a single RasterStack object.
bioclim_stack <- stack(raster_files)

# View basic information about the stack
print(bioclim_stack)

# Check the names of the layers (each corresponds to a BIOCLIM variable)
names(bioclim_stack)

# Creating two stacks with only temperature and precipitation variables in each 

temp_stack <- stack(bioclim_stack[[c(1,12,13,14,15,16,17,18,19,2,3)]])
names(temp_stack)
plot(temp_stack)

precip_stack <- stack(bioclim_stack[[c(4,5,6,7,8,9,10,11)]])
names(precip_stack)
plot(precip_stack)


#### Prepare the Data for PCA ####

#Before applying PCA, ensure that the data is in a format suitable for analysis (i.e., matrix or data frame). Convert the raster stack to a matrix, remove missing values, and scale the data.

# Convert the raster stack into a data frame (may be large depending on resolution)
temp_values <- values(temp_stack)
precip_values <- values(precip_stack)

# Remove rows with missing data
temp_values <- na.omit(temp_values)
precip_values <- na.omit(precip_values)

# Optionally, scale the data (PCA works best with scaled data)
temp_values_scaled <- scale(temp_values)
precip_values_scaled <- scale(precip_values)


#### Importing extracted environmental variables ####

BIOCLIM <- read.csv("C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/extracted_ALL_BIOCLIM_10min.csv", header = TRUE)
head(BIOCLIM)

# Removing the row that contains the coordinate for Morro Santana (not actually using it in the other analyses!)
BIOCLIM <- BIOCLIM[-6,]
BIOCLIM

bioclim_temperature <- BIOCLIM[,1:11]
bioclim_precipitation <- BIOCLIM[,12:19]

# Perform PCA with prcomp
pca_temperature <- prcomp(bioclim_temperature, center = TRUE, scale. = TRUE)
pca_precipitation <- prcomp(bioclim_precipitation, center = TRUE, scale. = TRUE)


#Inspect PCA Results 
#After running PCA, it's useful to examine the results, such as the proportion of variance explained by each principal component and the loadings.

# Summary of PCA
summary(pca_temperature)
summary(pca_precipitation)

# Scree plot to visualize variance explained by each component
plot(pca_temperature, type = "l")
plot(pca_temperature)

plot(pca_precipitation, type = "l")
plot(pca_precipitation)


# Biplot of the first two principal components
biplot(pca_temperature, scale = 0)
biplot(pca_precipitation, scale = 0)


# Plotting PCA with ggplot2
# Extract the PCA scores for the individuals (rows)
pca_scores_temp <- as.data.frame(pca_temperature$x)
pca_scores_precip <- as.data.frame(pca_precipitation$x)


# Plot the first two principal components (PC1 and PC2)
ggplot(pca_scores_precip, aes(x = PC1, y = PC2)) +
  geom_point() +
  xlab(paste0("PC1 (", round(summary(pca_scores_precip)$importance[2, 1] * 100, 2), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca_scores_precip)$importance[2, 2] * 100, 2), "%)")) +
  ggtitle("PCA: First Two Principal Components") +
  theme_minimal()

# Separating PCS 1 and 2 to use in other analyses as predictor variables 

PC1_temp <- pca_scores_temp$PC1
PC2_temp <- pca_scores_temp$PC2

PC1_precip <- pca_scores_precip$PC1
PC2_precip <- pca_scores_precip$PC2

# Save the PCs to use in other analyses 

PCs <- cbind(PC1_temp, PC2_temp, PC1_precip, PC2_precip, BIOCLIM$longitude, BIOCLIM$latitude)
PCs_df <- as.data.frame(PCs)

PCs_comma[] <- lapply(PCs_df, function(x) gsub("\\.", ",", x))
PCs_comma <- as.data.frame(PCs_comma)

write.csv(PCs_comma, "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Gene-Environment-Association/RDA/PCAs/environmental_PCs_withcoords.csv", row.names = FALSE)

# Saving as Rdata file too 

save(PC1_temp, PC2_temp, PC1_precip, PC1_precip, file = "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Gene-Environment-Association/RDA/PCAs/Environmental_PCs.RData")

#PC1 and PC2 are now separate vectors that you can use in regression models, clustering, or any other analysis.
#The PCA scores can also be visualized in scatterplots or used as features in machine learning algorithms.









