#Extracting environmental variables from WorldClim database#
#J. Mattos, 2024

library(raster)
library(sp)

path <- "C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim"

# Load population coordinates
coords <- read.csv2('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/Coords.csv')  # Assuming columns 'latitude' and 'longitude'

longitude <- coords$longitude
latitude <- coords$latitude

#checking if data is numeric 

is.numeric(coords$longitude)
is.numeric(coords$latitude)

# Convert latitude and longitude columns to numeric
coords$latitude <- as.numeric(as.character(coords$latitude))
coords$longitude <- as.numeric(as.character(coords$longitude))

# Rounding coordinates
coords$longitude <- round(coords$longitude, digits = 5)
coords$latitude <- round(coords$latitude, digits = 5)


# Convert to SpatialPointsDataFrame
coordinates(coords) <- ~longitude+latitude
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")  # Define the CRS


# Load WorldClim rasters
#BIO1 annual mean temperature
BIO1 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_1.tif')

#BIO2 mean diurnal range
BIO2 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_2.tif')

#BIO3 isothermality
BIO3 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_3.tif')

#BIO4 temperature seasonality
BIO4 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_4.tif')

#BIO5 max temp of warmest month
BIO5 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_5.tif')

#BIO6
BIO6 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_6.tif')

#BIO7
BIO7 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_7.tif')

#BIO8
BIO8 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_8.tif')

#BIO9
BIO9 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_9.tif')

#BIO10
BIO10 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_10.tif')

#BIO11
BIO11 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_11.tif')

#BIO12 annual precipitation
BIO12 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_12.tif')

#BIO13 precipitation of wettest month 
BIO13 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_13.tif')

#BIO14
BIO14 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_14.tif')

#BIO15
BIO15 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_15.tif')

#BIO16
BIO16 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_16.tif')

#BIO17
BIO17 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_17.tif')

#BIO18 precipitation of warmest quarter
BIO18 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_18.tif')

#BIO19 
BIO19 <- raster('C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/WorldClim/wc2.1_10m_bio/wc2.1_10m_bio_19.tif')



#checking for NA in the rasters
summary(BIO1)
summary(BIO2)
summary(BIO3)
summary(BIO4)
summary(BIO5)
summary(BIO6)
summary(BIO7)
summary(BIO8)
summary(BIO9)
summary(BIO10)
summary(BIO11)
summary(BIO12)
summary(BIO13)
summary(BIO14)
summary(BIO15)
summary(BIO16)
summary(BIO17)
summary(BIO18)
summary(BIO19)


#checking CRS - coordinate systems
crs(coords)
crs(BIO1)
crs(BIO2)
crs(BIO3)
crs(BIO4)
crs(BIO5)

# Check the extent of the raster
raster_extent <- extent(BIO1)
print(raster_extent)

# Compare with the range of your coordinates
range(coords$longitude)
range(coords$latitude)

# Plotting the raster to verify the points and raster extent
plot(BIO1)
points(coords, col = "black")

# Buffering raster (if applicable)
buffered_raster_temp <- extend(annual_temp_raster, 10)  # Extends the raster by 10 units
buffered_raster_precip <- extend(annual_precipitation_raster, 10)


# Extracting the environmental predictors for each individual/coordinate that we have

coords$bio1 <- extract(BIO1, coords)
coords$bio2 <- extract(BIO2, coords)
coords$bio3 <- extract(BIO3, coords)
coords$bio4 <- extract(BIO4, coords)
coords$bio5 <- extract(BIO5, coords)
coords$bio6 <- extract(BIO6, coords)
coords$bio7 <- extract(BIO7, coords)
coords$bio8 <- extract(BIO8, coords)
coords$bio9 <- extract(BIO9, coords)
coords$bio10 <- extract(BIO10, coords)
coords$bio11 <- extract(BIO11, coords)
coords$bio12 <- extract(BIO12, coords)
coords$bio13 <- extract(BIO13, coords)
coords$bio14 <- extract(BIO14, coords)
coords$bio15 <- extract(BIO15, coords)
coords$bio16 <- extract(BIO16, coords)
coords$bio17 <- extract(BIO17, coords)
coords$bio18 <- extract(BIO18, coords)
coords$bio19 <- extract(BIO19, coords)


# Save as dataframe
coords_df <- as.data.frame(coords)

# Identifying the NA values
na_coords <- coords_df[is.na(coords_df$temp) | is.na(coords_df$precip), ]
print(na_coords)


# Save the results
write.csv(as.data.frame(coords_df), 'C:/Users/JacquelineMattos/Documents/Docs_Jac/Doutorado/Analyses/PopulationGenomics/Population_Genomics/Enviromental_Variables/extracted_ALL_BIOCLIM_10min.csv', row.names = FALSE)









