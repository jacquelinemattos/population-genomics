##Script to plot environmental variables along the species/populations gradient ##

# Install necessary packages if not already installed
install.packages(c("raster", "rgdal", "ggplot2", "dplyr"))

# Load the packages
library(raster)
library(rgdal)
library(ggplot2)
library(dplyr)

sample_file <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/pop_ind_lat_long_env_data_ALL_BIOCLIM.csv")

#Converting sample coordinates into SpatialPoints object

coords <- sample_file[, c("longitude", "latitude")]
coords$latitude <- as.numeric(coords$latitude)
coords$longitude <- as.numeric(coords$longitude)
coordinates(coords) <- ~longitude + latitude

# Extract the bioclimatic variables for each sample location
#bioclim_values <- extract(sample_file[,5:23], coords)

#Na verdade eu ja tinha extraido os valores para cada localidade. 
#Os valores da minha tabela sample_file ja sao extraidos para minhas coordenadas. 

bioclim <- sample_file[,5:23]
summary(sample_file$bio1)

# Plot for Bio1 (Annual Mean Temperature) along the latitude gradient
ggplot(sample_file, aes(x = latitude, y = bio1, color = population)) +
  geom_point() +
  geom_smooth(method = 'loess') +
  labs(title = "Annual Mean Temperature (Bio1) vs. Latitude",
       x = "Latitude", 
       y = "Annual Mean Temperature (Bio1)") +
  theme_minimal()


# Plot for Bio12 (Annual Precipitation) along the latitude gradient
ggplot(sample_file, aes(x = latitude, y = bio12, color = population)) +
  geom_point() +
  geom_smooth(method = 'loess') +
  labs(title = "Annual Precipitation (Bio12) vs. Latitude",
       x = "Latitude", 
       y = "Annual Precipitation (Bio12)") +
  theme_minimal()


# Plot for Bio4 along the latitude gradient
ggplot(sample_file, aes(x = latitude, y = bio4, color = population)) +
  geom_point() +
  geom_smooth(method = 'loess') +
  labs(title = "Bio4 vs. Latitude",
       x = "Latitude", 
       y = "Bio4") +
  theme_minimal()


# Plot for Bio13 along the latitude gradient
ggplot(sample_file, aes(x = latitude, y = bio13, color = population)) +
  geom_point() +
  geom_smooth(method = 'loess') +
  labs(title = "Bio13 vs. Latitude",
       x = "Latitude", 
       y = "Bio13") +
  theme_minimal()


# Plot for Bio18 along the latitude gradient
ggplot(sample_file, aes(x = latitude, y = bio18, color = population)) +
  geom_point() +
  geom_smooth(method = 'loess') +
  labs(title = "Bio18 vs. Latitude",
       x = "Latitude", 
       y = "Bio18") +
  theme_minimal()


# Plot for Bio2 along the latitude gradient
ggplot(sample_file, aes(x = latitude, y = bio2, color = population)) +
  geom_point() +
  geom_smooth(method = 'loess') +
  labs(title = "Bio2 vs. Latitude",
       x = "Latitude", 
       y = "Bio2") +
  theme_minimal()


# Plot for Bio8 along the latitude gradient
ggplot(sample_file, aes(x = latitude, y = bio8, color = population)) +
  geom_point() +
  geom_smooth(method = 'loess') +
  labs(title = "Bio8 vs. Latitude",
       x = "Latitude", 
       y = "Bio8") +
  theme_minimal()


# Plot for Bio14 along the latitude gradient
ggplot(sample_file, aes(x = latitude, y = bio14, color = population)) +
  geom_point() +
  geom_smooth(method = 'loess') +
  labs(title = "Bio14 vs. Latitude",
       x = "Latitude", 
       y = "Bio14") +
  theme_minimal()








