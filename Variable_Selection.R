##Script to perform varible selection/reduction and remove collinearity ##

#### Multicolinearity analyses ####

library(raster)
library(usdm)

env_allBIOCLIM <- read.csv2("~/Documents/Jac/Population Genomics/Tabelas/pop_ind_lat_long_env_data_ALL_BIOCLIM.csv", header = TRUE)
bioclim <- env_allBIOCLIM[, 5:23]
bioclim <- as.matrix(bioclim)

corr <- vifcor(bioclim, th=0.7)
#kept bio1, bio2, bio8 and bio14

vif <- vifstep(bioclim, th = 7, keep = NULL, method = "pearson")
