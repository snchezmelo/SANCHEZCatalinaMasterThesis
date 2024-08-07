
##This script allows to calculate the beta diversity of raster files.
#Original by SANCHEZ Catalina.

library(BAT)
library(raster)
library(rgdal)

#Path to the models
raster_dir <- "./Resampled_F100_ReclassifiedPass"
#Pattern to recognize the models
raster_pattern <- "*\\.tif$"
#Read them as a list
raster_files <- list.files(path = raster_dir, pattern = raster_pattern, full.names = TRUE)
#Stack them
rasters <- stack(raster_files)
#Calculate beta diversity
beta_diversity <- raster.beta(rasters)
#Path to the output
output_file <- "./beta_diversity_map.tif"
#Save the beta diversity map to a file
writeRaster(beta_diversity, filename = output_file, format = "GTiff", overwrite = TRUE)