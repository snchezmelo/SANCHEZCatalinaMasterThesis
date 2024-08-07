#This code allows to model the potential distribution of species based on occurrence data.
#Adapted from a *Wallace* v2.1.1 session by SANCHEZ Catalina.
### Package installation

library(spocc)
library(spThin)
library(dismo)
library(sf)
library(ENMeval)
library(wallace)
library(dplyr)
library(tibble)


#Before we start, we need to create a vector with the names of the layers
layers <- c('bio1.tif','bio2.tif', 'bio3.tif', 'bio4.tif', 'bio5.tif','bio6.tif', 'bio7.tif', 'bio8.tif', 'bio9.tif', 'bio10.tif', 'bio11.tif', 'bio12.tif', 'bio13.tif', 'bio14.tif','bio15.tif', 'bio16.tif', 'bio17.tif','bio18.tif','bio19.tif')

#Then we also need to charge an R object containing a list of the spp abbreviation and complete names
PassifloraNames <- read.csv("PassNames.csv", sep = ";")

# NOTE: provide the folder path of the .csv file
occs_path <- "/Users/snchezmelo/Desktop/GEUR"
occs_path <- file.path(occs_path, "PassOccs.csv")
# get a list of species occurrence data
userOccs <- occs_userOccs(
  txtPath = occs_path,
  txtName = "PassOccs.csv",
  txtSep = ";",
  txtDec = ".")

### Obtain environmental data

## Specify the directory with the environmental variables
dir_envs <- "/Users/snchezmelo/Desktop/GEUR/AllBioClim"

envs_path <- file.path(dir_envs, layers)

envs <- envs_userEnvs(
  rasPath = envs_path,
  rasName = layers,
  doBrick = FALSE)

for (species in PassifloraNames$Species){
  tryCatch({
    # Construct the variable names dynamically
    # occs_envs_var <- paste0("envs_", species_abbr)
    #occs_cleaned_var <- paste0("occs_", species_abbr)
    
    # Extract xy data
    occs_xy<- userOccs[[species]]$cleaned[c('longitude', 'latitude')]
    #occs_xy <- get(paste0("occs_", species_abbr))[c('longitude', 'latitude')]
    
    # Extract values
    occs_vals <- as.data.frame(raster::extract(envs, occs_xy, cellnumbers = TRUE))
    
    # Remove duplicated cell values
    occs <- userOccs[[species]]$cleaned
    occs <- occs[!duplicated(occs_vals[, 1]), ]
    occs_vals <- occs_vals[!duplicated(occs_vals[, 1]), -1]
    
    # Remove occurrence records with NA environmental values
    occs <- occs[!(rowSums(is.na(occs_vals)) >= 1), ]
    
    # Also remove variable value rows with NA environmental values
    occs_vals <- na.omit(occs_vals)
    
    # Add columns for env variable values for each occurrence record
    occs <- cbind(occs, occs_vals)
    
    # Update the variables in the environment
    #assign(occs_cleaned_var, occs)
    
    # Thin occurrences 
    thinDist <- 1
    occs <- poccs_thinOccs(occs = occs, thinDist = thinDist)

    # Generate background extent
    bgExt <- penvs_bgExtent(occs = occs, bgSel = "point buffers", bgBuf = 1)

    # Mask environmental data to provided extent
    bgMask <- penvs_bgMask(occs = occs, envs = envs, bgExt = bgExt)

    # Sample background points from the provided area
    bgSample <- penvs_bgSample(occs = occs, bgMask = bgMask, bgPtsNum = 10000)

    # Extract values of environmental layers for each background point
    bgEnvsVals <- as.data.frame(raster::extract(bgMask, bgSample))
    
    # Add extracted values to background points table
    bgEnvsVals <- cbind(scientific_name = species, bgSample, occID = NA, year = NA,
                        institution_code = NA, country = NA, state_province = NA, locality = NA,
                        elevation = NA, record_type = NA, bgEnvsVals)

    # Check the count of occurrences
    if (nrow(occs) < 30) {
      method <- "jack"
    } else {
      method <- "block"
    }
    
    # R code to get partitioned data with the selected method
    groups <- part_partitionOccs(occs = occs, bg = bgSample, method = method, bgMask = bgMask, aggFact = 2)

    # Run maxent model for the selected species
    model <- model_maxent(occs = occs, bg = bgEnvsVals, user.grp = groups, bgMsk = bgMask, rms = c(1, 4),
                          rmsStep = 0.5, fcs = c('L', 'LQ', 'LQH'), clampSel = TRUE, algMaxent = "maxnet", parallel = TRUE,
                          numCores = 7)
    
    results_ordered <- model@results[order(model@results$or.mtp.avg, -model@results$auc.val.avg),]
    write.csv(results_ordered, file = paste0("/Users/snchezmelo/Desktop/GEUR/WallaceResults/results_", species, ".csv"), quote = FALSE)
    #Saving the model itself
    best_modelOR <- as.integer(rownames(results_ordered[1,]))
    m <- model@models[[best_modelOR]]
    #saveRDS(m, file = paste0("/Users/snchezmelo/Desktop/GEUR/WallaceResults/m_", species, ".rds"))
    #Saving the raster
    predSel <- predictMaxnet(m, envs,
                             type = "cloglog",
                             clamp = TRUE)
    writeRaster(predSel, file = paste0("/Users/snchezmelo/Desktop/GEUR/WallaceResults/predSel_", species, ".tif"),overwrite=TRUE)
    #Saving the used variables
    all_used_variables <- data.frame()
    ##From Andre's code
    all_used_variables <- t(as.data.frame(m$betas))
    write.csv(all_used_variables, file = paste0("/Users/snchezmelo/Desktop/GEUR/WallaceResults/usedVar_", species, ".csv"), quote = FALSE)
    
    
    
  }, error = function(e) {
    cat("Error processing species:", species, "\n")
    cat("Error message:", e$message, "\n")
  })
  
}
