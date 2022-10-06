# load packages

#install.packages("raster")
#install.packages("dplyr")
#install.packages("caret")
library(raster)
library(dplyr)
library(caret)

setwd("C:/Users/u0142858/OneDrive - KU Leuven/KUL/PhD/My Project/WP2_Map_SB/Rsources/GitHub/SDMs test")
# set directory where the cleaned occurrence data are stored (adapt this)
occ_clean_dir <- "./ConservePlants/DataCleaning"

# set directory where environmental data are stored (adapt this)
env_dir <- "./ConservePlants/EnvironmentalData"

# set directory where to store the models
SDM_dir <- "./ConservePlants/SDM"

# set directory where to store the maps
output_dir <- "./ConservePlants/Output"

##########################################################
######             Fitting an SDM (GLM)              #####
##########################################################

# load the occurrence data files
species_names <- "Campanula latifolia" # choose your own species
occurrences <-  list.files(path = occ_clean_dir, pattern = "*.csv", full.names = TRUE)

# load predictor variables in a rasterstack
setwd(dir = env_dir)
covariates_current <- stack(list.files())

# delete variables with Spearman correlation > 0.7
covariates_current <- covariates_current[[c("BIO_02",
                                            "BIO_07",
                                            "BIO_08",
                                            "BIO_09",
                                            "BIO_13",
                                            "BIO_14",
                                            "BIO_18",
                                            "BIO_19",
                                            "cfvo_final",
                                            "clay_final",
                                            "silt_final",
                                            "soc_final")]]

# model fitting 
setwd("C:/Users/u0142858/OneDrive - KU Leuven/KUL/PhD/My Project/WP2_Map_SB/Rsources/GitHub/SDMs test")

for (q in 1:length(species_names)) {
  
  print(q)
  
  # prepare the occurrence points
  occurrence_points <- 
    read.csv(occurrences[q], header = TRUE, sep = ",") %>%
    SpatialPointsDataFrame(coords = .[,  c("decimalLongitude", "decimalLatitude")], 
                           data = ., proj4string = CRS("+init=epsg:3035")) %>% # convert to SpatialPointsDataFrame
    raster::extract(covariates_current, .) %>% # extract environmental variables at presence locations
    as.data.frame() %>% # convert to dataframe
    mutate(presence = 1) # add column indicating presence
  
  # prepare pseudo-absence points by randomly sampling
  pseudo_abscence_points <-
    sampleRandom(covariates_current, size = 10000) %>%
    as.data.frame() %>%
    mutate(presence = 0)
  
  # put the presence and pseudo-absence points together
  species_data <- rbind(occurrence_points, pseudo_abscence_points)
  
  # convert presence/absence column to factor
  species_data$presence <- as.factor(species_data$presence)
  
  # prepare cross-validation folds (10-fold cross validation)
  fitControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
  
  # 10-fold cross-validation of binomial GLM using the caret package
  # note that this not the ideal method as we're working with presence-only data
  glmFit <- train(presence ~ BIO_02 + BIO_07 + BIO_08 + soc_final,
                  data = species_data,
                  method = "glm", 
                  family = "binomial", #Why it is binominal here?
                  trControl = fitControl,
                  na.action = na.exclude)
  
  # get a summary of the model:
  summary(glmFit) 
  
  # check accuracy of the model
  glmFit
  
  # save your model results
  save(glmFit, file = paste0(SDM_dir, "/", species_names[q], ".RData"))
  
  # make predictions
  species_prediction <- predict(covariates_current, glmFit)
  
  # plot
  plot(species_prediction)
  
  # thresholding to produce a binary map
  species_prediction[species_prediction >= 0.8] <- 1
  species_prediction[species_prediction < 0.8] <- 0
  plot(species_prediction)
  
  # save the map
  writeRaster(species_prediction, paste0(output_dir, "/", species_names[q], ".tif"), overwrite = TRUE)
  
}
