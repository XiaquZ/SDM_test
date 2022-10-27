library(dplyr) 
library(purrr) 
library(magrittr) 
library(rgbif) 
library(taxize) 
library(raster)
library(data.table)
library(CoordinateCleaner)
library(dismo)
library(mapr)
library(rgdal)
library(corrplot)

setwd("C:/Users/u0142858/OneDrive - KU Leuven/KUL/PhD/My Project/WP2_Map_SB/Rsources/GitHub/SDMs test")
# set directory where the downloaded occurrence data will be stored (adapt this)
occ_dir <- "./ConservePlants/GBIF/"

# set directory where the cleaned occurrence data will be stored (adapt this)
occ_clean_dir <- "./ConservePlants/DataCleaning/"

# set directory where environmental data are stored (adapt this)
env_dir <- "./ConservePlants/EnvironmentalData/"

##########################################################
######            Download data from GBIF            #####
##########################################################

# fill in your gbif.org credentials 
user <- "zhou068" # your gbif.org username 
pwd <- "Xiaqugbif123" # your gbif.org password
email <- "xiaqu.zhou@gmail.com" # your email 

# species for which you want to download
species_names <- "Campanula latifolia" # choose your own species here (alphabetical order)

# match species names with species names in GBIF database
gbif_taxon_keys <- 
  species_names %>% # use fewer names if you want to just test 
  taxize::get_gbifid_(method = "backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %T>% # combine all data.frames into one
  readr::write_csv(file = paste0(occ_dir, "/keys.csv")) %>% # save to inspect if you want
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% # get only accepted and matched names
  filter(kingdom == "Plantae") %>% # remove anything that might have matched to a non-plant
  pull(usagekey) # get the gbif taxonkeys
# check
gbif_taxon_keys 

# download occurrence data
occ_download(
  pred_in("taxonKey", gbif_taxon_keys), # important to use pred_in
  format = "SIMPLE_CSV",
  user = user,
  pwd = pwd,
  email = email
)

# download the records from the GBIF website and save them in the 01_GBIF directory
# under the name 'gbif.csv'
write.csv(d,"./ConservePlants/GBIF/gbif.csv", row.names = FALSE)

##########################################################
######          Cleaning occurrence records          #####
##########################################################

# load one of the environmental data rasters and reproject to WGS84
setwd(dir = env_dir)
BIO1 <- raster("BIO_01.tif") %>%
  projectRaster(crs = CRS("+init=epsg:4326"))

# plot to check
plot(BIO1)

# load the occurrence data
setwd(dir = occ_dir)
gbif_data <- read.csv("C:/Users/u0142858/OneDrive - KU Leuven/KUL/PhD/My Project/WP2_Map_SB/Rsources/GitHub/SDMs test/ConservePlants/GBIF/gbif.csv", header = TRUE, sep = ",")
# how many rows?
nrow(gbif_data)

# cleaning the occurrence data

# first remove impossible,missing and incomplete coordinates or coordinates outside Europe
# by extracting the BIO1 values at the coordinates and deleting the NA values
extr <- raster::extract(BIO1, cbind(gbif_data$decimalLongitude, gbif_data$decimalLatitude))
i <- which(is.na(extr))
gbif_data <- gbif_data[-i,]
# how many rows left?
nrow(gbif_data)

# further cleaning species by species (for-loop)
for (q in 1:length(species_names)) { 
  
  print(q)
  
  # several data cleaning steps
  species_of_interest <- 
    subset(gbif_data, species == species_names[q]) %>%
    dplyr::select("species", "decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters", "year", "basisOfRecord") %>% # select only relevant columns
    
    # I would here have coordinateUncertaintyInMeters so low that it ensures that the coordinate is likely within the raster cell it's supposed to represent
    # There are likely some good papers that you can use as references for this threshold.
    subset(coordinateUncertaintyInMeters <= 5000 | is.na(coordinateUncertaintyInMeters)) %>% # select only records with precision > 5 km
    
    # For plants it's good to use both preserved speciens and human observations. The former is often better than the latter, but to increase the number
    # of observations, it's smart to have both.
    subset(basisOfRecord == "HUMAN_OBSERVATION") %>% # delete occurrences based on fossil material, germplasm, literature
    
    # Do you have any good arguments for deleting occurrences before 1980? This threshold should be argued for
    subset(year >= 1980) %>% # delete occurrences from before 1980
    
    
    dplyr::select("species", "decimalLongitude", "decimalLatitude") %>% # select only relevant columns
    
    CoordinateCleaner::cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude", species = "species", additions = NULL, verbose = TRUE) %>% # remove duplicates
    CoordinateCleaner::cc_equ(lon = "decimalLongitude", lat = "decimalLatitude", test = "identical", verbose = TRUE) %>%
    CoordinateCleaner::cc_cap(lon = "decimalLongitude", lat = "decimalLatitude", verbose = TRUE) %>% # capital centroid detection
    CoordinateCleaner::cc_cen(lon = "decimalLongitude", lat = "decimalLatitude", test = "country", verbose = TRUE) %>% # country centroid detection
    CoordinateCleaner::cc_inst(lon = "decimalLongitude", lat = "decimalLatitude", verbose = TRUE) # check hyper-anthropogenic environment
  
  # transform to spatial data, retain one point per grid cell, transform to epsg:3035 and save to .csv
  species_of_interest_spatial <- 
    SpatialPointsDataFrame(coords = species_of_interest[, c("decimalLongitude", "decimalLatitude")], 
                           data = species_of_interest, proj4string = CRS("+init=epsg:4326")) %>%
    gridSample(BIO1, n = 1) %>% # one point per grid cell
    as.data.frame() %>%
    SpatialPointsDataFrame(.[, c("decimalLongitude", "decimalLatitude")], data = ., proj4string = CRS("+init=epsg:4326")) %>%
    spTransform(CRSobj = CRS("+init=epsg:3035")) %>% # transform to epsg:3035
    .@coords %>%
    as.data.frame() %>%
    mutate(species = species_names[q]) %>%
    write.csv(paste0(occ_clean_dir,"/", species_names[q], ".csv", sep = ""), row.names = FALSE)
}

# how many rows left?
nrow(species_of_interest)

##########################################################
######          Get environmental predictors         #####
##########################################################

# make a rasterstack of the environmental predictors
setwd(dir = env_dir)
env_stack <- stack(list.files())

# check the names
names(env_stack)

# convert to a dataframe
env_data <- as.data.frame(env_stack)

# remove NA values
env_data <- na.omit(env_data)

# check correlations between variables
corrplot.mixed(
  cor(env_data, method = c("spearman")),
  lower.col = "black",
  number.cex = 0.6,
  tl.pos = "lt",
  tl.col = "black",
  tl.srt = 45,
  cl.cex = 0.8
)

