#Calculate discrete maps
#HA/JRV, Dec 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
#options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation

#load packages
library(terra)
library(tidyverse)

#working directory
wd <- "~/common_data/atlas_hazards"

#list of GCMs
gcm_list <- c("CMIP6_ACCESS-ESM1-5",
              "CMIP6_MPI-ESM1-2-HR",
              "CMIP6_EC-Earth3",
              "CMIP6_INM-CM5-0",
              "CMIP6_MRI-ESM2-0",
              "CMIP6_ENSEMBLE")

#list of scenarios
sce_list <- c("historical" , "ssp245", "ssp585")

#periods
period_list <- c("hist", "near", "mid")

#indices
indx_list <- c("NDD", "NTx40", "NTx35", "HSM_NTx35", "HSH", "NDWS", "NDWL50", "THI", "TAI")

#stat
stat_list <- c("mean_year", "max_year", "median_year")

#source make_class_tb() function
source("~/Repositories/hazards/R/05_final_maps/makeClassTable.R")

#generate final categorical maps
category_map <- function(index="NDD", period="hist", scenario="historical", gcm="CMIP6_ENSEMBLE", stat="max_year"){
  #index="NDD"; period="hist"; scenario="historical"; gcm="CMIP6_ENSEMBLE"; stat="mean_year"
  
  #create class table
  class_tb <- make_class_tb()
  
  #some consistency checks
  if (scenario == "historical" & period != "hist") {stop("error, historical scenario should go with hist period")}
  if (scenario %in% c("ssp245", "ssp585") & !period %in% c("near", "mid")) {stop("error, use near or mid, for ssp scenarios")}
  
  #assign periods
  if (period == "hist") {years <- 1995:2014}
  if (period == "near") {years <- 2021:2040}
  if (period == "mid") {years <- 2041:2060}
  
  #list of folders to go through (1 for hist, 5 GCMs for future)
  #note folder structure
  #historical
  #ssp245_<GCM>_2021_2040
  #ssp245_<GCM>_2041_2060
  #ssp585_<GCM>_2021_2040
  #ssp585_<GCM>_2041_2060
  if (scenario == "historical") {
    dir_names <- scenario
  } else {
    dir_names <- paste0(scenario, "_", gsub("CMIP6_", "", gcm_list[grep(gcm, gcm_list)]), "_", min(years), "_", max(years))
  }
  
  #input file, both unmasked and masked
  infile <- paste0(wd, "/cmip6/indices/", dir_names, "/", index, "/", stat, c(".tif", "_masked.tif"))
  outfile <- paste0(wd, "/cmip6/indices/", dir_names, "/", index, "/", stat, c("", "_masked"), "_categorical.tif")
  
  #load raster
  r_in <- terra::rast(infile)
  
  #get category matrix for reclassification from class_tb
  cat_m <- class_tb %>%
    dplyr::filter(index_name == index) %>%
    dplyr::select(lower_lim:description)
  cls <- cat_m %>%
    dplyr::select(class, description) %>%
    dplyr::rename(id=class, level=description)
  cat_m <- cat_m %>%
    dplyr::select(-description) %>%
    as.matrix()
  
  #reclassify raster
  for (i in 1:terra::nlyr(r_in)) {
    rc <- terra::classify(x = r_in[[i]], rcl = cat_m, right = FALSE)
    levels(rc) <- cls
    terra::writeRaster(rc, outfile[i], overwrite = TRUE)
  }
  return("Done\n")
}

