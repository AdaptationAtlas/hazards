#Calculate discrete maps
#HA/JRV, Dec 2022

# R options
rm(list = ls()) # Remove objects
g <- gc(reset = T); rm(g) # Empty garbage collector

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
category_map <- function(index="NDD", HS.stat=NULL, period="hist", scenario="historical", gcm="CMIP6_ENSEMBLE", stat="max_year"){
  #index="TAI"; HS.stat=NULL; period="near"; scenario="ssp245"; gcm="CMIP6_ENSEMBLE"; stat="mean_year"
  
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
  if (gcm == "CMIP6_ENSEMBLE" & scenario != "historical") {
    if (index %in% c("HSH", "THI")) {
      infile <- paste0(wd, "/cmip6/indices/", dir_names, "/", index, "_", HS.stat, "/", stat, c(".tif", "_masked.tif"))
      outfile <- paste0(wd, "/cmip6/indices/", dir_names, "/", index, "_", HS.stat, "/", stat, c("", "_masked"), "_categorical.tif")
    } else {
      infile <- paste0(wd, "/cmip6/indices/", dir_names, "/", index, "/", stat, c(".tif", "_masked.tif"))
      outfile <- paste0(wd, "/cmip6/indices/", dir_names, "/", index, "/", stat, c("", "_masked"), "_categorical.tif")
    }
  } else {
    if (index %in% c("HSH", "THI")) {
      infile <- paste0(wd, "/cmip6/indices/", dir_names, "/", index, "/long_term_stats_", HS.stat, "/", stat, c(".tif", "_masked.tif"))
      outfile <- paste0(wd, "/cmip6/indices/", dir_names, "/", index, "/long_term_stats_", HS.stat, "/", stat, c("", "_masked"), "_categorical.tif")
    } else {
      infile <- paste0(wd, "/cmip6/indices/", dir_names, "/", index, "/long_term_stats/", stat, c(".tif", "_masked.tif"))
      outfile <- paste0(wd, "/cmip6/indices/", dir_names, "/", index, "/long_term_stats/", stat, c("", "_masked"), "_categorical.tif")
    }
  }
  
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

#selected index
indx <- indx_list[6]
cat("processing index=", indx, "\n")

#historical
for (stt in stat_list) {
  if (indx %in% c("HSH", "THI")) {
    category_map(index=indx, HS.stat="max", period=period_list[1], scenario=sce_list[1], gcm="CMIP6_ENSEMBLE", stat=stt)
    category_map(index=indx, HS.stat="mean", period=period_list[1], scenario=sce_list[1], gcm="CMIP6_ENSEMBLE", stat=stt)
  } else {
    if (indx == "TAI") {
      if (stt != "max_year") {category_map(index=indx, HS.stat=NULL, period=period_list[1], scenario=sce_list[1], gcm="CMIP6_ENSEMBLE", stat=stt)}
    } else {
      category_map(index=indx, HS.stat=NULL, period=period_list[1], scenario=sce_list[1], gcm="CMIP6_ENSEMBLE", stat=stt)
    }
  }
}

#future
for (i in 1:length(gcm_list)) {
  cat(gcm_list[i], "\n")
  for (prd in period_list[2:3]) {
    for (sce in sce_list[2:3]) {
      for (stt in stat_list) {
        if (indx %in% c("HSH", "THI")) {
          category_map(index=indx, HS.stat="max", period=prd, scenario=sce, gcm=gcm_list[i], stat=stt)
          category_map(index=indx, HS.stat="mean", period=prd, scenario=sce, gcm=gcm_list[i], stat=stt)
        } else {
          if (indx == "TAI") {
            if (stt != "max_year") {category_map(index=indx, HS.stat=NULL, period=prd, scenario=sce, gcm=gcm_list[i], stat=stt)}
          } else {
            category_map(index=indx, HS.stat=NULL, period=prd, scenario=sce, gcm=gcm_list[i], stat=stt)
          }
        }
      }
    }
  }
}
