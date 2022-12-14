#Calculate mean, median, max, as relevant for the different layers
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
              "CMIP6_MRI-ESM2-0")

#list of scenarios
sce_list <- c("historical" , "ssp245", "ssp585")

#periods
period_list <- c("hist", "near", "mid")

#indices
indx_list <- c("NDD", "NTx40", "NTx35", "HSM_NTx35", "HSH", "NDWS", "NDWL50", "THI", "TAI")

#load all masks
r_ref <- terra::rast(paste0(wd, "/roi/africa.tif")) # Africa base mask
r_wbd <- terra::rast(paste0(wd, "/roi/africa_wbd.tif")) # Africa water bodies
r_crop <- terra::rast(paste0(wd, "/mapspam_mask/spam2017V2r1_allcrop_mask.tif"))
r_lstk <- terra::rast(paste0(wd, "/livestock_mask/livestock_allspecies_mask.tif"))
r_pop <- terra::rast(paste0(wd, "/population_mask/pop_mask.tif"))

## Generate final continuous outputs
#domean=TRUE will calculate the climatological mean of each month, and then the mean of the 12 months
#domedian=TRUE will calculate the climatological median of each month, and then the median of the 12 months
#domax=TRUE will calculate the climatological mean of each month, and then the max of all months
#doensemble=TRUE will take all individual GCMs (in the case of rcp scenarios) and calculate the multi-model mean
#omitcalendar=TRUE will omit the use of the maize crop calendar when calculating the annual mean
continuous_map <- function(index="NDD", period="hist", scenario="historical", domean=TRUE, domedian=TRUE, domax=TRUE, doensemble=TRUE, omitcalendar=FALSE) {
  #index="NDD"; period="hist"; scenario="historical"; domean=TRUE; domedian=TRUE; domax=TRUE; doensemble=TRUE; omitcalendar=FALSE
  
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
    dir_names <- paste0(scenario, "_", gsub("CMIP6_", "", gcm_list), "_", min(years), "_", max(years))
  }
  
  #loop through the directories of interest
  for (i in 1:length(dir_names)) {
    #i <- 1
    #output file
    out_dir <- paste0(wd, "/cmip6/indices/", dir_names[i], "/", index, "/long_term_stats")
    if (!file.exists(out_dir)) {dir.create(out_dir, recursive=TRUE)}
    
    #input directory
    in_dir <- paste0(wd, "/cmip6/indices/", dir_names[i], "/", index)
    
    #list files
    fls <- list.files(path = in_dir, pattern = "\\.tif", full.names = TRUE)
    fls <- fls[grep(pattern = '.tif$', x = fls)]
    
    #names and indices for tapp
    bname <- basename(fls) %>% 
      gsub("NDD-", "", .) %>% 
      gsub("\\.tif", "", .)
    sumby <- bname %>%
      substr(., 6, 7) %>%
      as.numeric(.)
    
    #load all rasters
    r_data <- terra::rast(fls)
    names(r_data) <- bname
    
    #for the below use tapp with index being the month so that terra operates the layers by month
    #apply the various masks to these: r_ref, r_wbd, r_crop, r_lstk, r_pop
    #standard (all layers) apply r_ref, w_wbd
    #for THI also apply r_lstk
    #for HSM_NTx35 also apply r_crop
    #for NTx35 also apply r_crop
    #for HSH also apply r_pop
    
    #mean
    if (domean) {
      r_month <- tapp(r_data, sumby, fun=mean, na.rm=TRUE, filename=paste0(out_dir, "/mean_monthly.tif"), overwrite=TRUE)
      
      #for maize heat stress index that uses growing season, divide by growing season length instead of normal mean
      if (index == "HSM_NTx35" & !omitcalendar) {
        r_cal <- terra::rast("~/common_data/atlas_crop_calendar/intermediate/mai_rf_ggcmi_crop_calendar_phase3_v1.01_Africa.tif")
        r_cal <- r_cal[[3]] / 30
        r_mean <- sum(r_month, na.rm=TRUE) / r_cal
      } else {
        r_mean <- mean(r_month, na.rm=TRUE, filename=paste0(out_dir, "/mean_year.tif"), overwrite=TRUE)
      }
      
      r_monthk <- r_month %>%
        terra::mask(., r_ref) %>%
        terra::mask(., r_wbd, inverse=TRUE)
      if (index == "THI") {r_monthk <- terra::mask(r_monthk, r_lstk)}
      if (index %in% c("HSM_NTx35", "NTx35")) {r_monthk <- terra::mask(r_monthk, r_crop)}
      if (index == "HSH") {r_monthk <- terra::mask(r_monthk, r_pop)}
      terra::writeRaster(r_monthk, filename=paste0(out_dir, "/mean_monthly_masked.tif"), overwrite=TRUE)
      
      r_meank <- r_mean %>%
        terra::mask(., r_ref) %>%
        terra::mask(., r_wbd, inverse=TRUE)
      if (index == "THI") {r_meank <- terra::mask(r_meank, r_lstk)}
      if (index %in% c("HSM_NTx35", "NTx35")) {r_meank <- terra::mask(r_meank, r_crop)}
      if (index == "HSH") {r_meank <- terra::mask(r_meank, r_pop)}
      terra::writeRaster(r_meank, filename=paste0(out_dir, "/mean_year_masked.tif"), overwrite=TRUE)
    }
    
    #median
    if (domedian) {
      r_emonth <- tapp(r_data, sumby, fun=median, na.rm=TRUE, filename=paste0(out_dir, "/median_monthly.tif"), overwrite=TRUE)
      r_median <- median(r_emonth, na.rm=TRUE, filename=paste0(out_dir, "/median_year.tif"), overwrite=TRUE)
      
      r_emonthk <- r_emonth %>%
        terra::mask(., r_ref) %>%
        terra::mask(., r_wbd, inverse=TRUE)
      if (index == "THI") {r_emonthk <- terra::mask(r_emonthk, r_lstk)}
      if (index %in% c("HSM_NTx35", "NTx35")) {r_emonthk <- terra::mask(r_emonthk, r_crop)}
      if (index == "HSH") {r_emonthk <- terra::mask(r_emonthk, r_pop)}
      terra::writeRaster(r_emonthk, filename=paste0(out_dir, "/median_monthly_masked.tif"), overwrite=TRUE)
      
      r_mediank <- r_median %>%
        terra::mask(., r_ref) %>%
        terra::mask(., r_wbd, inverse=TRUE)
      if (index == "THI") {r_mediank <- terra::mask(r_mediank, r_lstk)}
      if (index %in% c("HSM_NTx35", "NTx35")) {r_mediank <- terra::mask(r_mediank, r_crop)}
      if (index == "HSH") {r_mediank <- terra::mask(r_mediank, r_pop)}
      terra::writeRaster(r_mediank, filename=paste0(out_dir, "/median_year_masked.tif"), overwrite=TRUE)
    }
    
    #max
    if (domax) {
      if (!file.exists(paste0(out_dir, "/mean_monthly.tif"))) {
        r_xmonth <- tapp(r_data, sumby, fun=mean, na.rm=TRUE, filename=paste0(out_dir, "/mean_monthly.tif"), overwrite=TRUE)
        r_xmonthk <- r_xmonth %>%
          terra::mask(., r_ref) %>%
          terra::mask(., r_wbd, inverse=TRUE)
        if (index == "THI") {r_xmonthk <- terra::mask(r_xmonthk, r_lstk)}
        if (index %in% c("HSM_NTx35", "NTx35")) {r_xmonthk <- terra::mask(r_xmonthk, r_crop)}
        if (index == "HSH") {r_xmonthk <- terra::mask(r_xmonthk, r_pop)}
        terra::writeRaster(r_xmonthk, filename=paste0(out_dir, "/mean_monthly_masked.tif"), overwrite=TRUE)
      } else {
        r_xmonth <- terra::rast(paste0(out_dir, "/mean_monthly.tif"))
      }
      r_max <- max(r_xmonth, na.rm=TRUE, filename=paste0(out_dir, "/max_year.tif"), overwrite=TRUE)
      
      r_maxk <- r_max %>%
        terra::mask(., r_ref) %>%
        terra::mask(., r_wbd, inverse=TRUE)
      if (index == "THI") {r_maxk <- terra::mask(r_maxk, r_lstk)}
      if (index %in% c("HSM_NTx35", "NTx35")) {r_maxk <- terra::mask(r_maxk, r_crop)}
      if (index == "HSH") {r_maxk <- terra::mask(r_maxk, r_pop)}
      terra::writeRaster(r_maxk, filename=paste0(out_dir, "/max_year_masked.tif"), overwrite=TRUE)
    }
  }
  return('Done\n')
}



