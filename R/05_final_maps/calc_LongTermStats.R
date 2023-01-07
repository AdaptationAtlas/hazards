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
  #index="TAI"; period="hist"; scenario="historical"; domean=TRUE; domedian=TRUE; domax=TRUE; doensemble=TRUE; omitcalendar=FALSE
  
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
  
  #if doensemble then create empty vectors to store the rasters
  if (doensemble & scenario %in% c("ssp245", "ssp585")) {
    if (domean) {ens_meanmonth <- ens_meanmonthk <- ens_meanyr <- ens_meanyrk <- c()}
    if (domedian) {ens_emonth <- ens_emonthk <- ens_medianyr <- ens_medianyrk <- c()}
    if (domax) {ens_xmonth <- ens_xmonthk <- ens_max <- ens_maxk <- c()}
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
    #get YEAR-MONTH name structure
    bname <- basename(fls) %>% 
      gsub(paste0(index,"-"), "", .) %>% 
      gsub("\\.tif", "", .)
    
    #all indices except TAI are done by month, hence requiring long-term monthly means
    #this extracts the month, which is used below in terra::tapp
    if (index != "TAI") {
      sumby <- bname %>%
        substr(., 6, 7) %>%
        as.numeric(.)
    }
    
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
      if (index != "TAI") {
        #calculate climatological mean of each month
        r_month <- tapp(r_data, sumby, fun=mean, na.rm=TRUE, filename=paste0(out_dir, "/mean_monthly.tif"), overwrite=TRUE)
        
        #mask the climatological mean
        r_monthk <- r_month %>%
          terra::mask(., r_ref) %>%
          terra::mask(., r_wbd, inverse=TRUE)
        if (index == "THI") {r_monthk <- terra::mask(r_monthk, r_lstk)}
        if (index %in% c("HSM_NTx35", "NTx35")) {r_monthk <- terra::mask(r_monthk, r_crop)}
        if (index == "HSH") {r_monthk <- terra::mask(r_monthk, r_pop)}
        terra::writeRaster(r_monthk, filename=paste0(out_dir, "/mean_monthly_masked.tif"), overwrite=TRUE)
        
        #calculate mean of all months (average of climatological means)
        #for maize heat stress index that uses growing season, divide by growing season length instead of normal mean
        if (index == "HSM_NTx35" & !omitcalendar) {
          r_cal <- terra::rast("~/common_data/atlas_crop_calendar/intermediate/mai_rf_ggcmi_crop_calendar_phase3_v1.01_Africa.tif")
          r_cal <- r_cal[[3]] / 30
          r_mean <- sum(r_month, na.rm=TRUE) / r_cal
        } else {
          r_mean <- mean(r_month, na.rm=TRUE, filename=paste0(out_dir, "/mean_year.tif"), overwrite=TRUE)
        }
      } else {
        r_mean <- mean(r_data, na.rm=TRUE, filename=paste0(out_dir, "/mean_year.tif"), overwrite=TRUE)
      }
      
      #mask the annual mean
      r_meank <- r_mean %>%
        terra::mask(., r_ref) %>%
        terra::mask(., r_wbd, inverse=TRUE)
      if (index == "THI") {r_meank <- terra::mask(r_meank, r_lstk)}
      if (index %in% c("HSM_NTx35", "NTx35")) {r_meank <- terra::mask(r_meank, r_crop)}
      if (index == "HSH") {r_meank <- terra::mask(r_meank, r_pop)}
      terra::writeRaster(r_meank, filename=paste0(out_dir, "/mean_year_masked.tif"), overwrite=TRUE)
      
      #do ensemble
      if (doensemble & scenario %in% c("ssp245", "ssp585")) {
        if (index != "TAI") {
          ens_meanmonth <- c(ens_meanmonth, r_month) 
          ens_meanmonthk <- c(ens_meanmonthk, r_monthk)
        }
        ens_meanyr <- c(ens_meanyr, r_mean)
        ens_meanyrk <- c(ens_meanyrk, r_meank)
      }
    }
    
    #median
    if (domedian) {
      if (index != "TAI") {
        #climatological median of each month
        r_emonth <- tapp(r_data, sumby, fun=median, na.rm=TRUE, filename=paste0(out_dir, "/median_monthly.tif"), overwrite=TRUE)
        
        #annual median of all months
        r_median <- median(r_emonth, na.rm=TRUE, filename=paste0(out_dir, "/median_year.tif"), overwrite=TRUE)
        
        #mask climatological monthly median
        r_emonthk <- r_emonth %>%
          terra::mask(., r_ref) %>%
          terra::mask(., r_wbd, inverse=TRUE)
        if (index == "THI") {r_emonthk <- terra::mask(r_emonthk, r_lstk)}
        if (index %in% c("HSM_NTx35", "NTx35")) {r_emonthk <- terra::mask(r_emonthk, r_crop)}
        if (index == "HSH") {r_emonthk <- terra::mask(r_emonthk, r_pop)}
        terra::writeRaster(r_emonthk, filename=paste0(out_dir, "/median_monthly_masked.tif"), overwrite=TRUE)
      } else {
        #climatological median of annual quantity (for TAI)
        r_median <- median(r_data, na.rm=TRUE, filename=paste0(out_dir, "/median_year.tif"), overwrite=TRUE)
      }
      
      r_mediank <- r_median %>%
        terra::mask(., r_ref) %>%
        terra::mask(., r_wbd, inverse=TRUE)
      if (index == "THI") {r_mediank <- terra::mask(r_mediank, r_lstk)}
      if (index %in% c("HSM_NTx35", "NTx35")) {r_mediank <- terra::mask(r_mediank, r_crop)}
      if (index == "HSH") {r_mediank <- terra::mask(r_mediank, r_pop)}
      terra::writeRaster(r_mediank, filename=paste0(out_dir, "/median_year_masked.tif"), overwrite=TRUE)
      
      #do ensemble
      if (doensemble & scenario %in% c("ssp245", "ssp585")) {
        if (index != "TAI") {
          ens_emonth <- c(ens_emonth, r_emonth)
          ens_emonthk <- c(ens_emonthk, r_emonthk)
        }
        ens_medianyr <- c(ens_medianyr, r_median)
        ens_medianyrk <- c(ens_medianyrk, r_mediank)
      }
    }
    
    #max. Note this statistic cannot be calculated for TAI since it does not have monthly values
    if (domax & index != "TAI") {
      if (!file.exists(paste0(out_dir, "/mean_monthly.tif"))) {
        #climatological mean of monthly values
        r_xmonth <- tapp(r_data, sumby, fun=mean, na.rm=TRUE, filename=paste0(out_dir, "/mean_monthly.tif"), overwrite=TRUE)
        
        #mask the climatological mean
        r_xmonthk <- r_xmonth %>%
          terra::mask(., r_ref) %>%
          terra::mask(., r_wbd, inverse=TRUE)
        if (index == "THI") {r_xmonthk <- terra::mask(r_xmonthk, r_lstk)}
        if (index %in% c("HSM_NTx35", "NTx35")) {r_xmonthk <- terra::mask(r_xmonthk, r_crop)}
        if (index == "HSH") {r_xmonthk <- terra::mask(r_xmonthk, r_pop)}
        terra::writeRaster(r_xmonthk, filename=paste0(out_dir, "/mean_monthly_masked.tif"), overwrite=TRUE)
      } else {
        r_xmonth <- terra::rast(paste0(out_dir, "/mean_monthly.tif"))
        r_xmonthk <- terra::rast(paste0(out_dir, "/mean_monthly_masked.tif"))
      }
      
      #calculate maximum of all months using the climatological monthly means
      r_max <- max(r_xmonth, na.rm=TRUE, filename=paste0(out_dir, "/max_year.tif"), overwrite=TRUE)
      
      #mask the annual maximum
      r_maxk <- r_max %>%
        terra::mask(., r_ref) %>%
        terra::mask(., r_wbd, inverse=TRUE)
      if (index == "THI") {r_maxk <- terra::mask(r_maxk, r_lstk)}
      if (index %in% c("HSM_NTx35", "NTx35")) {r_maxk <- terra::mask(r_maxk, r_crop)}
      if (index == "HSH") {r_maxk <- terra::mask(r_maxk, r_pop)}
      terra::writeRaster(r_maxk, filename=paste0(out_dir, "/max_year_masked.tif"), overwrite=TRUE)
      
      #do ensemble
      if (doensemble & scenario %in% c("ssp245", "ssp585")) {
        ens_xmonth <- c(ens_xmonth, r_xmonth)
        ens_xmonthk <- c(ens_xmonthk, r_xmonthk)
        ens_max <- c(ens_max, r_max)
        ens_maxk <- c(ens_maxk, r_maxk)
      }
    }
  }
  
  ###
  #calculate ensemble
  if (scenario %in% c("ssp245", "ssp585") & doensemble) {
    #output folder for ensemble calculation
    ens_dir <- paste0(wd, "/cmip6/indices/", scenario, "_ENSEMBLE_", min(years), "_", max(years), "/", index)
    if (!file.exists(ens_dir)) {dir.create(ens_dir, recursive=TRUE)}
    
    if (domean) {
      if (index != "TAI") {
        ens_meanmonth <- terra::rast(ens_meanmonth) %>%
          terra::tapp(., rep(1:12, 5), fun=mean, na.rm=TRUE, filename=paste0(ens_dir, "/mean_monthly.tif"), overwrite=TRUE)
        ens_meanmonthk <- terra::rast(ens_meanmonthk) %>%
          terra::tapp(., rep(1:12, 5), fun=mean, na.rm=TRUE, filename=paste0(ens_dir, "/mean_monthly_masked.tif"), overwrite=TRUE)
      }
      ens_meanyr <- terra::rast(ens_meanyr) %>%
        mean(., na.rm=TRUE, filename=paste0(ens_dir, "/mean_year.tif"), overwrite=TRUE)
      ens_meanyrk <- terra::rast(ens_meanyrk) %>%
        mean(., na.rm=TRUE, filename=paste0(ens_dir, "/mean_year_masked.tif"), overwrite=TRUE)
    }
    
    if (domedian) {
      if (index != "TAI") {
        ens_emonth <- terra::rast(ens_emonth) %>%
          terra::tapp(., rep(1:12, 5), fun=mean, na.rm=TRUE, filename=paste0(ens_dir, "/median_monthly.tif"), overwrite=TRUE)
        ens_emonthk <- terra::rast(ens_emonthk) %>%
          terra::tapp(., rep(1:12, 5), fun=mean, na.rm=TRUE, filename=paste0(ens_dir, "/median_monthly_masked.tif"), overwrite=TRUE)
      }
      ens_medianyr <- terra::rast(ens_medianyr) %>%
        mean(., na.rm=TRUE, filename=paste0(ens_dir, "/median_year.tif"), overwrite=TRUE)
      ens_medianyrk <- terra::rast(ens_medianyrk) %>%
        mean(., na.rm=TRUE, filename=paste0(ens_dir, "/median_year_masked.tif"), overwrite=TRUE)
    }
    
    if (domax & index != "TAI") {
      if (!domean) { #if domean we assume it has already been written
        ens_xmonth <- terra::rast(ens_xmonth) %>%
          terra::tapp(., rep(1:12, 5), fun=mean, na.rm=TRUE, filename=paste0(ens_dir, "/mean_monthly.tif"), overwrite=TRUE)
        ens_xmonthk <- terra::rast(ens_xmonthk) %>%
          terra::tapp(., rep(1:12, 5), fun=mean, na.rm=TRUE, filename=paste0(ens_dir, "/mean_monthly_masked.tif"), overwrite=TRUE)
      }
      ens_max <- terra::rast(ens_max) %>%
        mean(., na.rm=TRUE, filename=paste0(ens_dir, "/max_year.tif"), overwrite=TRUE)
      ens_maxk <- terra::rast(ens_maxk) %>%
        mean(., na.rm=TRUE, filename=paste0(ens_dir, "/max_year_masked.tif"), overwrite=TRUE)
    }
  }
  return("Done\n")
}



