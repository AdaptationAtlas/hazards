#Calculate heat stress maize
#JRV, Dec. 2022

#load packages
library(terra)
library(tidyverse)
library(lubridate)

#clean-up environment
rm(list=ls())
gc(verbose=FALSE, full=TRUE, reset=TRUE)

#working directory
wd <- "~/common_data"

#reference raster
r_ref <- terra::rast(paste0(wd,'/atlas_hazards/roi/africa.tif'))

# Calculate NTx40 function
calc_hsm <- function(yr, mn, thr, allyear=FALSE){
  #yr <- 1995
  #mn <- '02'
  #thr <- 35
  #give a file name
  if (allyear) {
    outfile <- paste0(out_dir, "/HSM_NTx", thr, "/AllYear_HSM_NTx", thr, "-", yr, "-", mn, ".tif")
  } else {
    outfile <- paste0(out_dir, "/HSM_NTx", thr, "/GSeason_HSM_NTx", thr, "-", yr, "-", mn, ".tif")
  }
  cat(outfile, "\n")
  
  #process only if file exists
  if (!file.exists(outfile)) {
    #if directory doesnt exist, create it
    if (!file.exists(dirname(outfile))) dir.create(dirname(outfile), recursive=TRUE)
    
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    
    # Files
    fls <- paste0(tx_dir, '/', yr, '/Tmax.', gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    fls <- fls[file.exists(fls)]
    
    # Read raster data
    tmx <- terra::rast(fls)
    tmx <- tmx %>% terra::crop(terra::ext(r_ref)) %>% terra::mask(r_ref)
    tmx <- terra::app(x=tmx, fun=function(x) {rout <- x; rout[rout[] < -9990] <- NA; return(rout)})
    
    #get dates back from names
    tdates <- names(tmx)
    tdates <- gsub(pattern="Tmax.", replacement="", x=tdates)
    tdates <- gsub(pattern=".", replacement="-", x=tdates, fixed=TRUE)
    names(tmx) <- tdates
    
    #calculate stack with day number
    tmx_doy <- tmx
    for (.x in 1:terra::nlyr(tmx)) {
      tmx_doy[[.x]][!is.na(tmx[[.x]][])] <- yday(as.Date(names(tmx[[.x]])))
    }
    
    # Read crop calendar data, if needed. Otherwise generate an all-year calendar
    if (allyear) {
      rc1 <- rc2 <- terra::rast(r_ref)
      rc1[!is.na(r_ref[])] <- 1
      rc2[!is.na(r_ref[])] <- 366
      r_cal <- c(rc1, rc2)
      rm(list=c("rc1", "rc2"))
      gc(verbose=FALSE, full=TRUE, reset=TRUE)
    } else {
      r_cal <- terra::rast(paste0(wd, "/atlas_crop_calendar/intermediate/mai_rf_ggcmi_crop_calendar_phase3_v1.01_Africa.tif"))
    }
    
    #calculate heat stress day yes/no
    r_hsm <- c()
    for (.x in 1:terra::nlyr(tmx)) {
      #determine if pixels exceed threshold
      rx <- terra::app(tmx[[.x]], fun=function(x) {ifelse(x >= thr, 1, 0)})
      
      #append layers together for lapp
      rx <- c(rx, tmx_doy[[.x]], r_cal[[1]], r_cal[[2]])
      
      #now determine if that day is in season
      ry <- terra::lapp(rx, fun=function(tt, doy, pd, md) {
        tout <- ifelse(md > pd, #normal season
                       ifelse(doy >= pd & doy <= md, #we are in season
                              ifelse(tt == 1, 1, 0), 0), #it is a hot day
                       ifelse(md < pd, #season starts and ends in diff. years
                              ifelse(doy <= md | doy >= pd, #we are not in season
                                     ifelse(tt == 1, 1, 0), 0), 0)) #it is a hot day
        return(tout)
      })
      
      #append
      r_hsm <- c(r_hsm, ry)
      
      #clean up
      rm(rx)
      gc(verbose=FALSE, full=TRUE, reset=TRUE)
    }
    r_hsm <- terra::rast(r_hsm)
    
    # Calculate heat stress generic crop
    terra::app(x   = r_hsm,
               fun = sum, na.rm=TRUE,
               filename = outfile)
    
    #clean up
    rm(r_hsm, ry, tmx, tmx_doy, tdates, r_cal)
    gc(verbose=FALSE, full=TRUE, reset=TRUE)
  }
}

# Setup
#scenario <- "ssp245" #c("historical", "ssp245", "ssp585")
#period <- "near" #c("hist", "near", "mid")
#gcm <- "MPI-ESM1-2-HR" #"ACCESS-ESM1-5","MPI-ESM1-2-HR", "EC-Earth3", "INM-CM5-0", "MRI-ESM2-0")

for (gcm in c("ACCESS-ESM1-5", "MPI-ESM1-2-HR", "EC-Earth3", "INM-CM5-0", "MRI-ESM2-0")) {
    for (scenario in c("ssp245", "ssp585")) {
        for (period in c("near", "mid")) {
            #assign periods
            if (period == "hist") {yrs <- 1995:2014}
            if (period == "near") {yrs <- 2021:2040}
            if (period == "mid") {yrs <- 2041:2060}

            #out_dir
            if (scenario == "historical") {
              out_dir <- paste0(wd, "/atlas_hazards/cmip6/indices/historical")
            } else {
              out_dir <- paste0(wd, "/atlas_hazards/cmip6/indices/", scenario, "_", gcm, "_", min(yrs), "_", max(yrs))
            }

            #chirts_dir
            if (scenario == "historical") {
              tx_dir <- paste0(wd, "/chirts/Tmax")
            } else {
              tx_dir <- paste0(wd, "/chirts_cmip6_africa/Tmax_", gcm, "_", scenario, "_", min(yrs), "_", max(yrs))
            }

            #data.frame for processing
            mns <- c(paste0('0',1:9),10:12)
            stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
            names(stp) <- c('yrs','mns')
            stp <- stp %>%
              dplyr::arrange(yrs, mns) %>%
              base::as.data.frame()

            1:nrow(stp) %>%
              purrr::map(.f = function(i){calc_hsm(yr = stp$yrs[i], 
                                                   mn = stp$mns[i],
                                                   thr = 35,
                                                   allyear = FALSE)})
        }
    }
}
