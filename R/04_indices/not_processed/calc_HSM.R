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

# Calculate HSM_NTx function
calc_hsm <- function(yr, mn, thr=35, allyear=FALSE){
  #yr <- 1995
  #mn <- '02'
  #thr <- 35
  #give a file name
  if (allyear) {
    outfile <- paste0(out_dir, "/HSM_NTx", thr, "/AllYear_HSM_NTx", thr, "-", yr, "-", mn, ".tif")
  } else {
    outfile <- paste0(out_dir, "/HSM_NTx", thr, "/GSeason_HSM_NTx", thr, "-", yr, "-", mn, ".tif")
  }
  thr <- thr[!file.exists(outfile)]
  outfile <- outfile[!file.exists(outfile)]
  cat("...processing n=", length(outfile), "files for yr=", yr, "/ mn=", mn, "\n")
  cat("...", out_dir, "\n")
  
  #process only if file exists
  if (length(outfile) > 0) {
    #if directory doesnt exist, create it
    1:length(outfile) %>%
      purrr::map(.f = function(j){dir.create(dirname(outfile[j]),F,T)})
    
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
    
    for (j in 1:length(thr)) {
      cat("...processing threshold thr=", thr[j], "/", out_dir, "\n")
      #calculate heat stress day yes/no
      r_hsm <- c()
      for (.x in 1:terra::nlyr(tmx)) {
        #determine if pixels exceed threshold
        rx <- terra::app(tmx[[.x]], fun=function(x) {ifelse(x >= thr[j], 1, 0)})
        
        #append layers together for lapp
        rx <- c(rx, tmx_doy[[.x]], r_cal[[1]], r_cal[[2]])
        
        #now determine if that day is in season
        ry <- terra::lapp(rx, fun=function(tt, doy, pd, md) {
          tout <- ifelse(md > pd, #normal season
                         ifelse(doy >= pd & doy <= md, #we are in season
                                ifelse(tt == 1, 1, 0), 0), #it is a hot day
                         ifelse(md < pd, #season starts and ends in diff. years
                                ifelse(doy <= md | doy >= pd, #we are in season
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
                 filename = outfile[j])
      
      #clean up
      rm(r_hsm, ry)
      gc(verbose=FALSE, full=TRUE, reset=TRUE)
    }
    rm(tmx, tmx_doy, tdates, r_cal)
    gc(verbose=FALSE, full=TRUE, reset=TRUE)
  }
}

# Setup
scenario <- "historical" #c("historical", "ssp245", "ssp585")
period <- "hist" #c("hist", "near", "mid")
gcm <- "ACCESS-ESM1-5" #"ACCESS-ESM1-5","MPI-ESM1-2-HR", "EC-Earth3", "INM-CM5-0", "MRI-ESM2-0")

#for (scenario in c("historical", "ssp245", "ssp585")) {
  if (scenario == "historical") {gcm_list <- "historical"} else {gcm_list <- c("ACCESS-ESM1-5", "MPI-ESM1-2-HR", "EC-Earth3", "INM-CM5-0", "MRI-ESM2-0")}
#  for (gcm in gcm_list) {
      if (scenario == "historical") {per_list <- "hist"} else {per_list <- c("near", "mid")}
#      for (period in per_list) {
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
                                                   thr = 30:50,
                                                   allyear = FALSE)})
#        }
#    }
#}


# # ----------------------------------------------------------------------
# # Data fixes
# # Get reruns file.
# source("~/Repositories/hazards/R/03_bias_correction/getReruns.R")
# other_bfiles <- c(paste0(wd, "/chirps_cmip6_africa/Prec_ACCESS-ESM1-5_ssp245_2041_2060/chirps-v2.0.2043.03.18.tif"), 
#                   paste0(wd, "/chirps_cmip6_africa/Prec_ACCESS-ESM1-5_ssp245_2041_2060/chirps-v2.0.2054.10.10.tif"),
#                   paste0(wd, "/chirps_cmip6_africa/Prec_ACCESS-ESM1-5_ssp245_2041_2060/chirps-v2.0.2056.03.29.tif"))
# reruns_df <- getReruns(newfiles=other_bfiles) %>%
#   dplyr::filter(varname == "Tmax") %>%
#   unique(.)
# 
# # Do the reruns
# for (j in 1:nrow(reruns_df)) {
#   gcm <- reruns_df$gcm[j]
#   ssp <- reruns_df$ssp[j]
#   prd <- reruns_df$prd[j]
#   cmb <- paste0(ssp,'_',gcm,'_',prd)
#   prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
#   mns <- c(paste0('0',1:9),10:12)
#   yrs <- prd_num[1]:prd_num[2]
#   stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
#   names(stp) <- c('yrs','mns')
#   
#   #folders
#   tx_dir <- paste0(wd, "/chirts_cmip6_africa/Tmax_", gcm, "_", ssp, "_", prd)
#   out_dir <- paste0(wd, "/atlas_hazards/cmip6/indices/", ssp, "_", gcm, "_", prd)
#   
#   #remove and redo file
#   this_file <- paste0(out_dir, "/HSM_NTx35/GSeason_HSM_NTx35-", reruns_df$yr[j], "-", reruns_df$mn[j], ".tif")
#   cat("redoing file", this_file, "\n")
#   system(paste0("rm -f ", this_file))
#   calc_hsm(yr = reruns_df$yr[j], 
#            mn = reruns_df$mn[j],
#            thr = 35,
#            allyear = FALSE)
# }

