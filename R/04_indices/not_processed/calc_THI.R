## Heat stress livestock (cattle) (THI)
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

sce_climate <- "historical" #"future"

# Calculate THI function
calc_thi <- function(yr, mn, dataset="CHIRPS"){
  outfile1 <- paste0(out_dir,'/THI_mean-',yr,'-',mn,'.tif')
  outfile2 <- paste0(out_dir,'/THI_max-',yr,'-',mn,'.tif')
  outfile3 <- paste0(out_dir,'/daily/THI_daily-',yr,'-',mn,'.tif')
  file.exists(c(outfile1,outfile2,outfile3))
  cat(outfile2, "\n")
  if(sum(file.exists(c(outfile1,outfile2,outfile3))) < 3){
    dir.create(dirname(outfile3),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    
    # Sequence of dates
    if(dataset == "CHIRPS" & as.numeric(yr) > 2020 & mn == '02'){
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-28')), by = 'day')
    } else {
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    }
    
    # Files
    if (dataset == "CHIRPS") {
      tx_fls <- paste0(tx_pth,'/',yr,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
      #this if makes sure historical RH is used for future scenarios. Should be changed once bias corrected RH exists
      if(as.numeric(yr) > 2020){
        dts_chr <- as.character(dts)
        dts_chr <- gsub(pattern = yr, replacement = yrs_mpg$Baseline[yrs_mpg$Future == yr], x = dts_chr)
        dts_h <- as.Date(dts_chr); rm(dts_chr)
        yr_h <- lubridate::year(dts_h)[1]
        rh_fls <- paste0(rh_pth,'/',yr_h,'/RH.',gsub(pattern='-', replacement='.', x=dts_h, fixed=T),'.tif')
      } else {
        rh_fls <- paste0(rh_pth,'/',yr,'/RH.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
      }
    } else if (dataset == "AgERA5") {
      #only for historical climate
      tx_fls <- paste0(ae5tx_pth,'/Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0.nc')
      rh_fls <- paste0(ae5rh_pth,'/Relative-Humidity-2m-12h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0.nc')
      if(as.numeric(yr) > 2020 & sce_climate == "future"){stop("dataset AgERA5 does not allow for future climate calculations\n")}
    }
    tx_fls <- tx_fls[file.exists(tx_fls)]
    rh_fls <- rh_fls[file.exists(rh_fls)]
    
    # Read variables
    tmx <- terra::rast(tx_fls)
    if (dataset == "AgERA5") {
      tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::resample(., ref) %>% terra::mask(ref)
      tmx <- tmx - 273.15
    } else {
      tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
      tmx[tmx == -9999] <- NA
    }
    
    rhm <- terra::rast(rh_fls)
    if (dataset == "AgERA5") {
      rhm <- rhm %>% terra::crop(terra::ext(ref)) %>% terra::resample(., ref) %>% terra::mask(ref)
    } else {
      rhm <- rhm %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    }
    
    thr_hum_idx <- function(tmax, rhum){
      thi = (1.8 * tmax + 32) - ((0.55 - 0.0055 * rhum) * (1.8 * tmax - 26.8))
      return(thi)
    }
    # Calculate human heat stress
    THI <- terra::lapp(x = terra::sds(tmx,rhm), fun = thr_hum_idx)
    terra::time(THI) <- dts
    THI_avg <- mean(THI, na.rm = T) %>% terra::mask(ref)
    THI_max <- max(THI, na.rm = T) %>% terra::mask(ref)
    
    # Write output
    terra::writeRaster(THI_avg, outfile1)
    terra::writeRaster(THI_max, outfile2)
    terra::writeRaster(THI, outfile3)
    
    # Clean up
    rm(tmx, rhm, THI, THI_avg, THI_max)
    gc(verbose=F, full=T, reset=T)
  }
}

if (sce_climate == "historical") {
  # Historical setup
  yrs <- 1981:2022 #1995:2014
  mns <- c(paste0('0',1:9),10:12)
  stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
  names(stp) <- c('yrs','mns')
  stp <- stp %>%
    dplyr::arrange(yrs, mns) %>%
    base::as.data.frame()
  
  #input data paths for chirps
  tx_pth <- paste0(root,'/chirts/Tmax') # Maximum temperature
  rh_pth <- paste0(root,'/chirts/RHum') # Relative humidity
  
  #input data paths for agera5
  ae5tx_pth <- paste0(root,'/ecmwf_agera5/2m_temperature-24_hour_maximum') # Maximum temperature
  ae5rh_pth <- paste0(root,'/ecmwf_agera5/2m_relative_humidity') # Relative humidity
  
  #output path
  out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/historical/THI_AgERA5')
  
  1:nrow(stp) %>%
    purrr::map(.f = function(i){
      calc_thi(yr = stp$yrs[i], mn = stp$mns[i], dataset="AgERA5")
      gc(verbose=F, full=T, reset=T)
      tmpfls <- list.files(tempdir(), full.names=TRUE)
      1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
    })
} else if (sce_climate == "future") {
  # Future setup
  #gcm <- 'ACCESS-ESM1-5'
  #ssp <- 'ssp245'
  #prd <- '2021_2040'
  
  for (gcm in c("ACCESS-ESM1-5", "MPI-ESM1-2-HR", "EC-Earth3", "INM-CM5-0", "MRI-ESM2-0")) {
      for (ssp in c('ssp245', 'ssp585')) {
          for (prd in c('2021_2040', '2041_2060')) {
              cmb <- paste0(ssp,'_',gcm,'_',prd)
              prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
              yrs <- prd_num[1]:prd_num[2]
              mns <- c(paste0('0',1:9),10:12)
              stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
              names(stp) <- c('yrs','mns')
              stp <- stp %>%
                dplyr::arrange(yrs, mns) %>%
                base::as.data.frame()
              tx_pth <- paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd) # Maximum temperature
              rh_pth <- paste0(root,'/chirts/RHum')                                   # Relative humidity
              out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/THI')
  
              yrs_mpg <- data.frame(Baseline = as.character(rep(1995:2014, 2)),
                                    Future = as.character(c(2021:2040,2041:2060)))
  
              1:nrow(stp) %>%
                purrr::map(.f = function(i){calc_thi(yr = stp$yrs[i], mn = stp$mns[i]); gc(verbose=F, full=T, reset=T)})
          }
      }
  }
} else {
  cat("select one of historical or future for sce_climate \n")
}



# # ----------------------------------------------------------------------
# # ----------------------------------------------------------------------
# # Data fixes
# # Get reruns file.
# source("~/Repositories/hazards/R/03_bias_correction/getReruns.R")
# other_bfiles <- c(paste0(root, "/chirps_cmip6_africa/Prec_ACCESS-ESM1-5_ssp245_2041_2060/chirps-v2.0.2043.03.18.tif"), 
#                   paste0(root, "/chirps_cmip6_africa/Prec_ACCESS-ESM1-5_ssp245_2041_2060/chirps-v2.0.2054.10.10.tif"),
#                   paste0(root, "/chirps_cmip6_africa/Prec_ACCESS-ESM1-5_ssp245_2041_2060/chirps-v2.0.2056.03.29.tif"))
# reruns_df <- getReruns(newfiles=other_bfiles) %>%
#   dplyr::filter(varname == "Tmax") %>%
#   unique(.)
# 
# # Do the reruns
# for (j in 1:nrow(reruns_df)) {
#   gcm <- reruns_df$gcm[j]
#   ssp <- reruns_df$ssp[j]
#   prd <- reruns_df$prd[j]
#   
#   cmb <- paste0(ssp,'_',gcm,'_',prd)
#   prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
#   yrs <- prd_num[1]:prd_num[2]
#   mns <- c(paste0('0',1:9),10:12)
#   stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
#   names(stp) <- c('yrs','mns')
#   stp <- stp %>%
#     dplyr::arrange(yrs, mns) %>%
#     base::as.data.frame()
#   
#   #folders
#   tx_pth <- paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd) # Maximum temperature
#   rh_pth <- paste0(root,'/chirts/RHum')                                   # Relative humidity
#   out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/THI')
#   
#   yrs_mpg <- data.frame(Baseline = as.character(rep(1995:2014, 2)),
#                         Future = as.character(c(2021:2040,2041:2060)))
#   
#   #remove files
#   this_file <- paste0(out_dir, "/THI_max-", reruns_df$yr[j], "-", reruns_df$mn[j], ".tif")
#   system(paste0("rm -f ", this_file))
#   this_file <- paste0(out_dir, "/THI_mean-", reruns_df$yr[j], "-", reruns_df$mn[j], ".tif")
#   system(paste0("rm -f ", this_file))
#   this_file <- paste0(out_dir, "/daily/THI_daily-", reruns_df$yr[j], "-", reruns_df$mn[j], ".tif")
#   system(paste0("rm -f ", this_file))
#   
#   #redo file
#   cat("redoing file", this_file, "\n")
#   calc_thi(yr = reruns_df$yr[j], mn = reruns_df$mn[j])
# }
