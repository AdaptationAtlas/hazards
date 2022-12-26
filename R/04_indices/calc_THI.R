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

# Calculate THI function
calc_thi <- function(yr, mn){
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
    if(as.numeric(yr) > 2020 & mn == '02'){
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-28')), by = 'day')
    } else {
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    }
    # Files
    tx_fls <- paste0(tx_pth,'/',yr,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    if(as.numeric(yr) > 2020){
      dts_chr <- as.character(dts)
      dts_chr <- gsub(pattern = yr, replacement = yrs_mpg$Baseline[yrs_mpg$Future == yr], x = dts_chr)
      dts_h <- as.Date(dts_chr); rm(dts_chr)
      yr_h <- lubridate::year(dts_h)[1]
      rh_fls <- paste0(rh_pth,'/',yr_h,'/RH.',gsub(pattern='-', replacement='.', x=dts_h, fixed=T),'.tif')
    } else {
      rh_fls <- paste0(rh_pth,'/',yr,'/RH.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    }
    rh_fls <- rh_fls[file.exists(rh_fls)]
    # Read variables
    tmx <- terra::rast(tx_fls)
    tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tmx[tmx == -9999] <- NA
    rhm <- terra::rast(rh_fls)
    rhm <- rhm %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    
    thr_hum_idx <- function(tmax, rhum){
      thi = (1.8 * tmax + 32) - ((0.55 - 0.0055 * rhum) * (1.8 * tmax - 26.8))
      return(thi)
    }
    # Calculate human heat stress
    THI <- terra::lapp(x = terra::sds(tmx,rhm), fun = thr_hum_idx)
    terra::time(THI) <- dts
    THI_avg <- mean(THI, na.rm = T) %>% terra::mask(ref)
    THI_max <- max(THI, na.rm = T) %>% terra::mask(ref)
    
    terra::writeRaster(THI_avg, outfile1)
    terra::writeRaster(THI_max, outfile2)
    terra::writeRaster(THI, outfile3)
  }
}

# # Historical setup
# yrs <- 1995:2014
# mns <- c(paste0('0',1:9),10:12)
# stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
# names(stp) <- c('yrs','mns')
# stp <- stp %>%
#   dplyr::arrange(yrs, mns) %>%
#   base::as.data.frame()
# tx_pth <- paste0(root,'/chirts/Tmax') # Maximum temperature
# rh_pth <- paste0(root,'/chirts/RHum') # Relative humidity
# out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/historical/THI')
# 1:nrow(stp) %>%
#   purrr::map(.f = function(i){calc_thi(yr = stp$yrs[i], mn = stp$mns[i]); gc(verbose=F, full=T, reset=T)})

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
