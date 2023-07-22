## Average monthly maximum temperature
## By: H. Achicanoy / Julian Ramirez-Villegas
## July, 2023

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

sce_climate <- "future" #"future"

# Calculate tmax function
calc_tmx <- function(yr, mn, dataset="CHIRPS"){
  outfile <- paste0(out_dir,'/TMAX-',yr,'-',mn,'.tif')
  if(!file.exists(outfile)){
    cat(outfile, "\n")
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    
    # Files
    if (dataset == "CHIRPS") {
      tx_fls <- paste0(tx_pth,'/',yr,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    } else if (dataset == "AgERA5") {
      tx_fls <- paste0(ae5tx_pth,'/Temperature-Air-2m-Max-24h_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0.nc')
    }
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tmx <- terra::rast(tx_fls)
    if (dataset == "AgERA5") {
      tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::resample(., ref) %>% terra::mask(ref)
      tmx <- tmx - 273.15
    } else {
      tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
      tmx[tmx == -9999] <- NA
    }
    
    # Calculate heat stress generic crop
    terra::app(x   = tmx,
               fun = function(x){ tmax = mean(x, na.rm = T); return(tmax) },
               filename = outfile)
    # Clean up
    rm(tmx)
    gc(verbose=FALSE, full=TRUE, reset=TRUE)
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
  tx_pth <- paste0(root,'/chirts/Tmax') # Daily maximum temperature
  ae5tx_pth <- paste0(root,'/ecmwf_agera5/2m_temperature-24_hour_maximum') # Maximum temperature
  out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/historical/TMAX_AgERA5')
  1:nrow(stp) %>%
    purrr::map(.f = function(i){
      calc_tmx(yr = stp$yrs[i], mn = stp$mns[i], dataset="AgERA5")
      tmpfls <- list.files(tempdir(), full.names=TRUE)
      1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
    })
} else if (sce_climate == "future") {
  # Future setup
  gcm <- 'MRI-ESM2-0' #'ACCESS-ESM1-5' 'MPI-ESM1-2-HR' 'EC-Earth3' 'INM-CM5-0' 'MRI-ESM2-0'
  for (ssp in c('ssp245', 'ssp585')) { #'ssp126' 'ssp370'
    for (prd in c('2021_2040', '2041_2060')) { #'2061_2080', '2081_2100'
      cat("...processing gcm=", gcm, "/ ssp=", ssp, "/ period=", prd, "\n")
      #ssp <- 'ssp245' #ssp585
      #prd <- '2041_2060' #2021_2040
      
      cmb <- paste0(ssp,'_',gcm,'_',prd)
      prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
      yrs <- prd_num[1]:prd_num[2]
      mns <- c(paste0('0',1:9),10:12)
      stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
      names(stp) <- c('yrs','mns')
      stp <- stp %>%
        dplyr::arrange(yrs, mns) %>%
        base::as.data.frame()
      tx_pth <- paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd) # Daily maximum temperatures
      out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/TMAX')
      
      1:nrow(stp) %>%
        purrr::map(.f = function(i){
          calc_tmx(yr = stp$yrs[i], mn = stp$mns[i])
          tmpfls <- list.files(tempdir(), full.names=TRUE)
          1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
        })
    }
  }
} else {
  cat("select one of historical or future for sce_climate \n")
}

