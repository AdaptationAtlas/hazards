## Number of soil water stress days (NDWS)
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))
source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

# Soil variables
scp <- terra::rast(paste0(root,'/atlas_hazards/soils/africa_scp.tif'))
scp <- scp %>% terra::resample(ref) %>% terra::mask(ref)
sst <- terra::rast(paste0(root,'/atlas_hazards/soils/africa_ssat.tif'))
sst <- sst %>% terra::resample(ref) %>% terra::mask(ref)

pr_pth <- paste0(root,'/chirps_wrld') # Precipitation
tm_pth <- paste0(root,'/chirts/Tmin') # Minimum temperature
tx_pth <- paste0(root,'/chirts/Tmax') # Maximum temperature
sr_pth <- paste0(root,'/ecmwf_agera5/solar_radiation_flux') # Solar radiation

# Calculate NDWS function
calc_ndws <- function(yr, mn){
  outfile <- paste0(root,'/atlas_hazards/cmip6/indices/historic/NDWS/NDWS-',yr,'-',mn,'.tif')
  if(!file.exists(outfile)){
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    # Files
    pr_fls <- paste0(pr_pth,'/chirps-v2.0.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    pr_fls <- pr_fls[file.exists(pr_fls)]
    tx_fls <- paste0(tx_pth,'/',yr,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tm_fls <- paste0(tm_pth,'/',yr,'/Tmin.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    sr_fls <- paste0(sr_pth,'/Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0.nc')
    sr_fls <- sr_fls[file.exists(sr_fls)]
    # Read variables
    prc <- terra::rast(pr_fls)
    prc <- prc %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    prc[prc == -9999] <- NA
    tmx <- terra::rast(tx_fls)
    tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tmx[tmx == -9999] <- NA
    tmn <- terra::rast(tm_fls)
    tmn <- tmn %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tmn[tmn == -9999] <- NA
    tav <- (tmx + tmn)/2
    srd <- terra::rast(sr_fls)
    srd <- srd %>% terra::crop(terra::ext(ref))
    srd <- srd/1000000
    srd <- srd %>% terra::resample(x = ., y = ref) %>% terra::mask(ref)
    
    # Maximum evapotranspiration
    ETMAX <- terra::lapp(x = terra::sds(srd,tmn,tav,tmx), fun = peest)
    
    # Compute water balance model
    date <- paste0(yr,'-',mn)
    if(date == '1995-01'){
      AVAIL <<- ref
      AVAIL[!is.na(AVAIL)] <- 0
    } else {
      AVAIL <<- terra::rast(paste0(dirname(outfile),'/AVAIL.tif'))
    }
    
    watbal <- 1:terra::nlyr(ETMAX) %>%
      purrr::map(.f = function(i){
        water_balance <- eabyep_calc(soilcp  = scp,
                                     soilsat = sst,
                                     avail   = AVAIL[[terra::nlyr(AVAIL)]],
                                     rain    = prc[[i]],
                                     evap    = ETMAX[[i]])
        AVAIL <<- water_balance$Availability
        return(water_balance)
      })
    ERATIO <- watbal %>% purrr::map('Eratio') %>% terra::rast()
    # Calculate number of soil water stress days
    NDWS   <- terra::app(x = ERATIO, fun = function(ERATIO){ifelse(ERATIO < 0.5, 1, 0)}) %>% sum()
    terra::writeRaster(NDWS, outfile)
    terra::writeRaster(AVAIL[[terra::nlyr(AVAIL)]], paste0(dirname(outfile),'/AVAIL.tif'), overwrite = T)
  }
}

# Setup
yrs <- 1995:2014
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()

1:nrow(stp) %>%
  purrr::map(.f = function(i){calc_ndws(yr = stp$yrs[i], mn = stp$mns[i])})
