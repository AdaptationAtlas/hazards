## Thornthwaite's Aridity Index (TAI)
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate,envirem))
source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

# Calculate TAI function
calc_tai <- function(yr){
  outfile <- paste0(out_dir,'/TAI-',yr,'.tif')
  if(!file.exists(outfile)){
    dir.create(dirname(outfile),F,T)
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-01-01')), to = as.Date(paste0(yr,'-12-31')), by = 'day')
    # Files
    pr_fls <- paste0(pr_pth,'/chirps-v2.0.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    pr_fls <- pr_fls[file.exists(pr_fls)]
    tx_fls <- paste0(tx_pth,'/',yr,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tm_fls <- paste0(tm_pth,'/',yr,'/Tmin.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tm_fls <- tm_fls[file.exists(tm_fls)]
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
    # Calculate average temperature
    tav <- (tmx + tmn)/2
    # Calculate temperature range
    rnge <- abs(tmx - tmn)
    # Calculate monthly summaries
    prc_month <- terra::tapp(x = prc, index = lubridate::month(dts), fun = sum)
    tav_month <- terra::tapp(x = tav, index = lubridate::month(dts), fun = mean)
    rng_month <- terra::tapp(x = rnge, index = lubridate::month(dts), fun = mean)
    # ET SRAD
    srf <- list.dirs(paste0(root,'/ET_SolRad'), full.names = T, recursive = F)
    srf <- srf[-length(srf)]
    srf <- srf %>% gtools::mixedsort()
    srd <- srf %>% raster::stack()
    srd <- srd %>% raster::resample(raster::raster(ref))
    srd <- srd %>% raster::crop(raster(ref))
    srd <- srd %>% raster::mask(raster(ref))
    names(srd) <- c(paste0('SRAD_0',1:9),paste0('SRAD_', 10:12))
    
    # Assign precipitation names in envirem environment
    envirem::assignNames(solrad='SRAD_##', tmean = 'TMEAN_##', precip = 'PREC_##')
    
    TMEAN <- tav_month %>% raster::stack()
    PREC  <- prc_month %>% raster::stack()
    TRNG  <- rng_month %>% raster::stack()
    
    names(TMEAN) <- c(paste0('TMEAN_0',1:9),paste0('TMEAN_', 10:12))
    names(PREC)  <- c(paste0('PREC_0',1:9),paste0('PREC_', 10:12))
    names(TRNG)  <- c(paste0('TRNG_0',1:9),paste0('TRNG_', 10:12))
    
    # Thornthwaite's Aridity Index
    PET <- envirem::monthlyPET(TMEAN, srd, TRNG) %>% raster::stack()
    names(PET)  <- c(paste0('PET_0',1:9), paste0('PET_', 10:12))
    PET <- raster::resample(PET, PREC[[1]])
    TAI <- envirem::aridityIndexThornthwaite(PREC, PET)
    TAI <- terra::rast(TAI)
    names(TAI) <- yr
    
    terra::writeRaster(TAI, outfile)
  }
}

# Historical setup
stp <- data.frame(yrs = 1995:2014)
pr_pth <- paste0(root,'/chirps_wrld') # Precipitation
tm_pth <- paste0(root,'/chirts/Tmin') # Minimum temperature
tx_pth <- paste0(root,'/chirts/Tmax') # Maximum temperature
# out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/historical/TAI')
out_dir <- '/home/jovyan/indices/historical/TAI'

# # Future setup
# gcm <- 'ACCESS-ESM1-5'
# ssp <- 'ssp245'
# prd <- '2021_2040'
# 
# cmb <- paste0(ssp,'_',gcm,'_',prd)
# prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
# stp <- data.frame(yrs = prd_num[1]:prd_num[2])
# pr_pth <- paste0(root,'/chirps_cmip6_africa/Prec_',gcm,'_',ssp,'_',prd) # Precipitation
# tm_pth <- paste0(root,'/chirts_cmip6_africa/Tmin_',gcm,'_',ssp,'_',prd) # Minimum temperature
# tx_pth <- paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd) # Maximum temperature
# out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/TAI')

1:nrow(stp) %>%
  purrr::map(.f = function(i){calc_tai(yr = stp$yrs[i]); gc(verbose=F, full=T, reset=T)})
