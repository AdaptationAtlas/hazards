#!/usr/bin Rscript
args = commandArgs(trailingOnly=TRUE)

gcm <- args[1]
ssp <- args[2]
prd <- args[3]
yr_i <- as.numeric(args[4])
yr_f <- as.numeric(args[5])

#verbose what is being processed
cat("gcm=", gcm, "/ ssp=", ssp, "/ period=", prd, "/ yr_i=", yr_i, "/ yr_f=", yr_f, "\n")

## Thornthwaite's Aridity Index (TAI)
## By: H. Achicanoy, modified by JRV so that it could be run via Rscript
## December, 2022

#how to run it: just paste the below in the console (modifying gcm, ssp, period, initial and end year)
#Rscript --vanilla ~/Repositories/hazards/R/04_indices/calc_TAI_sh.R ACCESS-ESM1-5 ssp245 2021_2040 1 10
#you can also write a bash (.sh) script with the list of jobs to run and execute it.

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate,envirem))

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

# ET SRAD, load and process only once
srf <- list.dirs(paste0(root,'/ET_SolRad'), full.names = T, recursive = F)
srf <- srf[-length(srf)]
srf <- srf %>% gtools::mixedsort()
srd <- srf %>% raster::stack()
srd <- srd %>% raster::resample(raster::raster(ref))
srd <- srd %>% raster::crop(raster(ref))
srd <- srd %>% raster::mask(raster(ref))
names(srd) <- c(paste0('SRAD_0',1:9),paste0('SRAD_', 10:12))

# Calculate TAI function
calc_tai <- function(yr){
  outfile <- paste0(out_dir,'/TAI-',yr,'.tif')
  cat(outfile, "\n")
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
    
    # Read variables, and calculate monthly summaries, do by month to reduce memory consumption
    prc_ls <- tav_ls <- rng_ls <- c()
    for (j in 1:12) {
        mth_fls <- paste0(pr_pth, '/chirps-v2.0.', yr, '.', sprintf("%02.0f", j), ".", sprintf("%02.0f", 1:31), '.tif')
        prc <- terra::rast(pr_fls[pr_fls %in% mth_fls])
        prc <- prc %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
        prc[prc < 0] <- NA
        prc_month <- sum(prc, na.rm=TRUE)
        rm(prc)
        gc(verbose=F, full=T, reset=T)
        
        mth_fls <- paste0(tx_pth, '/', yr, '/Tmax.', yr, '.', sprintf("%02.0f", j), ".", sprintf("%02.0f", 1:31), '.tif')
        tmx <- terra::rast(tx_fls[tx_fls %in% mth_fls])
        tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
        tmx[tmx == -9999] <- NA
        
        mth_fls <- paste0(tm_pth, '/', yr, '/Tmin.', yr, '.', sprintf("%02.0f", j), ".", sprintf("%02.0f", 1:31), '.tif')
        tmn <- terra::rast(tm_fls[tm_fls %in% mth_fls])
        tmn <- tmn %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
        tmn[tmn == -9999] <- NA

        # Calculate average temperature
        tav <- (tmx + tmn)/2
        tav_month <- mean(tav, na.rm=TRUE)
        rm(tav)
        gc(verbose=F, full=T, reset=T)

        # Calculate temperature range
        rnge <- abs(tmx - tmn)
        rng_month <- mean(rnge, na.rm=TRUE)
        rm(rnge, tmx, tmn)
        gc(verbose=F, full=T, reset=T)
        
        # Append
        prc_ls <- c(prc_ls, prc_month)
        tav_ls <- c(tav_ls, tav_month)
        rng_ls <- c(rng_ls, rng_month)
    }
    
    # Assign precipitation names in envirem environment
    envirem::assignNames(solrad='SRAD_##', tmean = 'TMEAN_##', precip = 'PREC_##')
    
    TMEAN <- tav_ls %>% terra::rast() %>% raster::stack()
    PREC  <- prc_ls %>% terra::rast() %>% raster::stack()
    TRNG  <- rng_ls %>% terra::rast() %>% raster::stack()
    
    names(TMEAN) <- c(paste0('TMEAN_0',1:9),paste0('TMEAN_', 10:12))
    names(PREC)  <- c(paste0('PREC_0',1:9),paste0('PREC_', 10:12))
    names(TRNG)  <- c(paste0('TRNG_0',1:9),paste0('TRNG_', 10:12))
    
    #clean up
    rm(tav_ls, prc_ls, rng_ls, tav_month, prc_month, rng_month)
    gc(verbose=F, full=T, reset=T)
    
    # Thornthwaite's Aridity Index
    PET <- envirem::monthlyPET(TMEAN, srd, TRNG) %>% raster::stack()
    names(PET)  <- c(paste0('PET_0',1:9), paste0('PET_', 10:12))
    PET <- raster::resample(PET, PREC[[1]])
    TAI <- envirem::aridityIndexThornthwaite(PREC, PET)
    TAI <- terra::rast(TAI)
    names(TAI) <- yr
    
    #write output
    terra::writeRaster(TAI, outfile)
    
    #clean up
    rm(PET, TMEAN, PREC, TRNG, TAI)
    gc(verbose=F, full=T, reset=T)
  }
}

# Options
#'ACCESS-ESM1-5' 'MPI-ESM1-2-HR' 'EC-Earth3' 'INM-CM5-0' 'MRI-ESM2-0'
#'ssp245' 'ssp585'
#'2021_2040' '2041_2060'

# Run the function
cmb <- paste0(ssp,'_',gcm,'_',prd)
prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
stp <- data.frame(yrs = prd_num[1]:prd_num[2])
pr_pth <- paste0(root,'/chirps_cmip6_africa/Prec_',gcm,'_',ssp,'_',prd) # Precipitation
tm_pth <- paste0(root,'/chirts_cmip6_africa/Tmin_',gcm,'_',ssp,'_',prd) # Minimum temperature
tx_pth <- paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd) # Maximum temperature
out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/TAI')
yr_i:yr_f %>%
  purrr::map(.f = function(i){calc_tai(yr = stp$yrs[i])})
