## Human heat stress (HSH)
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

tm_pth <- paste0(root,'/chirts/Tmin') # Minimum temperature
tx_pth <- paste0(root,'/chirts/Tmax') # Maximum temperature
rh_pth <- paste0(root,'/chirts/RHum') # Relative humidity

# Calculate HSH function
calc_hsh <- function(yr, mn){
  outfile1 <- paste0(root,'/atlas_hazards/cmip6/indices/historic/HSH/HSH_mean-',yr,'-',mn,'.tif')
  outfile2 <- paste0(root,'/atlas_hazards/cmip6/indices/historic/HSH/HSH_max-',yr,'-',mn,'.tif')
  outfile3 <- paste0(root,'/atlas_hazards/cmip6/indices/historic/HSH/daily/HSH_daily-',yr,'-',mn,'.tif')
  file.exists(c(outfile1,outfile2,outfile3))
  if(sum(file.exists(c(outfile1,outfile2,outfile3))) < 3){
    dir.create(dirname(outfile3),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    # Files
    tx_fls <- paste0(tx_pth,'/',yr,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tm_fls <- paste0(tm_pth,'/',yr,'/Tmin.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    rh_fls <- paste0(rh_pth,'/',yr,'/RH.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    rh_fls <- rh_fls[file.exists(rh_fls)]
    # Read variables
    tmx <- terra::rast(tx_fls)
    tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tmx[tmx == -9999] <- NA
    tmn <- terra::rast(tm_fls)
    tmn <- tmn %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tmn[tmn == -9999] <- NA
    tav <- (tmx + tmn)/2
    rhm <- terra::rast(rh_fls)
    rhm <- rhm %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    
    wet_bulb_temp <- function(tmean, rhum){
      twb <- tmean * atan(0.151977 * (rhum + 8.313659)^0.5) + atan(tmean + rhum) - atan(rhum - 1.676331) + 0.00391838 * rhum ^(3/2) * atan(0.023101 * rhum) - 4.686035
      return(twb)
    }
    # Calculate human heat stress
    Twb <- terra::lapp(x = terra::sds(tav,rhm), fun = wet_bulb_temp)
    terra::time(Twb) <- dts
    Twb_avg <- mean(Twb, na.rm = T) %>% terra::mask(ref)
    Twb_max <- max(Twb, na.rm = T) %>% terra::mask(ref)
    
    terra::writeRaster(Twb_avg, outfile1)
    terra::writeRaster(Twb_max, outfile2)
    terra::writeRaster(Twb, outfile3)
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
  purrr::map(.f = function(i){calc_hsh(yr = stp$yrs[i], mn = stp$mns[i]); gc(verbose=F, full=T, reset=T)})
