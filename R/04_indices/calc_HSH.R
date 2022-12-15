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

# Calculate HSH function
calc_hsh <- function(yr, mn){
  outfile1 <- paste0(out_dir,'/HSH_mean-',yr,'-',mn,'.tif')
  outfile2 <- paste0(out_dir,'/HSH_max-',yr,'-',mn,'.tif')
  outfile3 <- paste0(out_dir,'/daily/HSH_daily-',yr,'-',mn,'.tif')
  file.exists(c(outfile1,outfile2,outfile3))
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
    tm_fls <- paste0(tm_pth,'/',yr,'/Tmin.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    if(as.numeric(yr) > 2020){
      dts_chr <- as.character(dts)
      dts_chr <- gsub(pattern = yr, replacement = yrs_mpg$Baseline[yrs_mpg$Future == yr], x = dts_chr)
      dts_h <- as.Date(dts_chr); rm(dts_chr)
      rh_fls <- paste0(rh_pth,'/',yr,'/RH.',gsub(pattern='-', replacement='.', x=dts_h, fixed=T),'.tif')
    } else {
      rh_fls <- paste0(rh_pth,'/',yr,'/RH.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    }
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
    # Constants
    c1 = -8.78469475556
    c2 =  1.61139411
    c3 =  2.33854883889
    c4 = -0.14611605
    c5 = -0.012308094
    c6 = -0.0164248277778
    c7 =  2.211732 * 10^(-3)
    c8 =  7.2546 * 10^(-4)
    c9 = -3.582 * 10^(-6)
    heat_idx <- function(tmean, rhum){
      hi <- ifelse(tmean >= 25,
                   c1 + (c2*tmean) + (c3*rhum) + (c4*tmean*rhum) + (c5*tmean^2) + (c6*rhum^2) + (c7*tmean^2*rhum) + (c8*tmean*rhum^2) + (c9*tmean^2*rhum^2),
                   tmean)
      return(hi)
    }
    # Calculate human heat stress
    HI <- terra::lapp(x = terra::sds(tav,rhm), fun = heat_idx)
    terra::time(HI) <- dts
    HI_avg <- mean(HI, na.rm = T) %>% terra::mask(ref)
    HI_max <- max(HI, na.rm = T) %>% terra::mask(ref)
    
    terra::writeRaster(HI_avg, outfile1)
    terra::writeRaster(HI_max, outfile2)
    terra::writeRaster(HI, outfile3)
  }
}

# Historical setup
yrs <- 1995:2014
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()
tm_pth <- paste0(root,'/chirts/Tmin') # Minimum temperature
tx_pth <- paste0(root,'/chirts/Tmax') # Maximum temperature
rh_pth <- paste0(root,'/chirts/RHum') # Relative humidity
out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/historical/HSH')

# # Future setup
# gcm <- 'ACCESS-ESM1-5'
# ssp <- 'ssp245'
# prd <- '2021_2040'
# 
# cmb <- paste0(ssp,'_',gcm,'_',prd)
# prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
# yrs <- prd_num[1]:prd_num[2]
# mns <- c(paste0('0',1:9),10:12)
# stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
# names(stp) <- c('yrs','mns')
# stp <- stp %>%
#   dplyr::arrange(yrs, mns) %>%
#   base::as.data.frame()
# tm_pth <- paste0(root,'/chirts_cmip6_africa/Tmin_',gcm,'_',ssp,'_',prd) # Minimum temperature
# tx_pth <- paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd) # Maximum temperature
# rh_pth <- paste0(root,'/chirts/RHum')                                   # Relative humidity
# out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/HSH')
# 
# yrs_mpg <- data.frame(Baseline = as.character(rep(1995:2014, 2)),
#                       Future = as.character(c(2021:2040,2041:2060)))

1:nrow(stp) %>%
  purrr::map(.f = function(i){calc_hsh(yr = stp$yrs[i], mn = stp$mns[i]); gc(verbose=F, full=T, reset=T)})
