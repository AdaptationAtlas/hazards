## Reverse engineering: fix daily data files
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

# Calculate NDD function
fix_day <- function(date = '2041-10-19', var = 'pr', gcm = 'ACCESS-ESM1-5', ssp = 'ssp585', prd = '2041_2060'){
  
  prd  <- as.character(prd)
  date <- as.Date(date)
  yr   <- lubridate::year(date)
  mn   <- lubridate::month(date)
  
  if(var == 'pr'){prfx <- 'Prec'}
  if(var == 'tasmax'){prfx <- 'Tmax'}
  if(var == 'tasmin'){prfx <- 'Tmin'}
  
  if(var == 'pr'){
    fl2fix <- paste0(root,'/chirps_cmip6_africa/',prfx,'_',gcm,'_',ssp,'_',prd,'/chirps-v2.0.',gsub('-','.',as.character(date)),'.tif')
  } else {
    if(var %in% c('tasmax','tasmin')){
      fl2fix <- paste0(root,'chirts_cmip6_africa/',prfx,'_',gcm,'_',ssp,'_',prd,'/',yr,'/',prfx,'.',gsub('-','.',as.character(date)),'.tif')
    }
  }
  
  # Interpolated monthly anomalies directory
  anm_pth <- paste0(root,'/esfg_cmip6/intermediate/interpolated_mthly_anomaly')
  file <- paste0('CMIP6_',gcm,'_',ssp,'_r1i1p1f1_',var,'_Africa_monthly_intp_anomaly_',prd,'.tif')
  
  # Load deltas
  dlts <- terra::rast(paste0(anm_pth,'/',file))
  dlts <- dlts %>% terra::resample(ref)
  dlts <- dlts %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
  his_yrs <- 1995:2014
  prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
  
  # Temporal mapping
  Baseline = seq(from = as.Date('1995-01-01'),
                 to   = as.Date('2014-12-31'),
                 by   = 'day')
  Future   = seq(from = as.Date(paste0(prd_num[1],'-01-01')),
                 to   = as.Date(paste0(prd_num[2],'-12-31')),
                 by   = 'day')
  # Remove feb-29 from all leap years (not coincidence between years)
  Baseline <- Baseline[!(format(Baseline,"%m") == "02" & format(Baseline, "%d") == "29"), drop = FALSE]
  Future   <- Future[!(format(Future,"%m") == "02" & format(Future, "%d") == "29"), drop = FALSE]
  mpg <- data.frame(Baseline, Future)
  mpg$year  <- lubridate::year(mpg$Baseline)
  mpg$year_fut <- lubridate::year(mpg$Future)
  mpg$month <- lubridate::month(mpg$Baseline)
  
  mpg <- mpg[mpg$Future == date,]
  
  delta <- dlts[[mn]]
  
  if(var == 'pr'){
    
    his_pth  <- paste0(root,'/chirps_wrld')
    his_file <- paste0('chirps-v2.0.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
    
    r <- terra::rast(paste0(his_pth,'/',his_file))
    r <- r %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    r[r == -9999] <- NA
    r <- r * (1 + delta)
    terra::writeRaster(r, fl2fix, overwrite = T)
    
    return(cat(paste0(basename(fl2fix),' fixed\n')))
    
  }
  if(var %in% c('tasmax','tasmin')){
    
    his_pth <- ifelse(var == 'tasmax',
                      paste0(root,'/chirts/Tmax/',mpg$year),
                      paste0(root,'/chirts/Tmin/',mpg$year))
    
    if(var == 'tasmax'){
      his_file <- paste0('Tmax.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
    } else {
      his_file <- paste0('Tmin.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
    }
    
    r <- terra::rast(paste0(his_pth,'/',his_file))
    r <- r %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    r[r == -9999] <- NA
    r <- r + delta
    terra::writeRaster(r, fl2fix, overwrite = T)
    
    return(cat(paste0(basename(fl2fix),' fixed\n')))
  }
  
}

fix_day(date = '2041-10-19', var = 'pr', gcm = 'ACCESS-ESM1-5', ssp = 'ssp585', prd = '2041_2060')
