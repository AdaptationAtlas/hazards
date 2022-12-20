## Get daily future data
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate,furrr))

root <- '/home/jovyan/common_data'

# Africa reference raster
ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

# Interpolated monthly anomalies directory
anm_pth <- paste0(root,'/esfg_cmip6/intermediate/interpolated_mthly_anomaly')

# Setup
gcms <- c('ACCESS-ESM1-5','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0') # 'EC-Earth3'
ssps <- c('ssp245','ssp585')
vrss <- c('pr','tasmax','tasmin')
prds <- c('2021_2040','2041_2060')
stp <- base::expand.grid(gcms,ssps,vrss,prds) %>% base::as.data.frame()
names(stp) <- c('gcm','ssp','var','prd'); rm(gcms, ssps, vrss, prds)
stp <- stp %>%
  dplyr::arrange(gcm,ssp,prd,var) %>%
  base::as.data.frame()

# gcm <- gcms[1]
# ssp <- ssps[1]
# var <- vrss[1]
# prd <- prds[1]

# Read monthly deltas
get_daily_future_data <- function(gcm, ssp, var, prd){
  prd <- as.character(prd)
  file <- paste0('CMIP6_',gcm,'_',ssp,'_r1i1p1f1_',var,'_Africa_monthly_intp_anomaly_',prd,'.tif')
  # Load deltas
  dlts <- terra::rast(paste0(anm_pth,'/',file))
  dlts <- dlts %>% terra::resample(ref)
  dlts <- dlts %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
  prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
  his_yrs <- 1995:2014
  fut_yrs <- prd_num[1]:prd_num[2]
  
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
  
  if(var == 'pr'){
    # Paths
    his_pth <- paste0(root,'/chirps_wrld')
    fut_pth <- paste0(root,'/chirps_cmip6_africa/Prec_',gcm,'_',ssp,'_',prd); dir.create(fut_pth, F, T)
    if(length(list.files(fut_pth)) < 7300){
      # File structure
      his_str <- paste0('chirps-v2.0.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
      fut_str <- paste0('chirps-v2.0.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      # Split by months
      his_lst <- split(his_str, mpg$month)
      fut_lst <- split(fut_str, mpg$month)
      1:length(his_lst) %>%
        purrr::map(.f = function(j){
          delta <- dlts[[j]]
          his_daily <- his_lst[[j]]
          fut_daily <- fut_lst[[j]]
          plan(multicore, workers = 20)
          1:length(his_daily) %>%
            furrr::future_map(.f = function(k){
              outfile <- paste0(fut_pth,'/',fut_daily[k])
              if(!file.exists(outfile)){
                r <- terra::rast(paste0(his_pth,'/',his_daily[k]))
                r <- r %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
                r[r == -9999] <- NA
                r <- r * (1 + delta)
                terra::writeRaster(r, outfile)
              }
            })
          future:::ClusterRegistry("stop"); gc(reset = T)
        })
    }
  }
  if(var %in% c('tasmax','tasmin')){
    # Paths
    his_pth <- ifelse(var == 'tasmax',
                      paste0(root,'/chirts/Tmax'),
                      paste0(root,'/chirts/Tmin'))
    fut_pth <- ifelse(var == 'tasmax',
                      paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd),
                      paste0(root,'/chirts_cmip6_africa/Tmin_',gcm,'_',ssp,'_',prd))
    if(length(list.files(fut_pth)) < 7300){
      # File structure
      if(var == 'tasmax'){
        his_str <- paste0('Tmax.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('Tmax.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      } else {
        his_str <- paste0('Tmin.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('Tmin.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      }
      yrs_str   <- mpg$year
      yrs_f_str <- mpg$year_fut
      # Split by months
      his_lst   <- split(his_str, mpg$month)
      fut_lst   <- split(fut_str, mpg$month)
      yrs_lst   <- split(yrs_str, mpg$month)
      yrs_f_lst <- split(yrs_f_str, mpg$month)
      1:length(his_lst) %>%
        purrr::map(.f = function(j){
          delta <- dlts[[j]]
          his_daily   <- his_lst[[j]]
          fut_daily   <- fut_lst[[j]]
          yrs_daily   <- yrs_lst[[j]]
          yrs_f_daily <- yrs_f_lst[[j]]
          plan(multicore, workers = 20)
          1:length(his_daily) %>%
            furrr::future_map(.f = function(k){
              outfile <- paste0(fut_pth,'/',yrs_f_daily[k],'/',fut_daily[k]); dir.create(dirname(outfile),F,T)
              if(!file.exists(outfile)){
                r <- terra::rast(paste0(his_pth,'/',yrs_daily[k],'/',his_daily[k]))
                r <- r %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
                r[r == -9999] <- NA
                r <- r + delta
                terra::writeRaster(r, outfile)
              }
            })
          future:::ClusterRegistry("stop"); gc(reset = T)
        })
    }
  }
  
  return(cat(paste0(var,'_',gcm,'_',ssp,'_',prd,' ready.\n')))
  
}

1:nrow(stp) %>%
  purrr::map(.f = function(i){get_daily_future_data(gcm = stp$gcm[i],
                                                    ssp = stp$ssp[i],
                                                    var = stp$var[i],
                                                    prd = stp$prd[i])
    gc(verbose=F, full=T, reset=T)
  })
