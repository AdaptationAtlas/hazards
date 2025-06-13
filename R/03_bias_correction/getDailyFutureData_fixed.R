## Get daily future data
## By: H. Achicanoy
## December, 2022

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate,furrr))

root <- '/home/jovyan/common_data'

# Africa reference raster
ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

# Interpolated monthly anomalies directory
anm_pth <- paste0(root,'/esfg_cmip6/intermediate/interpolated_mthly_anomaly')

# Read monthly deltas
get_daily_future_data <- function(gcm, ssp, var, prd, mn){
  cat(paste0("processing ",var,'_',gcm,'_',ssp,'_',prd,'_',mn,' \n'))
  prd <- as.character(prd)
  file <- paste0('CMIP6_',gcm,'_',ssp,'_r1i1p1f1_',var,'_Africa_monthly_intp_anomaly_',prd,'_fixed.tif')
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
  mpg <- mpg |> dplyr::filter(month == as.numeric(mn))
  
  if(var %in% c('pr','rsds')){
    # Paths
    if (var == "pr") {
      his_pth <- paste0(root,'/chirps_wrld')
      fut_pth <- paste0(root,'/chirps_cmip6_africa/Prec_',gcm,'_',ssp,'_',prd); dir.create(fut_pth, F, T)
    } else if (var == "rsds") {
      his_pth <- paste0(root,'/ecmwf_agera5/solar_radiation_flux')
      fut_pth <- paste0(root,'/ecmwf_agera5_cmip6_africa/solar_radiation_flux_',gcm,'_',ssp,'_',prd)
      dir.create(fut_pth, F, T)
    }
    #if(length(list.files(fut_pth)) < 7300){
      # File structure
      if (var == "pr") {
        his_str <- paste0('chirps-v2.0.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('chirps-v2.0.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      } else if (var == "rsds") {
        his_str <- paste0('Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=mpg$Baseline, fixed=T),'_final-v1.0.nc')
        fut_str <- paste0('Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=mpg$Future, fixed=T),'_final-v1.0.nc')
      }
      # Split by months
      his_lst <- split(his_str, mpg$month)
      fut_lst <- split(fut_str, mpg$month)
      1:length(his_lst) %>%
        purrr::map(.f = function(j){
          delta <- dlts[[j]]
          his_daily <- his_lst[[j]]
          fut_daily <- fut_lst[[j]]
          plan(multicore, workers = 5)
          1:length(his_daily) %>%
            furrr::future_map(.f = function(k){
              outfile <- paste0(fut_pth,'/',fut_daily[k])
              r <- terra::rast(paste0(his_pth,'/',his_daily[k]))
              if (var == "pr") {
                r <- r %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
                r <- terra::classify(r, rcl = cbind(-Inf, -9990, NA))
              } else if (var == "rsds") {
                r <- r %>% terra::crop(terra::ext(ref)) %>% terra::resample(., ref) %>% terra::mask(ref)
                r <- r * 1e-6
              }
              r <- r * (1 + delta)
              r <- terra::classify(r, cbind(-Inf, 0, 0))
              terra::writeRaster(r, outfile, overwrite = T)
            })
          plan(sequential); gc(reset = T)
        })
    #}
  }
  if(var %in% c('tasmax','tasmin','hurs')){
    # Paths
    his_pth <- ifelse(var == 'tasmax',
                      paste0(root,'/chirts/Tmax'),
                      ifelse(var == 'tasmin',
                             paste0(root,'/chirts/Tmin'), 
                             paste0(root,'/chirts/RHum')))
    fut_pth <- ifelse(var == 'tasmax',
                      paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd),
                      ifelse(var == 'tasmin',
                             paste0(root,'/chirts_cmip6_africa/Tmin_',gcm,'_',ssp,'_',prd),
                             paste0(root,'/chirts_cmip6_africa/RHum_',gcm,'_',ssp,'_',prd)))
    if(length(list.files(fut_pth)) < 7300){
      # File structure
      if(var == 'tasmax'){
        his_str <- paste0('Tmax.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('Tmax.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      } else if (var == 'tasmin'){
        his_str <- paste0('Tmin.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('Tmin.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
      } else {
        his_str <- paste0('RH.',gsub(pattern='-', replacement='.', x=mpg$Baseline, fixed=T),'.tif')
        fut_str <- paste0('RH.',gsub(pattern='-', replacement='.', x=mpg$Future, fixed=T),'.tif')
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
          plan(multicore, workers = 5)
          1:length(his_daily) %>%
            furrr::future_map(.f = function(k){
              outfile <- paste0(fut_pth,'/',yrs_f_daily[k],'/',fut_daily[k]); dir.create(dirname(outfile),F,T)
              if(!file.exists(outfile)){
                r <- terra::rast(paste0(his_pth,'/',yrs_daily[k],'/',his_daily[k]))
                r <- r %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
                r[r == -9999] <- NA
                r <- r + delta
                #for hurs limit min to 0, max to 100
                if (var == "hurs") {
                  r <- min(r, 100)
                  r <- max(r, 0)
                }
                terra::writeRaster(r, outfile)
              }
            })
          plan(sequential); gc(reset = T)
        })
    }
  }
  
  return(cat(paste0(var,'_',gcm,'_',ssp,'_',prd,' ready.\n')))
  
}

get_daily_future_data(gcm = gcm, ssp = ssp, var = var, prd = prd, mn = mn)
gc(verbose = F, full = T, reset = T)