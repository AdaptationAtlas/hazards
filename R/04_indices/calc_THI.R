## Heat stress livestock (cattle) (THI)
## By: H. Achicanoy
## July, 2025

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

# Extent CHIRPS
msk <- terra::rast(file.path(root,'chirps_wrld/chirps-v2.0.1981.01.01.tif'))
xtd <- terra::ext(msk)

sce_climate <- 'future' # historical, future

# Calculate THI function
calc_thi <- function(yr, mn){
  
  outfile1 <- paste0(out_dir,'/THI_mean-',yr,'-',mn,'.tif')
  outfile2 <- paste0(out_dir,'/THI_max-',yr,'-',mn,'.tif')
  outfile3 <- paste0(out_dir,'/daily/THI_daily-',yr,'-',mn,'.tif')
  
  file.exists(c(outfile1,outfile2,outfile3))
  
  cat(outfile2,'\n')
  
  if(sum(file.exists(c(outfile1,outfile2,outfile3))) < 3){
    
    dir.create(dirname(outfile3),F,T)
    
    # Sequence of dates
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01'))) # Last day of the month
    if(as.numeric(yr) > 2020 & mn == '02'){
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-28')), by = 'day')
    } else {
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    }
    
    # Files
    tx_fls <- paste0(tx_pth,'/tasmax_',dts,'.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    rh_fls <- paste0(rh_pth,'/hurs_',dts,'.tif')
    rh_fls <- rh_fls[file.exists(rh_fls)]
    
    # Read variables
    tmx <- terra::rast(tx_fls) |> terra::crop(xtd)
    rhm <- terra::rast(rh_fls) |> terra::crop(xtd)
    
    thr_hum_idx <- function(tmax, rhum){
      thi = (1.8 * tmax + 32) - ((0.55 - 0.0055 * rhum) * (1.8 * tmax - 26.8))
      return(thi)
    }
    cthr_hum_idx <- compiler::cmpfun(thr_hum_idx)
    # Calculate human heat stress
    THI <- cthr_hum_idx(tmax = tmx, rhum = rhm)
    THI_avg <- mean(THI)
    THI_max <- max(THI)
    
    # Write output
    terra::writeRaster(THI_avg, outfile1)
    terra::writeRaster(THI_max, outfile2)
    terra::writeRaster(THI, outfile3)
    
    # Clean up
    rm(tmx, rhm, THI, THI_avg, THI_max)
    gc(verbose=F, full=T, reset=T)
  }
}

# Runs
if (sce_climate == 'future') {
  
  ssps <- c('ssp126','ssp245','ssp370','ssp585') 
  gcms <- c('ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CMCC-ESM2','EC-Earth3','EC-Earth3-Veg-LR','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','TaiESM1')
  
  for (ssp in ssps) {
    for (gcm in gcms) {
      
      ## Parameters
      cmb <- paste0(ssp, '_', gcm)
      yrs <- 2021:2100
      mnt <- sprintf('%02.0f',1:12)
      stp <- base::expand.grid(yrs, mnt, stringsAsFactors = F) |> setNames(c('yrs', 'mnt')) |> dplyr::arrange(yrs, mnt) |> base::as.data.frame(); rm(yrs, mnt)
      
      ## Setup in/out files
      tx_pth <- paste0(root, '/nex-gddp-cmip6/tasmax/', ssp, '/', gcm) # Daily maximum temperatures
      rh_pth <- paste0(root, '/nex-gddp-cmip6/hurs/', ssp, '/', gcm) # Daily relative humidity
      out_dir <- paste0(root,'/nex-gddp-cmip6_indices/',ssp,'_',gcm,'/THI')
      
      1:nrow(stp) |> purrr::map(.f = function(i){calc_thi(yr = stp$yrs[i], mn = stp$mnt[i])})
      tmpfls <- list.files(tempdir(), full.names = T)
      1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ', tmpfls[k]))})
      ##
      cat('----Finish----\n')
      
    }
  }
  
} else {
  if (sce_climate == 'historical') {
    
    ssp  <- 'historical'
    gcms <- c('ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CMCC-ESM2','EC-Earth3','EC-Earth3-Veg-LR','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','TaiESM1')
    
    for (gcm in gcms) {
      
      ## Parameters
      cmb <- paste0(ssp, '_', gcm)
      yrs <- 1995:2014
      mnt <- sprintf('%02.0f',1:12)
      stp <- base::expand.grid(yrs, mnt, stringsAsFactors = F) |> setNames(c('yrs', 'mnt')) |> dplyr::arrange(yrs, mnt) |> base::as.data.frame(); rm(yrs, mnt)
      
      ## Setup in/out files
      tx_pth <- paste0(root, '/nex-gddp-cmip6/tasmax/', ssp, '/', gcm) # Daily maximum temperatures
      rh_pth <- paste0(root, '/nex-gddp-cmip6/hurs/', ssp, '/', gcm) # Daily relative humidity
      out_dir <- paste0(root,'/nex-gddp-cmip6_indices/',ssp,'_',gcm,'/THI')
      
      1:nrow(stp) |> purrr::map(.f = function(i){calc_thi(yr = stp$yrs[i], mn = stp$mnt[i])})
      tmpfls <- list.files(tempdir(), full.names = T)
      1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ', tmpfls[k]))})
      ##
      cat('----Finish----\n')
      
    }
    
  }
  
}
