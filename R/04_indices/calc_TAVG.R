## Average monthly temperature
## By: H. Achicanoy / J. Ramirez-Villegas
## July, 2025

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

# Arguments
# args <- commandArgs(trailingOnly = T)

# Root folder
root <- '/home/jovyan/common_data'

# Extent CHIRPS
msk <- terra::rast('/home/jovyan/common_data/chirps_wrld/chirps-v2.0.1981.01.01.tif')
xtd <- terra::ext(msk)

sce_climate <- 'future' # historical, future

# Function
calc_tav <- function(yr, mn){
  
  # yr <- '2021'; mn <- '01'
  
  outfile <- paste0(out_dir,'TAVG-',yr,'-',mn,'.tif')
  cat(basename(outfile),'\n')
  
  if(!file.exists(outfile)){
    
    ## Create output directory
    dir.create(dirname(outfile), F, T)
    
    ## Tidy the dates
    cat('Processing -------> ', yr, mn, '\n')
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    
    txfls <- paste0(tx_pth,'/tasmax_',dts,'.tif')
    txfls <- txfls[file.exists(txfls)]
    tnfls <- paste0(tn_pth,'/tasmin_',dts,'.tif')
    tnfls <- tnfls[file.exists(tnfls)]
    
    ## Read maximum temperature data
    tmx <- terra::rast(txfls)
    tmx <- terra::crop(tmx, xtd)
    
    ## Read minimum temperature data
    tmn <- terra::rast(tnfls)
    tmn <- terra::crop(tmn, xtd)
    
    #calculate average temp, clean-up
    tav <- (tmx + tmn)/2
    rm(tmx, tmn); gc(verbose=F, full=T, reset=T)
    tav <- mean(tav)
    terra::writeRaster(x = tav, filename = outfile, overwrite = T)
    rm(tav); gc(verbose=F, full=T, reset=T)
    
  }
  
}

# Future setup
# ssps  <- args[1]
# gcms  <- args[2]

# Runs
if (sce_climate == 'future') {
  
  ssps <- c('ssp126','ssp245','ssp370','ssp585') 
  gcms <- c('ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CMCC-ESM2','EC-Earth3','EC-Earth3-Veg-LR','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','TaiESM1')
  
  for(ssp in ssps){
    for(gcm in gcms){
      
      ## Parameters
      cmb <- paste0(ssp, '_', gcm)
      yrs <- 2021:2100
      mnt <- sprintf('%02.0f',1:12)
      stp <- base::expand.grid(yrs, mnt, stringsAsFactors = F) |> setNames(c('yrs', 'mnt')) |> dplyr::arrange(yrs, mnt) |> base::as.data.frame(); rm(yrs, mnt)
      
      ## Setup in/out files
      tn_pth <- paste0(root, '/nex-gddp-cmip6/tasmin/', ssp, '/', gcm) # Daily minimum temperatures
      tx_pth <- paste0(root, '/nex-gddp-cmip6/tasmax/', ssp, '/', gcm) # Daily maximum temperatures
      out_dir <- paste0(root, '/nex-gddp-cmip6_indices/', ssp, '_', gcm, '/TAVG/')
      
      1:nrow(stp) |> purrr::map(.f = function(i){calc_tav(yr = stp$yrs[i], mn = stp$mnt[i])})
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
      tn_pth <- paste0(root, '/nex-gddp-cmip6/tasmin/', ssp, '/', gcm) # Daily minimum temperatures
      tx_pth <- paste0(root, '/nex-gddp-cmip6/tasmax/', ssp, '/', gcm) # Daily maximum temperatures
      out_dir <- paste0(root, '/nex-gddp-cmip6_indices/', ssp, '_', gcm, '/TAVG/')
      
      1:nrow(stp) |> purrr::map(.f = function(i){calc_tav(yr = stp$yrs[i], mn = stp$mnt[i])})
      tmpfls <- list.files(tempdir(), full.names = T)
      1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ', tmpfls[k]))})
      ##
      cat('----Finish----\n')
      
    }
    
  }
}
