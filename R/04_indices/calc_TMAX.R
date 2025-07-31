# Compute maximum monthly temperature (TMAX)
# By: H. Achicanoy / J. Ramirez-Villegas / F. Castro-Llanos
# Alliance Bioversity-International & CIAT, 2025

# R options
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

# TMAX function
calc_tmax <- function(yr, mn){
  
  outfile <- paste0(out_dir,'/TMAX-',yr,'-',mn,'.tif')
  cat(outfile,'\n')
  
  if(!file.exists(outfile)){
    
    dir.create(dirname(outfile), F, T)
    
    # Sequence of dates
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    
    # Files
    txfls <- paste0(tx_pth,'/tasmax_',dts,'.tif')
    txfls <- txfls[file.exists(txfls)]
    
    ## Read daily maximum temperature data
    tmx <- terra::rast(txfls)
    tmx <- terra::crop(tmx, xtd)
    
    # Calculate maximum temperature
    tmax <- mean(tmx)
    terra::writeRaster(x = tmax, filename = outfile, overwrite = T)
    
    # Clean-up
    rm(tmx, tmax); gc(F, T, T)
    
  }
  
}

# Future setup
# ssps  <- args[1]
# gcms  <- args[2]

# Runs
scenario <- 'historical' # historical, future
if (scenario == 'future') {
  ssps <- c('ssp126', 'ssp245', 'ssp370', 'ssp585')
  yrs <- 2021:2100
} else {
  if (scenario == 'historical') {
    ssps <- 'historical'
    yrs <- 1995:2014
  }
}
gcms <- c('ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CMCC-ESM2','EC-Earth3','EC-Earth3-Veg-LR','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','TaiESM1')

for (gcm in gcms) {
  
  for (ssp in ssps) {
    
    ## Parameters
    cmb <- paste0(ssp, '_', gcm)
    mnt <- sprintf('%02.0f',1:12)
    stp <- base::expand.grid(yrs, mnt, stringsAsFactors = F) |> setNames(c('yrs', 'mnt')) |> dplyr::arrange(yrs, mnt) |> base::as.data.frame(); rm(mnt)
    
    ## Setup in/out files
    tx_pth <- paste0(root, '/nex-gddp-cmip6/tasmax/', ssp, '/', gcm) # Daily maximum temperatures
    out_dir <- paste0(root, '/nex-gddp-cmip6_indices/', ssp, '_', gcm, '/TMAX')
    
    1:nrow(stp) |> purrr::map(.f = function(i){calc_tmax(yr = stp$yrs[i], mn = stp$mnt[i])}); gc(F, T, T)
    tmpfls <- list.files(tempdir(), full.names = T)
    1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ', tmpfls[k]))})
    cat('----Finish----\n')
    
  }
  
}
