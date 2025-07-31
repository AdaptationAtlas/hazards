# Heat stress generic crop and threshold (i.e., NTx40)
# By: H. Achicanoy, F. Castro-Llanos
# Alliance Bioversity International & CIAT, 2025

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

# Get CHIRPS extent
msk <- terra::rast('/home/jovyan/common_data/chirps_wrld/chirps-v2.0.1981.01.01.tif')
xtd <- terra::ext(msk); rm(msk)

# NTx function
calc_ntx <- function(yr, mn, thr = 40) {
  
  outfile <- paste0(out_dir,'/NTx',thr,'/NTx',thr,'-',yr,'-',mn,'.tif') 
  thr <- thr[!file.exists(outfile)]
  outfile <- outfile[!file.exists(outfile)]
  
  if (length(outfile) > 0) {
    
    cat('...processing n=', length(outfile), 'files for yr=', yr, '/ mn=', mn, '\n')
    
    # Create directories
    1:length(outfile) |> purrr::map(.f = function(j){dir.create(dirname(outfile[j]),F,T)})
    
    # Sequence of dates
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    
    # Files
    fls <- paste0(tx_pth,'/tasmax','_',dts,'.tif')
    fls <- fls[file.exists(fls)]
    
    # Read daily maximum temperature data
    tmx <- terra::rast(fls)
    tmx <- terra::crop(tmx, xtd)
    
    # Calculate heat stress generic crop
    for (j in 1:length(thr)) {
      cat('...processing threshold thr=',thr[j],'\n')
      ntx <- sum(tmx > thr[j])
      terra::writeRaster(x = ntx, filename = outfile[j], overwrite = T)
    }
    
    # Clean-up
    rm(tmx, ntx); gc(F, T, T)
    
  }
  
}

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
    stp <- base::expand.grid(yrs, mnt, stringsAsFactors = F) |> base::as.data.frame() |> setNames(c('yrs', 'mnt')) |> dplyr::arrange(yrs, mnt) |> base::as.data.frame(); rm(mnt)
    
    ## Setup in/out files
    tx_pth  <- paste0(root,'/nex-gddp-cmip6/tasmax/',ssp,'/',gcm)
    thr <- 20:50
    out_dir <- paste0(root,'/nex-gddp-cmip6_indices/',ssp,'_',gcm) 
    
    1:nrow(stp) |> purrr::map(.f = function(i){calc_ntx(yr = stp$yrs[i], mn = stp$mnt[i], thr = thr)}); gc(F, T, T)
    tmpfls <- list.files(tempdir(), full.names = T)
    1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ', tmpfls[k]))})
    cat('----Finish----\n')
    
  }
  
}
