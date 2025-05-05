## Heat stress generic crop and threshold (i.e., NTx40)
## By: H. Achicanoy, F. Castro
## May, 2025

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

# Get CHIRPS extent
msk <- terra::rast('/home/jovyan/common_data/chirps_wrld/chirps-v2.0.1981.01.01.tif')
xtd <- terra::ext(msk); rm(msk)

# Calculate NDD function
calc_ntx <- function(yr, mn, thr = 40){
  outfile <- paste0(out_dir,'/NTx',thr,'/NTx',thr,'-',yr,'-',mn,'.tif') 
  thr <- thr[!file.exists(outfile)]
  outfile <- outfile[!file.exists(outfile)]
  if (length(outfile) > 0) {
    cat('...processing n=', length(outfile), 'files for yr=', yr, '/ mn=', mn, '\n')
    # Create directories
    1:length(outfile) |> purrr::map(.f = function(j){dir.create(dirname(outfile[j]),F,T)})
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    # Files
    fls <- paste0(tx_pth,'/tasmax','_',dts,'.tif')
    fls <- fls[file.exists(fls)]
    # Read maximum temperature data
    tmx <- terra::rast(fls)
    tmx <- terra::crop(tmx, xtd)
    # Calculate heat stress generic crop
    for (j in 1:length(thr)) {
      cat('...processing threshold thr=',thr[j],'\n')
      ntx <- sum(tmx > thr[j])
      terra::writeRaster(x = ntx, filename = outfile[j])
    }
    rm(tmx); gc(F, T, T) # Clean up
  }
}


# Future setup
for(ssp in c('ssp126', 'ssp245', 'ssp370', 'ssp585')){
  for(gcm in c('ACCESS-ESM1-5', 'MPI-ESM1-2-HR', 'EC-Earth3', 'INM-CM5-0', 'MRI-ESM2-0')){
    
    
    ## To start the process
    cat('----------------------', ssp, ' ', gcm, '--------------------------\n')
    
    ## Parameters
    cmb <- paste0(ssp, '_', gcm)
    yrs <- 2021:2100
    mnt <- c(paste0('0', 1:9), 10:12)
    stp <- base::expand.grid(yrs, mnt, stringsAsFactors = F) |> base::as.data.frame() |> setNames(c('yrs', 'mnt')) |> dplyr::arrange(yrs, mnt) |> base::as.data.frame(); rm(yrs, mnt)
    
    ## Setup in/out files
    tx_pth  <- paste0(root, '/nex-gddp-cmip6/tasmax/', ssp, '/', gcm)
    thr = 30:50
    out_dir <- paste0(root, '/nex-gddp-cmip6_indices/', ssp, '_', gcm) 
    
    1:nrow(stp) |>
      purrr::map(.f = function(i){calc_ntx(yr = stp$yrs[i], mn = stp$mnt[i], thr = 30)})
    
    ##
    cat('----Finish----\n')
    
  }
}
