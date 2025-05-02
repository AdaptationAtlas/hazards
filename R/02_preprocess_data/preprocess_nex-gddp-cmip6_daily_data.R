# Nex-GDDP CMIP6 to Atlas data structure
# By: H. Achicanoy
# ABC, March 2025

options(warn = -1, scipen = 999)
if (!require(pacman)) {install.packages('pacman'); library(pacman)} else {library(pacman)}
pacman::p_load(tidyverse, terra, furrr, future, xts, tsbox)

root <- '/home/jovyan/common_data'

## Conversion factors to apply
# Precipitation: pr * 86400
# Temperatures: tasmax, tasmin - 273.15
# Solar radiation: rsds * 86400 / 1000000
# rotate

# Set up parallel processing
future::plan(future::multisession, workers = 6) # parallel::detectCores() - 1

get_daily_data <- function (vr, ssp, gcm) {
  
  # Input directory
  indir <- file.path(root,'nex-gddp-cmip6_raw',vr,ssp,gcm)
  # Output directory
  outdir <- file.path(root,'nex-gddp-cmip6',vr,ssp,gcm)
  dir.create(path = outdir, F, recursive = T)
  
  # Files in input directory
  fls <- list.files(path = indir, pattern = '.nc$', full.names = T)
  
  # Process files in parallel
  furrr::future_map(.x = fls, .f = function(fl) {
    
    # Read annual raster
    r <- terra::rast(fl)
    # Get daily dates
    dts <- as.character(terra::time(r))
    
    # Create output filenames
    out_files <- file.path(outdir, paste0(vr,'_',dts,'.tif'))
    
    # Check which files need to be processed
    to_process <- !file.exists(out_files)
    
    if (any(to_process)) {
      # Wrap raster for parallel processing
      r_wrapped <- terra::wrap(r)
      
      # Apply unit transformations
      r_unwrapped <- terra::unwrap(r_wrapped)
      if (vr == 'pr') {
        r_unwrapped <- r_unwrapped * 86400
      } else if (vr %in% c('tasmax','tasmin')) {
        r_unwrapped <- r_unwrapped - 273.15
      } else if (vr == 'rsds') {
        r_unwrapped <- r_unwrapped * 86400 / 1000000
      }
      
      # Rotate rasters
      r_unwrapped <- terra::rotate(r_unwrapped)
      
      # Write only the files that don't exist
      terra::writeRaster(x = r_unwrapped[[to_process]], filename = out_files[to_process], overwrite=T)
    }
    
  }, .progress = TRUE)
  
  return(cat('Done.\n'))
  
}

vrs <- c('pr','tasmax','tasmin','rsds','hurs') # Weather variables
ssps <- c('ssp126','ssp245','ssp370','ssp585') # SSP scenarios
gcms <- c('ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0') # GCM models

stp <- base::expand.grid(gcm = gcms, ssp = ssps, vr = vrs, stringsAsFactors = F) |>
  base::as.data.frame(); rm(vrs, ssps, gcms)

1:nrow(stp) |>
  purrr::map(.f = function(j) {
    vr  <- stp$vr[j]; ssp <- stp$ssp[j]; gcm <- stp$gcm[j]
    get_daily_data(vr = vr, ssp = ssp, gcm = gcm)
    cat(vr,ssp,gcm,'ready...\n')
  })
