# Nex-GDDP-CMIP6 to Atlas data structure
# By: H. Achicanoy
# Alliance Bioversity-International & CIAT, 2025

# R options
options(warn = -1, scipen = 999)
if (!require(pacman)) {install.packages('pacman'); library(pacman)} else {library(pacman)}
pacman::p_load(tidyverse, terra, furrr, future, xts, tsbox)

grep2 <- Vectorize(grep, 'pattern')

# Root directory
root <- '/home/jovyan/common_data'

## Conversion factors
# Precipitation: pr * 86400
# Temperatures: tasmax, tasmin - 273.15
# Solar radiation: rsds * 86400 / 1000000
# Rotate

# Get daily transformed data
get_daily_data <- function (vr, ssp, gcm) {
  
  # Input directory
  indir <- file.path(root,'nex-gddp-cmip6_raw',vr,ssp,gcm)
  # Output directory
  outdir <- file.path(root,'nex-gddp-cmip6',vr,ssp,gcm)
  dir.create(path = outdir, F, recursive = T)
  
  # Files in input directory
  fls <- list.files(path = indir, pattern = '.nc$', full.names = T)
  
  # Process files in parallel
  # Set up parallel processing
  future::plan(future::multisession, workers = 6) # parallel::detectCores() - 1
  furrr::future_map(.x = fls, .f = function(fl) {
    
    # Read annual raster
    r <- terra::rast(fl)
    # Get daily dates
    if (gcm == 'KACE-1-0-G') {
      
      dts <- as.character(terra::time(r))
      dts <- dts[-which(duplicated(dts))] # Remove duplicated dates
      yr  <- unique(lubridate::year(dts))
      leap_yr <- grep(pattern = paste0(yr,'-02-29'), x = dts)
      if (length(leap_yr) > 0) {dts <- dts[-leap_yr]} # Remove Feb 29 if leap year
      r <- r[[match(x = dts, table = as.character(time(r)))]]
      
    } else {
      
      dts <- as.character(terra::time(r))
      
    }
    
    # Create output filenames
    out_files <- file.path(outdir, paste0(vr,'_',dts,'.tif'))
    
    # Check which files need to be processed
    to_process <- !file.exists(out_files)
    
    if (any(to_process)) {
      # Wrap raster for parallel processing
      r_wrapped <- terra::wrap(r)
      
      # Apply unit transformations
      r_unwrapped <- terra::unwrap(r_wrapped)
      if (vr == 'pr2') {
        r_unwrapped <- r_unwrapped * 86400
      } else if (vr %in% c('tasmax','tasmin')) {
        r_unwrapped <- r_unwrapped - 273.15
      } else if (vr == 'rsds') {
        r_unwrapped <- r_unwrapped * 86400 / 1000000
      }
      
      # Rotate rasters
      r_unwrapped <- terra::rotate(r_unwrapped)
      
      # Write only the files that don't exist
      terra::writeRaster(x = r_unwrapped[[to_process]], filename = out_files[to_process], overwrite = T)
    }
    
  }, .progress = T)
  future::plan(future::sequential)
  gc(F, T, T)
  
  return(cat('Done.\n'))
  
}

scenario <- 'historical'

if (scenario == 'future'){
  ssps <- c('ssp126','ssp245','ssp370','ssp585')
} else {
  if (scenario == 'historical') {
    ssps <- 'historical'
  }
}
vrs <- c('pr') # ,'hurs' 'tasmax','tasmin','rsds'
gcms <- c('ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CMCC-ESM2','EC-Earth3','EC-Earth3-Veg-LR','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','TaiESM1')

stp <- base::expand.grid(gcm = gcms, ssp = ssps, vr = vrs, stringsAsFactors = F) |>
  base::as.data.frame(); rm(vrs, ssps, gcms)

1:nrow(stp) |>
  purrr::map(.f = function(j) {
    vr  <- paste0(stp$vr[j],'2'); ssp <- stp$ssp[j]; gcm <- stp$gcm[j]
    get_daily_data(vr = vr, ssp = ssp, gcm = gcm)
    cat(vr,ssp,gcm,'ready...\n')
  })
