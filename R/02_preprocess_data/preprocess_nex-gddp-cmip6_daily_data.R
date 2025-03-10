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

get_daily_data <- function (vr, ssp, gcm) {
  
  # Input directory
  indir <- file.path(root,'nex-gddp-cmip6_raw',vr,ssp,gcm)
  # Output directory
  outdir <- file.path(root,'nex-gddp-cmip6',vr,ssp,gcm)
  dir.create(path = outdir, F, recursive = T)
  
  # Files in input directory
  fls <- list.files(path = indir, pattern = '.nc$', full.names = T)
  
  1:length(fls) |>
    purrr::map(.f = function(i) {
      
      # Read annual raster
      r <- terra::rast(fls[i])
      # Get daily dates
      dts <- as.character(terra::time(r))
      
      if (!all(file.exists(paste0(outdir,'/',vr,'_',dts,'.tif')))) {
        # Apply unit transformations
        if (vr == 'pr') { # to mm/day
          r <- r * 86400
        } else {
          if (vr %in% c('tasmax','tasmin')) { # to Celsius degrees
            r <- r - 273.15
          } else {
            if (vr == 'rsds') { # to MJ/m^2/day (megajoules per square meter per day)
              r <- r * 86400 / 1000000
            }
          }
        }
        # Rotate rasters
        r <- terra::rotate(r)
        terra::writeRaster(r, filename = paste0(outdir,'/',vr,'_',dts,'.tif'))
      }
      
    })
  
  return(cat('Done.\n'))
  
}

vrs <- c('pr','tasmax','tasmin','rsds','hurs') # Weather variables
ssps <- c('ssp126','ssp245','ssp370','ssp585') # SSP scenarios
gcms <- c('ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0') # GCM models

stp <- base::expand.grid(vr = vrs, ssp = ssps, gcm = gcms, stringsAsFactors = F) |>
  base::as.data.frame() |>
  dplyr::arrange(vr,ssp,gcm); rm(vrs, ssps, gcms)

1:nrow(stp) |>
  purrr::map(.f = function(j) {
    vr  <- stp$vr[j]; ssp <- stp$ssp[j]; gcm <- stp$gcm[j]
    get_daily_data(vr = vr, ssp = ssp, gcm = gcm)
    cat(vr,ssp,gcm,'ready...\n')
  })
