## Number of waterlogging days (days) at start of saturation (NDWL0)
## By: H. Achicanoy & A. Mendez
## April, 2025

# R options
args <- commandArgs(trailingOnly = T)
options(warn = -1, scipen = 999) # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate,compiler,raster,ncdf4))

peest2 <- function(srad, tmin, tmean, tmax){
  
  # Convert rasters to matrices for faster processing
  srad_m  <- raster::values(srad) # this variable comes in .nc format
  tmin_m  <- terra::values(tmin)
  tmean_m <- terra::values(tmean)
  tmax_m  <- terra::values(tmax)
  
  # Constants
  albedo  <- 0.2
  vpd_cte <- 0.7
  a_eslope <- 611.2
  b_eslope <- 17.67
  c_eslope <- 243.5
  
  # Pre-allocate matrices
  n_cells <- ncol(srad_m)
  n_layers <- nrow(srad_m)
  rn <- matrix(0, nrow = n_layers, ncol = n_cells)
  eslope <- matrix(0, nrow = n_layers, ncol = n_cells)
  et_max <- matrix(0, nrow = n_layers, ncol = n_cells)
  
  # Optimized calculations
  rn <- (1-albedo) * srad_m
  rm(srad_m)
  
  temp_denom <- tmean_m + c_eslope
  eslope <- (a_eslope * b_eslope * c_eslope * exp(b_eslope * tmean_m/temp_denom)) / (temp_denom^2)
  rm(tmean_m, temp_denom)
  
  vpd <- vpd_cte * (0.61120 * (exp((17.67*tmax_m)/(tmax_m+243.5)) - 
                                 exp((17.67*tmin_m)/(tmin_m+243.5))))
  rm(tmin_m, tmax_m)
  
  # Constants
  pt_const <- 1.26
  pt_fact  <- 1
  vpd_ref  <- 1
  psycho   <- 62
  rho_w    <- 997
  rlat_ht  <- 2.26E6
  
  # Final calculation
  pt_coef <- 1 + (pt_fact*pt_const-1) * vpd / vpd_ref
  conversion_factor <- 1E6 * 100 / (rlat_ht * rho_w) * 10
  et_max <- pt_coef * rn * eslope/(eslope+psycho) * conversion_factor
  
  # Convert back to raster
  et_max_rast <- tmin
  terra::values(et_max_rast) <- et_max
  
  return(et_max_rast)
  
}

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

# Soil variables
scp <- terra::rast(paste0(root,'/atlas_hazards/soils/africa_scp.tif'))
scp <- scp |> terra::resample(ref) |> terra::mask(ref)
sst <- terra::rast(paste0(root,'/atlas_hazards/soils/africa_ssat.tif'))
sst <- sst |> terra::resample(ref) |> terra::mask(ref)

# Calculate NDWL0 function
calc_ndwl0 <- function(yr, mn){
  outfile <- paste0(out_dir,'/NDWL0-',yr,'-',mn,'.tif')
  cat(outfile,'\n')
  if (!file.exists(outfile)) {
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    if(as.numeric(yr) > 2020 & mn == '02'){
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-28')), by = 'day')
    } else {
      dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    }
    
    cat(">>> Iniciando proceso:", yr, "-", mn, " \n")
    
    pr_fls <- paste0(pr_pth,'/chirps-v2.0.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    pr_fls <- pr_fls[file.exists(pr_fls)]
    tx_fls <- paste0(tx_pth,'/',yr,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tm_fls <- paste0(tm_pth,'/',yr,'/Tmin.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    if(as.numeric(yr) > 2020){
      dts_chr <- as.character(dts)
      dts_chr <- gsub(pattern = yr, replacement = yrs_mpg$Baseline[yrs_mpg$Future == yr], x = dts_chr)
      dts <- as.Date(dts_chr); rm(dts_chr)
    }
    sr_fls <- paste0(sr_pth,'/Solar-Radiation-Flux_C3S-glob-agric_AgERA5_',gsub(pattern='-', replacement='', x=dts, fixed=T),'_final-v1.0.nc')
    sr_fls <- sr_fls[file.exists(sr_fls)]
    
    # Read variables
    prc <<- terra::rast(pr_fls)
    tmx <<- terra::rast(tx_fls)
    tmn <<- terra::rast(tm_fls)
    tav <<- (tmx + tmn)/2
    srd <- raster::stack(sr_fls)
    srd <- raster::crop(srd, raster::extent(raster::raster(ref) ))
    srd <- srd/1000000
    srd <- srd |> raster::resample(raster(ref)) |> raster::mask(raster(ref))
    
    # Maximum evapotranspiration
    ETMAX <<- peest2(srad = srd, tmin = tmn, tmean = tav, tmax = tmx)
    rm(list = c('tmn','tmx','tav','srd','sptd'))
    
    # Compute water balance model
    date <- paste0(yr,'-',mn)
    if(date %in% c('1995-01','2021-01','2041-01','2061-01','2081-01')){
      AVAIL <<- ref
      AVAIL[!is.na(AVAIL)] <- 0
    } else {
      AVAIL <<- terra::rast(paste0(dirname(outfile),'/AVAIL.tif'))
    }
    
    eabyep_calc <- function(soilcp = scp, soilsat = ssat, avail = AVAIL, rain = prc[[1]], evap = ETMAX[[1]]){
      
      avail   <- min(avail, soilcp)
      
      # ERATIO
      percwt <- min(avail/soilcp*100, 100)
      percwt <- max(percwt, 1)
      eratio <- min(percwt/(97-3.868*sqrt(soilcp)), 1)
      
      demand  <- eratio * evap
      result  <- avail + rain - demand
      logging <- result - soilcp
      logging <- max(logging, 0)
      logging <- min(logging, soilsat)
      # runoff  <- result - logging + soilcp
      avail   <- min(soilcp, result)
      avail   <- max(avail, 0)
      # runoff  <- max(runoff, 0)
      
      out     <- list(Availability = c(AVAIL, avail),
                      # Demand       = demand,
                      # Eratio       = eratio,
                      Logging      = logging
      )
      return(out)
    }
    ceabyep_calc <- compiler::cmpfun(eabyep_calc)
    
    watbal <- vector("list", terra::nlyr(ETMAX))
    for(j in 1:terra::nlyr(ETMAX)){
      water_balance <- ceabyep_calc(soilcp  = scp,
                                    soilsat = sst,
                                    avail   = AVAIL[[terra::nlyr(AVAIL)]],
                                    rain    = prc[[j]],
                                    evap    = ETMAX[[j]])
      # Update AVAIL with deep copy to avoid memory leaks
      AVAIL <<- terra::deepcopy(water_balance$Availability)
      # Store result and clean temporary objects
      watbal[[j]] <- water_balance
      rm(water_balance)
    }
    
    LOGGING <- watbal |> purrr::map('Logging') |> terra::rast()
    NDWL0  <- sum(LOGGING > 0)
    terra::writeRaster(NDWL0, outfile, overwrite = T)
    terra::writeRaster(AVAIL[[terra::nlyr(AVAIL)]], paste0(dirname(outfile),'/AVAIL.tif'), overwrite = T)
    
  }
}

ssps  <- args[1]
gcms  <- args[2]
prds  <- args[3]

#c('ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0')
# Rscript fast_calc_NDWL0.R 'ssp245' 'ACCESS-ESM1-5' 'EC-Earth3'                - screen id = ndwl0
# Rscript fast_calc_NDWL0.R 'ssp245' 'INM-CM5-0' 'MPI-ESM1-2-HR' 'MRI-ESM2-0'     screen id = ndwl0_2

#'2021_2040','2041_2060',
# Future setup
for (gcm in gcms ) {
  for (ssp in ssps) {
    for (prd in prds) {
      
      cmb <- paste0(ssp,'_',gcm,'_',prd)
      prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
      yrs <- prd_num[1]:prd_num[2]
      mns <- c(paste0('0',1:9),10:12)
      stp <- base::expand.grid(yrs, mns) |> base::as.data.frame(); rm(yrs,mns)
      names(stp) <- c('yrs','mns')
      stp <- stp |>
        dplyr::arrange(yrs, mns) |>
        base::as.data.frame()
      pr_pth <- paste0(root,'/chirps_cmip6_africa/Prec_',gcm,'_',ssp,'_',prd) # Precipitation
      tx_pth <- paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd) # Maximum temperature
      tm_pth <- paste0(root,'/chirts_cmip6_africa/Tmin_',gcm,'_',ssp,'_',prd) # Minimum temperature
      sr_pth <- paste0(root,'/ecmwf_agera5/solar_radiation_flux')             # Solar radiation
      out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/NDWL0')
      
      yrs_mpg <- data.frame(Baseline = as.character(rep(1995:2014, 4)),
                            Future = as.character(c(2021:2040,2041:2060,2061:2080,2081:2100)))
      
      1:nrow(stp) |>
        purrr::map(.f = function(i){
          calc_ndwl0(yr = stp$yrs[i], mn = stp$mns[i])
          if (i%%5 == 0) {
            tmpfls <- list.files(tempdir(), full.names=TRUE)
            1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
          }
        })
      
    }
  }
}
