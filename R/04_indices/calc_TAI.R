# Compute Thornthwaite's Aridity Index (TAI)
# By: H. Achicanoy, F. Castro
# Alliance Bioversity International & CIAT, 2025

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate,envirem))

root <- '/home/jovyan/common_data'

# Get CHIRPS extent
msk <- terra::rast('/home/jovyan/common_data/chirps_wrld/chirps-v2.0.1981.01.01.tif')
xtd <- terra::ext(msk); rm(msk)

# ET SRAD, load and process only once
srf <- list.dirs(paste0(root,'/ET_SolRad'), full.names = T, recursive = F)
srf <- srf[-length(srf)]
srf <- srf |> gtools::mixedsort()
srd <- srf |> terra::rast()
names(srd) <- c(paste0('SRAD_0',1:9),paste0('SRAD_', 10:12))
gc(F, T, T)

# TAI function
calc_tai <- function(yr){
  outfile <- paste0(out_dir,'/TAI-',yr,'.tif')
  cat(outfile,'\n')
  if (!file.exists(outfile)) {
    dir.create(dirname(outfile),F,T)
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-01-01')), to = as.Date(paste0(yr,'-12-31')), by = 'day')
    # Files
    pr_fls <- paste0(pr_pth,'/pr_', dts,'.tif')
    pr_fls <- pr_fls[file.exists(pr_fls)]
    tx_fls <- paste0(tx_pth, '/tasmax_', dts, '.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tm_fls <- paste0(tm_pth, '/tasmin_', dts, '.tif')
    tm_fls <- tm_fls[file.exists(tm_fls)]
    # Read variables, and calculate monthly summaries, do by month to reduce memory consumption
    prc_ls <- tav_ls <- rng_ls <- c()
    for (j in 1:12) {
      cat('month=',j,'\n')
      # Load precipitation
      mth_fls <- paste0(pr_pth, '/pr_', yr, '-', sprintf('%02.0f', j), '-', sprintf('%02.0f', 1:31), '.tif')
      prc <- terra::rast(pr_fls[pr_fls %in% mth_fls])
      prc_month <- sum(prc)
      prc_month <- prc_month |> terra::crop(xtd); rm(prc); gc(F, T, T)
      
      # Load maximum temperature
      mth_fls <- paste0(tx_pth, '/tasmax_', yr, '-', sprintf('%02.0f', j), '-', sprintf('%02.0f', 1:31), '.tif')
      tmx <- terra::rast(tx_fls[tx_fls %in% mth_fls])
      ## Load minimum temperature
      mth_fls <- paste0(tm_pth, '/tasmin_', yr, '-', sprintf('%02.0f', j), '-', sprintf('%02.0f', 1:31), '.tif')
      tmn <- terra::rast(tm_fls[tm_fls %in% mth_fls])
      
      ## Calculate average temperature 
      tav <- (tmx + tmn) / 2
      tav_month <- mean(tav)
      tav_month <- tav_month |> terra::crop(xtd); rm(tav); gc(F, T, T)
      
      ## Calculate temperature range
      rnge <- abs(tmx - tmn)
      rng_month <- mean(rnge)
      rng_month <- rng_month |> terra::crop(xtd); rm(rnge, tmx, tmn); gc(F, T, T)
      
      # Append
      prc_ls <- c(prc_ls, prc_month)
      tav_ls <- c(tav_ls, tav_month)
      rng_ls <- c(rng_ls, rng_month)
    }; gc(F, T, T)
    
    # Assign precipitation names in envirem environment
    envirem::assignNames(solrad='SRAD_##', tmean = 'TMEAN_##', precip = 'PREC_##')
    
    TMEAN <- tav_ls |> terra::rast()
    PREC  <- prc_ls |> terra::rast()
    TRNG  <- rng_ls |> terra::rast()
    
    names(TMEAN) <- c(paste0('TMEAN_0',1:9),paste0('TMEAN_', 10:12))
    names(PREC)  <- c(paste0('PREC_0',1:9),paste0('PREC_', 10:12))
    names(TRNG)  <- c(paste0('TRNG_0',1:9),paste0('TRNG_', 10:12))
    
    # Clean up
    rm(tav_ls, prc_ls, rng_ls, tav_month, prc_month, rng_month)
    gc(F, T, T)
    
    # Resample extraterrestrial solar radiation to 0.25
    srd <- srd |> terra::resample(TMEAN[[1]])
    srd <- srd |> terra::crop(xtd)
    gc(F, T, T)
    
    # Thornthwaite's Aridity Index
    PET <- envirem::monthlyPET(Tmean = TMEAN, RA = srd, TD = TRNG)
    names(PET)  <- c(paste0('PET_0',1:9), paste0('PET_', 10:12))
    TAI <- envirem::aridityIndexThornthwaite(precipStack = PREC, PETstack = PET)
    names(TAI) <- yr
    
    # Write output
    terra::writeRaster(x = TAI, filename = outfile)
    
    # Clean up
    rm(PET, TMEAN, PREC, TRNG, TAI)
    gc(F, T, T)
    
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
    
    cmb <- paste0(ssp, '_', gcm)
    cat(cmb, '\n')
    
    pr_pth <- paste0(root,'/nex-gddp-cmip6/pr/',ssp,'/',gcm) # Precipitation
    tm_pth <- paste0(root,'/nex-gddp-cmip6/tasmin/',ssp,'/',gcm) # Minimum temperature
    tx_pth <- paste0(root,'/nex-gddp-cmip6/tasmax/',ssp,'/',gcm) # Maximum temperature
    out_dir <- paste0(root,'/nex-gddp-cmip6_indices/',cmb,'/TAI')
    
    yrs |> purrr::map(.f = function(yr){ calc_tai(yr = yr) }); gc(F, T, T)
    tmpfls <- list.files(tempdir(), full.names = T)
    1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ', tmpfls[k]))})
    
  }
  
}
