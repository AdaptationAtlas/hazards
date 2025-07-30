# Compute human heat stress (HSH)
# By: H. Achicanoy
# Alliance Bioversity International & CIAT, 2025

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

# Extent CHIRPS
msk <- terra::rast(file.path(root,'chirps_wrld/chirps-v2.0.1981.01.01.tif'))
xtd <- terra::ext(msk)

# HSH function
calc_hsh <- function(yr, mn){
  
  outfile1 <- paste0(out_dir,'/HSH_mean-',yr,'-',mn,'.tif')
  outfile2 <- paste0(out_dir,'/HSH_max-',yr,'-',mn,'.tif')
  outfile3 <- paste0(out_dir,'/daily/HSH_daily-',yr,'-',mn,'.tif')
  
  cat(outfile2,'\n')
  
  if(sum(file.exists(c(outfile1,outfile2,outfile3))) < 3){
    
    dir.create(dirname(outfile3),F,T)
    
    # Sequence of dates
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01'))) # Last day of the month
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    
    # Files
    tx_fls <- paste0(tx_pth,'/tasmax_',dts,'.tif')
    tx_fls <- tx_fls[file.exists(tx_fls)]
    tn_fls <- paste0(tn_pth,'/tasmin_',dts,'.tif')
    tn_fls <- tn_fls[file.exists(tn_fls)]
    rh_fls <- paste0(rh_pth,'/hurs_',dts,'.tif')
    rh_fls <- rh_fls[file.exists(rh_fls)]
    
    # Read daily minimum and maximum temperatures, and relative humidity data
    tmx <- terra::rast(tx_fls) |> terra::crop(xtd)
    tmn <- terra::rast(tn_fls) |> terra::crop(xtd)
    rhm <- terra::rast(rh_fls) |> terra::crop(xtd)
    tav <- (tmx + tmn)/2
    
    # Constants
    c1 = -8.78469475556
    c2 =  1.61139411
    c3 =  2.33854883889
    c4 = -0.14611605
    c5 = -0.012308094
    c6 = -0.0164248277778
    c7 =  2.211732 * 10^(-3)
    c8 =  7.2546 * 10^(-4)
    c9 = -3.582 * 10^(-6)
    # Compute human stress index
    heat_idx <- function(tmean, rhum){
      hi <- terra::ifel(test = tmean >= 25,
                        yes  = c1 + (c2*tmean) + (c3*rhum) + (c4*tmean*rhum) + (c5*tmean^2) + (c6*rhum^2) + (c7*tmean^2*rhum) + (c8*tmean*rhum^2) + (c9*tmean^2*rhum^2),
                        no   = tmean)
      return(hi)
    }
    cheat_idx <- compiler::cmpfun(heat_idx)
    HI <- cheat_idx(tmean = tav, rhum = rhm)
    HI_avg <- mean(HI)
    HI_max <- max(HI)
    
    terra::writeRaster(HI_avg, outfile1, overwrite = T)
    terra::writeRaster(HI_max, outfile2, overwrite = T)
    terra::writeRaster(HI, outfile3, overwrite = T)
    
    # Clean-up
    rm(tmx, tmn, tav, rhm, HI, HI_avg, HI_max); gc(F, T, T)
    
  }
}

# Runs
scenario <- 'future' # historical, future
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
    stp <- base::expand.grid(yrs, mnt, stringsAsFactors = F) |> setNames(c('yrs', 'mnt')) |> dplyr::arrange(yrs, mnt) |> base::as.data.frame(); rm(yrs, mnt)
    
    ## Setup in/out files
    tx_pth <- paste0(root,'/nex-gddp-cmip6/tasmax/',ssp,'/',gcm) # Daily maximum temperatures
    tn_pth <- paste0(root,'/nex-gddp-cmip6/tasmin/',ssp,'/',gcm) # Daily minimum temperatures
    rh_pth <- paste0(root,'/nex-gddp-cmip6/hurs/',ssp,'/',gcm)   # Daily relative humidity
    out_dir <- paste0(root,'/nex-gddp-cmip6_indices/',ssp,'_',gcm,'/HSH')
    
    1:nrow(stp) |> purrr::map(.f = function(i){calc_hsh(yr = stp$yrs[i], mn = stp$mnt[i])}); gc(F, T, T)
    tmpfls <- list.files(tempdir(), full.names = T)
    1:length(tmpfls) |> purrr::map(.f = function(k) {system(paste0('rm -f ', tmpfls[k]))})
    cat('----Finish----\n')
    
  }
  
}
