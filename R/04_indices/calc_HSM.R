#Calculate heat stress maize
#JRV, Dec. 2022

#load packages
library(terra)
library(tidyverse)
library(lubridate)

#clean-up environment
rm(list=ls())
gc(verbose=FALSE, full=TRUE, reset=TRUE)

#working directory
wd <- "~/common_data"

#reference raster
r_ref <- terra::rast(paste0(wd,'/atlas_hazards/roi/africa.tif'))

# Load daily maximum temperature
tx_dir <- paste0(wd, "/chirts/Tmax")

# Calculate NTx40 function
calc_hsm_thr <- function(yr, mn, thr){
  #yr <- 1995
  #mn <- 01
  #thr <- 35
  outfile <- paste0(wd, "/atlas_hazards/cmip6/indices/historic/HSM_NTx", thr, "/NTx", thr, "-", yr, "-", mn, ".tif")
  if (!file.exists(outfile)) {
    if (!file.exists(dirname(outfile))) dir.create(dirname(outfile), recursive=TRUE)
    
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    
    # Files
    fls <- paste0(tx_dir, '/', yr, '/Tmax.', gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    fls <- fls[file.exists(fls)]
    
    # Read raster data
    tmx <- terra::rast(fls)
    tmx <- tmx %>% terra::crop(terra::ext(r_ref)) %>% terra::mask(r_ref)
    tmx <- terra::app(x=tmx, fun=function(x) {rout <- x; rout[rout[] < -9990] <- NA; return(rout)})
    
    #get dates back from names
    tdates <- names(tmx)
    tdates <- gsub(pattern="Tmax.", replacement="", x=tdates)
    tdates <- gsub(pattern=".", replacement="-", x=tdates, fixed=TRUE)
    names(tmx) <- tdates
    
    #calculate stack with day number
    tmx_doy <- terra::app(x=tmx[[1:10]], fun=function(x) {
      rout <- x; rout[!is.na(x[])] <- yday(as.Date(names(x))); return(rout)
      })
    
    # Read crop calendar data
    r_cal <- terra::rast(paste0(wd, "/atlas_crop_calendar/intermediate/mai_rf_ggcmi_crop_calendar_phase3_v1.01_Africa.tif"))
    
    # Calculate heat stress generic crop
    terra::app(x   = tmx,
               fun = function(x){ ntx40 = sum(x > 40, na.rm = T); return(ntx40) },
               filename = outfile)
  }
}

# Setup
yrs <- 1995:2014
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()

1:nrow(stp) %>%
  purrr::map(.f = function(i){calc_ntx40(yr = stp$yrs[i], mn = stp$mns[i])})
