## Average monthly minimum temperature
## By: H. Achicanoy / Julian Ramirez-Villegas
## July, 2023

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

sce_climate <- "historical" #"future"

# Calculate tmin function
calc_tmn <- function(yr, mn){
  outfile <- paste0(out_dir,'/TMIN/TMIN-',yr,'-',mn,'.tif')
  if(!file.exists(outfile)){
    cat(outfile, "\n")
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    
    # Files
    tnfls <- paste0(tn_pth,'/',yr,'/Tmin.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    tnfls <- tnfls[file.exists(tnfls)]
    
    # Read temperature data
    tmn <- terra::rast(tnfls)
    tmn <- tmn %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tmn[tmn <= -9990] <- NA
    
    # Calculate heat stress generic crop
    terra::app(x   = tmn,
               fun = function(x){ tmin = mean(x, na.rm = T); return(tmin) },
               filename = outfile)
    # Clean up
    rm(tmn)
    gc(verbose=FALSE, full=TRUE, reset=TRUE)
  }
}

if (sce_climate == "historical") {
  # Historical setup
  yrs <- 1995:2014
  mns <- c(paste0('0',1:9),10:12)
  stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
  names(stp) <- c('yrs','mns')
  stp <- stp %>%
    dplyr::arrange(yrs, mns) %>%
    base::as.data.frame()
  tn_pth <- paste0(root,'/chirts/Tmin') # Daily minimum temperature
  out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/historical')
  1:nrow(stp) %>%
    purrr::map(.f = function(i){
      calc_tmn(yr = stp$yrs[i], mn = stp$mns[i])
      tmpfls <- list.files(tempdir(), full.names=TRUE)
      1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
    })
} else if (sce_climate == "future") {
  # Future setup
  gcm <- 'MRI-ESM2-0' #'ACCESS-ESM1-5' 'MPI-ESM1-2-HR' 'EC-Earth3' 'INM-CM5-0' 'MRI-ESM2-0'
  for (ssp in c('ssp245', 'ssp585')) { #'ssp126' 'ssp370'
    for (prd in c('2021_2040', '2041_2060')) { #'2061_2080', '2081_2100'
      cat("...processing gcm=", gcm, "/ ssp=", ssp, "/ period=", prd, "\n")
      #ssp <- 'ssp245' #ssp585
      #prd <- '2041_2060' #2021_2040
      
      cmb <- paste0(ssp,'_',gcm,'_',prd)
      prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
      yrs <- prd_num[1]:prd_num[2]
      mns <- c(paste0('0',1:9),10:12)
      stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
      names(stp) <- c('yrs','mns')
      stp <- stp %>%
        dplyr::arrange(yrs, mns) %>%
        base::as.data.frame()
      tn_pth <- paste0(root,'/chirts_cmip6_africa/Tmin_',gcm,'_',ssp,'_',prd) # Daily minimum temperatures
      out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb)
      
      1:nrow(stp) %>%
        purrr::map(.f = function(i){
          calc_tmn(yr = stp$yrs[i], mn = stp$mns[i])
          tmpfls <- list.files(tempdir(), full.names=TRUE)
          1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
        })
    }
  }
} else {
  cat("select one of historical or future for sce_climate \n")
}

