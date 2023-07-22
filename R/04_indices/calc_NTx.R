## Heat stress generic crop (NTx40)
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

sce_climate <- "future" #"future"

# Calculate NTx function
calc_ntx <- function(yr, mn, thr=40){
  outfile <- paste0(out_dir,'/NTx',thr,'/NTx',thr,'-',yr,'-',mn,'.tif') 
  thr <- thr[!file.exists(outfile)]
  outfile <- outfile[!file.exists(outfile)]
  if(length(outfile) > 0){
    cat("...processing n=", length(outfile), "files for yr=", yr, "/ mn=", mn, "\n")
    # Create directories
    1:length(outfile) %>%
      purrr::map(.f = function(j){dir.create(dirname(outfile[j]),F,T)})
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    # Files
    fls <- paste0(tx_pth,'/',yr,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    fls <- fls[file.exists(fls)]
    # Read maximum temperature data
    tmx <- terra::rast(fls)
    tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    tmx[tmx == -9999] <- NA
    # Calculate heat stress generic crop
    for (j in 1:length(thr)) {
      cat("...processing threshold thr=", thr[j], "\n")
      terra::app(x   = tmx,
                 fun = function(x){ ntxval = sum(x >= thr[j], na.rm = T); return(ntxval) },
                 filename = outfile[j])
    }
    # Clean up
    rm(tmx)
    gc(verbose=FALSE, full=TRUE, reset=TRUE)
  }
}

#perform calculation
if (sce_climate == "historical") {
  # Historical setup
  yrs <- 1995:2014
  mns <- c(paste0('0',1:9),10:12)
  stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
  names(stp) <- c('yrs','mns')
  stp <- stp %>%
    dplyr::arrange(yrs, mns) %>%
    base::as.data.frame()
  tx_pth <- paste0(root,'/chirts/Tmax') # Daily maximum temperature
  out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/historical')
  1:nrow(stp) %>%
    purrr::map(.f = function(i){
        calc_ntx(yr = stp$yrs[i], mn = stp$mns[i], thr=30:50)
        tmpfls <- list.files(tempdir(), full.names=TRUE)
        1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
      })
} else if (sce_climate == "future") {
  # Future setup
  gcm <- 'MRI-ESM2-0' #'ACCESS-ESM1-5' 'MPI-ESM1-2-HR' 'EC-Earth3' 'INM-CM5-0' 'MRI-ESM2-0'
  for (ssp in c('ssp245', 'ssp585')) {
      for (prd in c('2021_2040', '2041_2060')) {
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
          tx_pth <- paste0(root,'/chirts_cmip6_africa/Tmax_',gcm,'_',ssp,'_',prd) # Daily maximum temperatures
          out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb)
  
          1:nrow(stp) %>%
            purrr::map(.f = function(i){
              calc_ntx(yr = stp$yrs[i], mn = stp$mns[i], thr=30:50)
              tmpfls <- list.files(tempdir(), full.names=TRUE)
              1:length(tmpfls) %>% purrr::map(.f = function(k) {system(paste0("rm -f ", tmpfls[k]))})
            })
      }
  }
} else {
  cat("select one of historical or future for sce_climate \n")
}
