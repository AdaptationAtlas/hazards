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

# Load daily maximum temperature
tx_pth <- paste0(root,'/chirts/Tmax')

# Calculate NTx40 function
calc_ntx40 <- function(yr, mn){
  outfile <- paste0(root,'/atlas_hazards/indices/historic/NTx40/NTx40-',yr,'-',mn,'.tif')
  if(!file.exists(outfile)){
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    # Files
    fls <- paste0(tx_pth,'/',yr,'/Tmax.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    fls <- fls[file.exists(fls)]
    # Read precipitation data
    tmx <- terra::rast(fls)
    tmx <- tmx %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
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
