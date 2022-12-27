## Number of dry days
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

# Calculate NDD function
calc_ndd <- function(yr, mn){
  outfile <- paste0(out_dir,'/NDD-',yr,'-',mn,'.tif')
  if(!file.exists(outfile)){
    dir.create(dirname(outfile),F,T)
    # Last day of the month
    last_day <- lubridate::days_in_month(as.Date(paste0(yr,'-',mn,'-01')))
    # Sequence of dates
    dts <- seq(from = as.Date(paste0(yr,'-',mn,'-01')), to = as.Date(paste0(yr,'-',mn,'-',last_day)), by = 'day')
    # Files
    fls <- paste0(pr_pth,'/chirps-v2.0.',gsub(pattern='-', replacement='.', x=dts, fixed=T),'.tif')
    fls <- fls[file.exists(fls)]
    # Read precipitation data
    prc <- terra::rast(fls)
    prc <- prc %>% terra::crop(terra::ext(ref)) %>% terra::mask(ref)
    prc[prc == -9999] <- NA
    # Calculate number of dry days
    terra::app(x   = prc,
               fun = function(x){ ndd = sum(x < 1, na.rm = T); return(ndd) },
               filename = outfile)
  }
}

# # Historical setup
# yrs <- 1995:2014
# mns <- c(paste0('0',1:9),10:12)
# stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
# names(stp) <- c('yrs','mns')
# stp <- stp %>%
#   dplyr::arrange(yrs, mns) %>%
#   base::as.data.frame()
# pr_pth <- paste0(root,'/chirps_wrld') # Daily precipitation
# out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/historical/NDD')

# Future setup
gcm <- 'EC-Earth3'
ssp <- 'ssp585'
prd <- '2021_2040'

cmb <- paste0(ssp,'_',gcm,'_',prd)
prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
yrs <- prd_num[1]:prd_num[2]
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp %>%
  dplyr::arrange(yrs, mns) %>%
  base::as.data.frame()
pr_pth <- paste0(root,'/chirps_cmip6_africa/Prec_',gcm,'_',ssp,'_',prd) # Precipitation
out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/NDD')

1:nrow(stp) %>%
  purrr::map(.f = function(i){calc_ndd(yr = stp$yrs[i], mn = stp$mns[i])})
