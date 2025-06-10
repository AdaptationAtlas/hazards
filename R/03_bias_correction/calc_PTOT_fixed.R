## Total precipitation per month
## By: H. Achicanoy / J. Ramirez-Villegas
## December, 2022

# R options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate))

root <- '/home/jovyan/common_data'

ref <- terra::rast(paste0(root,'/atlas_hazards/roi/africa.tif'))

# Calculate PTOT function
calc_ptot <- function(yr, mn){
  outfile <- paste0(out_dir,'/PTOT-',yr,'-',mn,'.tif')
  cat(outfile, "\n")
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
  prc <- prc |> terra::crop(terra::ext(ref)) |> terra::mask(ref)
  prc <- terra::classify(prc, rcl = cbind(-Inf,-9990, NA))
  prc <- terra::classify(prc, rcl = cbind(-Inf, 0, 0))
  # Calculate total monthly precipitation
  ptot <- sum(prc)
  terra::writeRaster(x = ptot, filename = outfile, overwrite = T)
  return(ptot)
}

# # # Historical setup
# yrs <- 1995:2014
# mns <- c(paste0('0',1:9),10:12)
# stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
# names(stp) <- c('yrs','mns')
# stp <- stp %>%
#   dplyr::arrange(yrs, mns) %>%
#   base::as.data.frame()
# stp <- stp |> dplyr::filter(mns == '02')
# pr_pth <- paste0(root,'/chirps_wrld') # Daily precipitation
# out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/historical/PTOT')
# 
# obs_clm <- 1:nrow(stp) %>%
#   purrr::map(.f = function(i){calc_ptot(yr = stp$yrs[i], mn = stp$mns[i])}) |>
#   terra::rast() |>
#   max()

###
# Future setup
cmb <- paste0(ssp,'_',gcm,'_',prd)
prd_num <- as.numeric(unlist(strsplit(x = prd, split = '_')))
yrs <- prd_num[1]:prd_num[2]
mns <- c(paste0('0',1:9),10:12)
stp <- base::expand.grid(yrs, mns) %>% base::as.data.frame(); rm(yrs,mns)
names(stp) <- c('yrs','mns')
stp <- stp |>
  dplyr::arrange(yrs, mns) |>
  base::as.data.frame()
stp <- stp |> dplyr::filter(mns == mn)
pr_pth <- paste0(root,'/chirps_cmip6_africa/Prec_',gcm,'_',ssp,'_',prd) # Precipitation
out_dir <- paste0(root,'/atlas_hazards/cmip6/indices/',cmb,'/PTOT')

1:nrow(stp) |> purrr::map(.f = function(i){calc_ptot(yr = stp$yrs[i], mn = stp$mns[i])})
