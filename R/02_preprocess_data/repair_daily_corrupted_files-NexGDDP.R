# Repair daily corrupted files Nex-GDDP-CMIP6
# By: H. Achicanoy
# Alliance Bioversity International & CIAT, 2025

options(warn = -1, scipen = 999)
library(pacman)
pacman::p_load(terra)

# corrupted_file <- '/home/jovyan/common_data/nex-gddp-cmip6/hurs/ssp126/NorESM2-MM/hurs_2072-08-29.tif'
# evl <- tryCatch(expr = {terra::rast(corrupted_file)}, error = function(e) {msg <- 'corrupted'})
# tst <- evaluate::evaluate('rhm <- terra::rast(rh_fls) |> terra::crop(xtd)')
# fle <- strsplit(x = as.character(tst[[2]]), split = 'cannot open this file as a SpatRaster: ')[[1]][2]
# fle <- gsub('\n', '', fle, fixed = T)

# The corrupted file is manually identified during the
# computation of monthly indices. Please copy and paste
# the path of the corrupted file. This script goes to
# the raw data and rewrite it again applying the proper
# transformations.

# Corrupted file
fle <- '/home/jovyan/common_data/nex-gddp-cmip6/hurs/ssp126/NorESM2-MM/hurs_2072-08-29.tif'

crr_ssp <- strsplit(fle, split = '/') |> purrr::map(7) |> unlist()
crr_gcm <- strsplit(fle, split = '/') |> purrr::map(8) |> unlist()
crr_dte <- strsplit(fle, split = '/') |> purrr::map(9) |> unlist()
crr_var <- strsplit(x = crr_dte, split = '_')[[1]][1]
crr_dte <- gsub(paste0(crr_var,'_'), '', crr_dte)
crr_dte <- gsub('.tif', '', crr_dte)

dte <- as.Date(crr_dte)
raw_pth <- paste0('~/common_data/nex-gddp-cmip6_raw/hurs/',crr_ssp,'/',crr_gcm)
raw_fls <- list.files(path = raw_pth, pattern = '.nc$', full.names = T, recursive = F)
raw_fle <- raw_fls[grep(lubridate::year(dte), raw_fls)]

raw_r <- terra::rast(raw_fle)
# Precipitation: pr * 86400
# Temperatures: tasmax, tasmin - 273.15
# Solar radiation: rsds * 86400 / 1000000
if (crr_var == 'pr') { raw_r <- raw_r * 86400 }
if (crr_var == 'rsds') { raw_r <- raw_r * 86400 / 1000000 }
if (crr_var %in% c('tasmax','tasmin')) { raw_r <- raw_r - 273.15 }
ppr_r <- raw_r[[which(terra::time(raw_r) == dte)]] |> terra::rotate()
terra::writeRaster(x = ppr_r, filename = fle, overwrite = T)