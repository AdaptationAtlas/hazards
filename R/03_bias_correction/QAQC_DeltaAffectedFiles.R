# Quality control for delta, monthly, and daily precipitation files derived from delta method
# H. Achicanoy
# Alliance of Bioversity-CIAT
# July, 2025

# R options
options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
pacman::p_load(terra, tidyverse, FactoMineR)
# Set working directory
setwd('~')

## Delta's quality control ----
# Get delta statistics
dlt_pth <- '~/common_data/esfg_cmip6/intermediate/interpolated_mthly_anomaly'
system(paste0('~/common_data/cloud-convert run-qaqc ',dlt_pth,' --quantiles --pct-check 100'))
# Read delta statistics
delta_sts <- utils::read.csv(file.path(dlt_pth,'qaqc.csv'))
delta_sts <- delta_sts[grep('_fixed',delta_sts$file),]
delta_sts <- delta_sts[-grep('_thr',delta_sts$file),]
# Check minimum and maximum values
delta_sts |>
  ggplot2::ggplot(aes(x = min, y = max)) +
  ggplot2::geom_point(alpha = 0.1) +
  ggplot2::theme_bw()
# PCA
delta.pca.res <- delta_sts |> dplyr::select(mean, min, max, stdev, q1, median, q3) |> FactoMineR::PCA(scale.unit = T, ncp = 7, graph = T)

## Monthly's quality control ----
# Setup table
ssps <- c('ssp126','ssp245','ssp370','ssp585')
gcms <- c('ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0')
prds <- c('2021_2040','2041_2060','2061_2080','2081_2100')
stp_tbl <- expand.grid(ssp = ssps, gcm = gcms, prd = prds, stringsAsFactors = F) |>
  base::as.data.frame() |>
  dplyr::arrange(gcm, ssp, prd)
rm(ssps, gcms, prds)
stp_tbl$folder <- paste0(stp_tbl$ssp,'_',stp_tbl$gcm,'_',stp_tbl$prd)
# Get monthly statistics
root <- '~/common_data'
app_call <- '~/common_data/cloud-convert run-qaqc '
app_args <- ' --quantiles --pct-check 100'
for (i in 1:nrow(stp_tbl)) {
  
  cat('Processing', stp_tbl$folder[i],'\n')
  system(command =
           paste0(app_call,
                  file.path(root,'atlas_hazards/cmip6/indices',stp_tbl$folder[i],'PTOT'),
                  app_args)
  )
  
}; rm(i)
# Read monthly statistics
mthly_sts <- 1:nrow(stp_tbl) |>
  purrr::map(.f = function(i) {
    
    sts <- utils::read.csv(file.path(root,'atlas_hazards/cmip6/indices',stp_tbl$folder[i],'PTOT/qaqc.csv'))
    return(sts)
    
  }) |> dplyr::bind_rows()
mthly_sts <- mthly_sts[grep('PTOT-[0-9][0-9][0-9][0-9]-[0-9][0-9].tif', mthly_sts$file),]
# Check median and maximum values
mthly_sts |>
  ggplot2::ggplot(aes(x = median, y = max)) +
  ggplot2::geom_point(alpha = 0.1) +
  ggplot2::theme_bw()
# PCA
mthly.pca.res <- mthly_sts |> dplyr::select(mean, min, max, stdev, q1, median, q3) |> FactoMineR::PCA(scale.unit = T, ncp = 7, graph = T)

## Daily's quality control ----
# Get daily statistics
stp_tbl$daily_folder <- paste0('Prec_',stp_tbl$gcm,'_',stp_tbl$ssp,'_',stp_tbl$prd)
for (i in 1:nrow(stp_tbl)) {
  
  cat('Processing', stp_tbl$daily_folder[i],'\n')
  system(command =
           paste0(app_call,
                  file.path(root,'chirps_cmip6_africa',stp_tbl$daily_folder[i]),
                  app_args)
  )
  
}; rm(i)
# Read daily statistics
daily_sts <- 1:nrow(stp_tbl) |>
  purrr::map(.f = function(i) {
    
    sts <- utils::read.csv(file.path(root,'chirps_cmip6_africa',stp_tbl$daily_folder[i],'qaqc.csv'))
    return(sts)
    
  }) |> dplyr::bind_rows()
daily_sts <- daily_sts[grep('chirps-v2.0*.*.tif', daily_sts$file),]
# Check Q3 and maximum values
daily_sts |>
  ggplot2::ggplot(aes(x = q3, y = max)) +
  ggplot2::geom_point(alpha = 0.1) +
  ggplot2::theme_bw()
# PCA
daily.pca.res <- daily_sts |> dplyr::select(mean, min, max, stdev, q1, median, q3) |> FactoMineR::PCA(scale.unit = T, ncp = 7, graph = T)

# library(terra)
# r <- terra::rast(file.path(root,'atlas_hazards/cmip6/indices',stp_tbl$folder[i],'PTOT/PTOT-2021-01.tif'))
# terra::extract(x = r, y = data.frame(lon = -6.40778, lat = 9.29731))