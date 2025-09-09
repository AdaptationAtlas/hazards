# Quality control for daily files derived from Nex-GDDP-CMIP6
# H. Achicanoy
# Alliance of Bioversity-CIAT
# September, 2025

# R options
options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
pacman::p_load(terra, tidyverse, FactoMineR, arrow)
list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'pattern')

# Setup table
scenario <- 'historical' # historical, future
gcms <- c('ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CMCC-ESM2','EC-Earth3','EC-Earth3-Veg-LR','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','TaiESM1')
if (scenario == 'future') {
  ssps <- c('ssp126','ssp245','ssp370','ssp585')
  prds <- c('2021_2100')
  stp_tbl <- expand.grid(ssp = ssps, gcm = gcms, prd = prds, stringsAsFactors = F) |> base::as.data.frame() |> dplyr::arrange(gcm, ssp, prd)
  stp_tbl$folder <- paste0(stp_tbl$ssp,'_',stp_tbl$gcm,'_',stp_tbl$prd)
  stp_tbl$ini_year <- strsplit(stp_tbl$prd, '_') |> purrr::map(1) |> unlist() |> as.numeric()
  stp_tbl$end_year <- strsplit(stp_tbl$prd, '_') |> purrr::map(2) |> unlist() |> as.numeric()
} else {
  if (scenario == 'historical') {
    ssps <- 'historical'
    stp_tbl <- expand.grid(ssp = ssps, gcm = gcms, stringsAsFactors = F) |> base::as.data.frame() |> dplyr::arrange(gcm, ssp)
    stp_tbl$folder <- paste0(stp_tbl$ssp,'_',stp_tbl$gcm)
    stp_tbl$ini_year <- 1981 # 1995
    stp_tbl$end_year <- 2014
  }
}

root <- '/home/jovyan/common_data'

vrs <- c('hurs','pr','rsds','tasmax','tasmin')
# vrs <- c('pr','hurs')
# vrs <- paste0(vrs,'2')

# Applying QAQC tool to all folders
dfm_qaqc <- purrr::map(.x = vrs, .f = function(index) {
  
  cat('.... Validating files for index:',index,'\n')
  
  dfm_qaqc <- purrr::map(.x = 1:nrow(stp_tbl), .f = function(i) {
    
    cat('Over the scenario:',stp_tbl$folder[i],'\n')
    # Target
    trg_pth <- paste0(root,'/nex-gddp-cmip6/',index,'/',stp_tbl$ssp[i],'/',stp_tbl$gcm[i])
    
    app_call <- '~/common_data/cloud-convert run-qaqc '
    app_args <- ' --quantiles --pct-check 100'
    
    if (!file.exists(file.path(trg_pth,'qaqc.csv'))) {
      cat('Processing', trg_pth,'\n')
      system(command = paste0(app_call, trg_pth, app_args))
    }
    
    dfm_qaqc <- utils::read.csv(file.path(trg_pth,'qaqc.csv'))
    dfm_qaqc$index <- index
    dfm_qaqc$folder <- trg_pth
    dfm_qaqc$name <- as.character(dfm_qaqc$name)
    return(dfm_qaqc)
    
  }) |> dplyr::bind_rows()
  
  return(dfm_qaqc)
  
}) |> dplyr::bind_rows()

pca.res <- dfm_qaqc |>
  dplyr::filter(index == 'PTOT') |>
  dplyr::select(mean, min, max, stdev, q1, median, q3) |> FactoMineR::PCA(scale.unit = T, ncp = 7, graph = T)

ptot_issue <- dfm_qaqc |>
  dplyr::filter(index == 'PTOT')

ptot_issue$folder[which.max(ptot_issue$max)]

# Counting of number of files
summary(
  table(
    dfm_qaqc$folder[!(dfm_qaqc$index %in% c('TAI','HSH','THI'))]
  ) |> as.numeric()
)

summary(
  table(
    dfm_qaqc$folder[(dfm_qaqc$index %in% c('TAI'))]
  ) |> as.numeric()
)

summary(
  table(
    dfm_qaqc$folder[(dfm_qaqc$index %in% c('THI'))]
  ) |> as.numeric()
)

# Verifying if there are corrupted files
table(dfm_qaqc$folder[is.na(dfm_qaqc$mean)]) |> sort(decreasing = T)


# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## Indices quality control ----
# # Get indices statistics
# dlt_sts_file <- '~/common_data/affected_geographies/corrected_delta_stats.parquet'
# if (!file.exists(dlt_sts_file)) {
#   dlt_pth <- '~/common_data/esfg_cmip6/intermediate/interpolated_mthly_anomaly'
#   system(paste0('~/common_data/cloud-convert run-qaqc ',dlt_pth,' --quantiles --pct-check 100'))
#   # Read delta statistics
#   delta_sts <- utils::read.csv(file.path(dlt_pth,'qaqc.csv'))
#   delta_sts <- delta_sts[grep('_fixed',delta_sts$file),]
#   delta_sts <- delta_sts[-grep('_thr',delta_sts$file),]
#   arrow::write_parquet(delta_sts, dlt_sts_file)
# } else {
#   delta_sts <- arrow::read_parquet(dlt_sts_file) |> base::as.data.frame()
# }; rm(dlt_sts_file)
# # Check minimum and maximum values
# delta_sts |>
#   ggplot2::ggplot(aes(x = min, y = max)) +
#   ggplot2::geom_point(alpha = 0.1) +
#   ggplot2::theme_bw()
# # PCA
# delta.pca.res <- delta_sts |> dplyr::select(mean, min, max, stdev, q1, median, q3) |> FactoMineR::PCA(scale.unit = T, ncp = 7, graph = T)
# 
# ## Monthly's quality control ----
# # Setup table
# ssps <- c('ssp126','ssp245','ssp370','ssp585')
# gcms <- c('ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0')
# prds <- c('2021_2040','2041_2060','2061_2080','2081_2100')
# stp_tbl <- expand.grid(ssp = ssps, gcm = gcms, prd = prds, stringsAsFactors = F) |>
#   base::as.data.frame() |>
#   dplyr::arrange(gcm, ssp, prd)
# rm(ssps, gcms, prds)
# stp_tbl$folder <- paste0(stp_tbl$ssp,'_',stp_tbl$gcm,'_',stp_tbl$prd)
# # Get monthly statistics
# mnt_sts_file <- '~/common_data/affected_geographies/corrected_monthly_stats.parquet'
# if (!file.exists(mnt_sts_file)) {
#   root <- '~/common_data'
#   app_call <- '~/common_data/cloud-convert run-qaqc '
#   app_args <- ' --quantiles --pct-check 100'
#   for (i in 1:nrow(stp_tbl)) {
#     
#     cat('Processing', stp_tbl$folder[i],'\n')
#     system(command =
#              paste0(app_call,
#                     file.path(root,'atlas_hazards/cmip6/indices',stp_tbl$folder[i],'PTOT'),
#                     app_args)
#     )
#     
#   }; rm(i)
#   # Read monthly statistics
#   mthly_sts <- 1:nrow(stp_tbl) |>
#     purrr::map(.f = function(i) {
#       
#       sts <- utils::read.csv(file.path(root,'atlas_hazards/cmip6/indices',stp_tbl$folder[i],'PTOT/qaqc.csv'))
#       return(sts)
#       
#     }) |> dplyr::bind_rows()
#   mthly_sts <- mthly_sts[grep('PTOT-[0-9][0-9][0-9][0-9]-[0-9][0-9].tif', mthly_sts$file),]
#   arrow::write_parquet(mthly_sts, mnt_sts_file)
# } else {
#   mthly_sts <- arrow::read_parquet(mnt_sts_file) |> base::as.data.frame()
# }; rm(mnt_sts_file)
# # Check median and maximum values
# mthly_sts |>
#   ggplot2::ggplot(aes(x = median, y = max)) +
#   ggplot2::geom_point(alpha = 0.1) +
#   ggplot2::theme_bw()
# # PCA
# mthly.pca.res <- mthly_sts |> dplyr::select(mean, min, max, stdev, q1, median, q3) |> FactoMineR::PCA(scale.unit = T, ncp = 7, graph = T)
# 
# ## Daily's quality control ----
# # Get daily statistics
# dly_sts_file <- '~/common_data/affected_geographies/corrected_daily_stats.parquet'
# if (!file.exists(dly_sts_file)) {
#   stp_tbl$daily_folder <- paste0('Prec_',stp_tbl$gcm,'_',stp_tbl$ssp,'_',stp_tbl$prd)
#   for (i in 1:nrow(stp_tbl)) {
#     
#     cat('Processing', stp_tbl$daily_folder[i],'\n')
#     system(command =
#              paste0(app_call,
#                     file.path(root,'chirps_cmip6_africa',stp_tbl$daily_folder[i]),
#                     app_args)
#     )
#     
#   }; rm(i)
#   # Read daily statistics
#   daily_sts <- 1:nrow(stp_tbl) |>
#     purrr::map(.f = function(i) {
#       
#       sts <- utils::read.csv(file.path(root,'chirps_cmip6_africa',stp_tbl$daily_folder[i],'qaqc.csv'))
#       return(sts)
#       
#     }) |> dplyr::bind_rows()
#   daily_sts <- daily_sts[grep('chirps-v2.0*.*.tif', daily_sts$file),]
#   arrow::write_parquet(daily_sts, dly_sts_file)
# } else {
#   daily_sts <- arrow::read_parquet(dly_sts_file) |> base::as.data.frame()
# }; rm(dly_sts_file)
# # Check Q3 and maximum values
# daily_sts |>
#   ggplot2::ggplot(aes(x = q3, y = max)) +
#   ggplot2::geom_point(alpha = 0.1) +
#   ggplot2::theme_bw()
# # PCA
# daily.pca.res <- daily_sts |> dplyr::select(mean, min, max, stdev, q1, median, q3) |> FactoMineR::PCA(scale.unit = T, ncp = 7, graph = T)
# 
# # library(terra)
# # r <- terra::rast(file.path(root,'atlas_hazards/cmip6/indices',stp_tbl$folder[i],'PTOT/PTOT-2021-01.tif'))
# # terra::extract(x = r, y = data.frame(lon = -6.40778, lat = 9.29731))
# 
# # Solar radiation
# if (!file.exists(dly_sts_file)) {
#   stp_tbl$daily_folder <- paste0('solar_radiation_flux_',stp_tbl$gcm,'_',stp_tbl$ssp,'_',stp_tbl$prd)
#   for (i in 1:nrow(stp_tbl)) {
#     
#     cat('Processing', stp_tbl$daily_folder[i],'\n')
#     system(command =
#              paste0(app_call,
#                     file.path(root,'ecmwf_agera5_cmip6_africa',stp_tbl$daily_folder[i]),
#                     app_args)
#     )
#     
#   }; rm(i)
#   # Read daily statistics
#   daily_sts <- 1:nrow(stp_tbl) |>
#     purrr::map(.f = function(i) {
#       
#       sts <- utils::read.csv(file.path(root,'ecmwf_agera5_cmip6_africa',stp_tbl$daily_folder[i],'qaqc.csv'))
#       return(sts)
#       
#     }) |> dplyr::bind_rows()
#   #daily_sts <- daily_sts[grep('chirps-v2.0*.*.tif', daily_sts$file),]
#   arrow::write_parquet(daily_sts, dly_sts_file)
# } else {
#   daily_sts <- arrow::read_parquet(dly_sts_file) |> base::as.data.frame()
# }; rm(dly_sts_file)
# 
# daily_sts |>
#   ggplot2::ggplot(aes(x = q3, y = max)) +
#   ggplot2::geom_point(alpha = 0.1) +
#   ggplot2::theme_bw()
# # PCA
# daily.pca.res <- daily_sts |> dplyr::select(mean, min, max, stdev, q1, median, q3) |> FactoMineR::PCA(scale.unit = T, ncp = 7, graph = T)
