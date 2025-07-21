options(warn = -1, scipen = 999)
library(pacman)
pacman::p_load(terra, geodata, tidyverse, future, furrr)
# Own functions'
list.files2 <- Vectorize(list.files, vectorize.args = 'pattern')

root <- '/home/jovyan/common_data'

## Mask creation ----
# Region of interest
country <- 'BFA' # KEN, BFA, GHA, AGO
bbox <- geodata::gadm(country, 0, tempdir())
# Load production data in tons from MapSPAM
mapspam_all <- terra::rast(file.path(root,'hazards_prototype/Data/mapspam/2020V1r2_SSA/processed/variable=prod_t/spam_prod_t_all.tif'))
mapspam_maize <- mapspam_all[['maize']]
mapspam_maize <- mapspam_maize |> terra::crop(terra::ext(bbox)) |> terra::mask(bbox)
mapspam_maize <- terra::classify(mapspam_maize, cbind(-Inf, 100, NA))
# Template raster from Nex-GDDP-CMIP6
tmp25 <- terra::rast(file.path(root,'atlas_nex-gddp_hazards/cmip6/indices/ssp245_ACCESS-ESM1-5_2021_2040/NDWS/NDWS-2021-01.tif'))
tmp25 <- tmp25 |> terra::crop(terra::ext(bbox)) |> terra::mask(bbox)
# Resampling maize area at 0.25 degrees
msk <- terra::resample(x = mapspam_maize, y = tmp25)
terra::writeRaster(msk, '~/maize_mask.tif', overwrite = T)

## Setup table ----
index <- 'NDWS'
if (index == 'TAI') {
  ssps <- c('ssp126','ssp245','ssp370','ssp585')
  gcms <- c('ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0')
  prds <- c('2021_2040','2041_2060','2061_2080','2081_2100')
  stp_tbl <- expand.grid(ssp = ssps, gcm = gcms, prd = prds, stringsAsFactors = F) |>
    base::as.data.frame() |>
    dplyr::arrange(gcm, ssp, prd)
  rm(ssps, gcms, prds)
} else {
  ssps <- c('ssp126','ssp245','ssp370','ssp585')
  gcms <- c('ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0')
  prds <- c('2021_2040','2041_2060','2061_2080','2081_2100')
  mns  <- sprintf('%02.0f',1:12)
  stp_tbl <- expand.grid(ssp = ssps, gcm = gcms, prd = prds, mn = mns, stringsAsFactors = F) |>
    base::as.data.frame() |>
    dplyr::arrange(gcm, ssp, prd)
  rm(ssps, gcms, prds, mns)
}

future::plan(future::multisession, workers = 15)
indices_qaqc <- furrr::future_map(.x = 1:nrow(stp_tbl), .f = function(i) {
  
  ssp <- stp_tbl$ssp[i]
  gcm <- stp_tbl$gcm[i]
  prd <- stp_tbl$prd[i]
  yini <- strsplit(prd, '_')[[1]][1] |> as.numeric()
  yend <- strsplit(prd, '_')[[1]][2] |> as.numeric()
  if (index != 'TAI') { mn <- stp_tbl$mn[i] }
  
  # Wrap and unwrapped mask
  msk <- terra::rast('~/maize_mask.tif')
  msk_wrp <- terra::wrap(msk)
  msk_uwr <- terra::unwrap(msk_wrp); rm(msk_wrp)
  bbox <- geodata::gadm(country, 0, tempdir())
  bbox_wrp <- terra::wrap(bbox)
  bbox_uwr <- terra::unwrap(bbox_wrp); rm(bbox_wrp)
  
  ## Atlas delta indices ----
  if (index == 'TAI') {
    fls_dlt <- paste0(root,'/atlas_hazards/cmip6/indices/',ssp,'_',gcm,'_',prd,'/',index,'/',index,'-',yini:yend,'.tif')
  } else {
    fls_dlt <- paste0(root,'/atlas_hazards/cmip6/indices/',ssp,'_',gcm,'_',prd,'/',index,'/',index,'-',yini:yend,'-',mn,'.tif')
  }
  idx_dlt <- terra::rast(fls_dlt)
  idx_dlt_wrp <- terra::wrap(idx_dlt)
  idx_dlt_uwr <- terra::unwrap(idx_dlt_wrp); rm(idx_dlt_wrp)
  if (index %in% c('NDD','NDWS','NDWL0','NDWL50')) {mthd = 'near'} else {mthd = 'bilinear'}
  idx_dlt_uwr <- idx_dlt_uwr |>
    terra::crop(terra::ext(bbox_uwr)) |>
    terra::resample(msk_uwr, method = mthd) |>
    terra::mask(msk_uwr)
  mdn_dlt <- terra::values(idx_dlt_uwr, na.rm = T) |> apply(2, median) |> median()
  
  ## Atlas nexgddp indices ----
  if (index == 'TAI') {
    fls_nxg <- paste0(root,'/atlas_nex-gddp_hazards/cmip6/indices/',ssp,'_',gcm,'_',prd,'/',index,'/',index,'-',yini:yend,'.tif')
  } else {
    fls_nxg <- paste0(root,'/atlas_nex-gddp_hazards/cmip6/indices/',ssp,'_',gcm,'_',prd,'/',index,'/',index,'-',yini:yend,'-',mn,'.tif')
  }
  idx_nxg <- terra::rast(fls_nxg)
  idx_nxg_wrp <- terra::wrap(idx_nxg)
  idx_nxg_uwr <- terra::unwrap(idx_nxg_wrp); rm(idx_nxg_wrp)
  idx_nxg_uwr <- idx_nxg_uwr |>
    terra::crop(terra::ext(bbox_uwr)) |>
    terra::mask(msk_uwr)
  mdn_nxg <- terra::values(idx_nxg_uwr, na.rm = T) |> apply(2, median) |> median()
  
  res <- stp_tbl[i,]
  res$delta <- mdn_dlt
  res$nexgddp <- mdn_nxg
  return(res)
}, .progress = T) |> dplyr::bind_rows()
future::plan(future::sequential)

indices_qaqc_lng <- indices_qaqc |> tidyr::pivot_longer(cols = delta:nexgddp, names_to = 'dataset', values_to = 'value')

if (index == 'TAI') {
  indices_qaqc_lng |>
    ggplot2::ggplot(aes(x = prd, y = value, colour = gcm, group = gcm)) +
    ggplot2::geom_line() +
    ggplot2::xlab('Month') +
    ggplot2::ylab(index) +
    ggplot2::facet_grid(ssp ~ dataset) +
    ggplot2::theme(legend.position = 'bottom') +
    ggplot2::theme_bw()
} else {
  indices_qaqc_lng |>
    dplyr::filter(prd == '2081_2100') |>
    dplyr::mutate(mn = factor(mn, levels = sprintf('%02d', 1:12), ordered = T)) |>
    ggplot2::ggplot(aes(x = mn, y = value, colour = gcm, group = gcm)) +
    ggplot2::geom_line() +
    ggplot2::xlab('Month') +
    ggplot2::ylab(index) +
    ggplot2::facet_grid(ssp ~ dataset) +
    ggplot2::theme(legend.position = 'bottom') +
    ggplot2::theme_bw()
}

## Posibles diferencias ----
# Different baseline: Atlas-delta uses CHIRPS/CHIRTS, Nex-GDDP-CMIP6 uses 
# Different bias-correction method
# Different spatial resolutions
# Indices affected by bad daily files

indices_qaqc <- 1:nrow(stp_tbl) |>
  purrr::map(.f = function(i) {
    
    ssp <- stp_tbl$ssp[i]
    gcm <- stp_tbl$gcm[i]
    prd <- stp_tbl$prd[i]
    yini <- strsplit(prd, '_')[[1]][1] |> as.numeric()
    yend <- strsplit(prd, '_')[[1]][2] |> as.numeric()
    mn <- stp_tbl$mn[i]
    
    ## Atlas delta indices ----
    fls_dlt <- paste0(root,'/atlas_hazards/cmip6/indices/',ssp,'_',gcm,'_',prd,'/',index,'/',index,'-',yini:yend,'-',mn,'.tif')
    idx_dlt <- terra::rast(fls_dlt)
    if (index %in% c('NDD','NDWS','NDWL0','NDWL50')) {mthd = 'near'} else {mthd = 'bilinear'}
    idx_dlt <- idx_dlt |>
      terra::crop(terra::ext(bbox)) |>
      terra::resample(tmp25, method = mthd) |>
      terra::mask(msk)
    mdn_dlt <- terra::values(idx_dlt, na.rm = T) |> apply(2, median) |> median()
    
    ## Atlas nexgddp indices ----
    fls_nxg <- paste0(root,'/atlas_nex-gddp_hazards/cmip6/indices/',ssp,'_',gcm,'_',prd,'/',index,'/',index,'-',yini:yend,'-',mn,'.tif')
    idx_nxg <- terra::rast(fls_nxg)
    idx_nxg <- idx_nxg |>
      terra::crop(terra::ext(bbox)) |>
      terra::mask(msk)
    mdn_nxg <- terra::values(idx_nxg, na.rm = T) |> apply(2, median) |> median()
    
    res <- stp_tbl[i,]
    res$delta <- mdn_dlt
    res$nexgddp <- mdn_nxg
    return(res)
  }) |> dplyr::bind_rows()
