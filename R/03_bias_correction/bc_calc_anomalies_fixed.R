#calculate and interpolate anomalies from GCM data
#JRV, Dec 2022
#HAE, Jun 2025

#load libraries
if(!require(pacman)){install.packages('pacman');library(pacman)} else library(pacman)
pacman::p_load(tidyverse, terra, fields, furrr, future)

#options
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation

#working directory
wd <- "~/common_data/esfg_cmip6"
raw_dir <- paste0(wd, "/raw")
rsds_dir <- paste0(wd, "/raw_rsds")
hurs_dir <- paste0(wd, "/raw_hurs")
int_dir <- paste0(wd, "/intermediate")
if (!file.exists(int_dir)) {dir.create(int_dir)}

mth_dir <- paste0(int_dir, "/monthly_annual")
if (!file.exists(mth_dir)) {dir.create(mth_dir)}

clm_dir <- paste0(int_dir, "/monthly_climatology")
if (!file.exists(clm_dir)) {dir.create(clm_dir)}

anom_dir <- paste0(int_dir, "/interpolated_mthly_anomaly")
if (!file.exists(anom_dir)) {dir.create(anom_dir)}

#list of GCMs
gcm_list <- c("CMIP6_ACCESS-ESM1-5",
              "CMIP6_MPI-ESM1-2-HR",
              "CMIP6_EC-Earth3",
              "CMIP6_INM-CM5-0",
              "CMIP6_MRI-ESM2-0")

####
#function to calculate climatology
calc_climatology <- function(data_file, period, sce_lab, gcm_name, varname, mth_dir, clm_dir) {
  #load data
  #data_file <- rcp_file
  #period <- rcp_period
  #sce_lab <- rcp #"historical"
  
  #load raster
  r_data <- terra::rast(data_file)
  if (varname %in% c("rsds","hurs") & sce_lab != "historical") {names(r_data) <- paste(terra::time(r_data))}
  
  #date data.frame
  date_df <- names(r_data) %>%
    purrr::map_chr(.f=function(.x) {gsub(pattern=" 12:00:00 GMT", replacement="", x=.x)})
  date_df <- data.frame(fullname=names(r_data), fulldate=date_df)
  date_df <- date_df %>%
    dplyr::mutate(year=purrr::map_chr(.x=1:nrow(date_df), .f=function(.x) {substr(.$fulldate[.x], 1, 4)})) %>%
    dplyr::mutate(month=purrr::map_chr(.x=1:nrow(date_df), .f=function(.x) {substr(.$fulldate[.x], 6, 7)})) %>%
    dplyr::mutate(day=purrr::map_chr(.x=1:nrow(date_df), .f=function(.x) {substr(.$fulldate[.x], 9, 10)})) %>%
    dplyr::mutate(year=as.numeric(year), month=as.numeric(month), day=as.numeric(day)) %>%
    dplyr::mutate(index=1:nrow(.)) %>%
    dplyr::filter(year %in% period)
  
  #monthly file data.frame
  month_df <- date_df %>%
    dplyr::select(year, month) %>%
    dplyr::distinct(.) %>%
    dplyr::arrange(year, month) %>%
    dplyr::mutate(index=1:nrow(.))
  
  mth_fname <- paste0(mth_dir, "/", gcm_name, "_", sce_lab, "_r1i1p1f1_", varname, "_Africa_monthly_annual_", min(period), "_", max(period), ".tif")
  if (!file.exists(mth_fname)) {
    r_monthly <- c()
    for (i in 1:nrow(month_df)) {
      #month and year
      #i <- 1
      mth <- month_df$month[i]
      yr <- month_df$year[i]
      #cat("processing i=", i, "/ year=", yr, "/ month=", mth, "\n")
      tdates <- date_df %>%
        dplyr::filter(year==yr, month==mth)
      
      #subset raster based on relevant indices
      r_mth <- r_data[[tdates$index]]
      if (varname == "pr") {r_mth <- sum(r_mth, na.rm=TRUE)} else {r_mth <- mean(r_mth, na.rm=TRUE)}
      names(r_mth) <- paste0(yr, "-", sprintf("%02.0f",mth))
      
      #append
      r_monthly <- c(r_monthly, r_mth)
    }
    #create raster, write
    r_monthly <- terra::rast(r_monthly)
    terra::writeRaster(r_monthly, mth_fname)
  } else {
    r_monthly <- terra::rast(mth_fname)
  }
  
  #now calculate climatological means
  clm_fname <- paste0(clm_dir, "/", gcm_name, "_", sce_lab, "_r1i1p1f1_", varname, "_Africa_monthly_climatology_", min(period), "_", max(period), ".tif")
  if (!file.exists(clm_fname)) {
    r_clm <- c()
    for (i in 1:12) {
      #i <- 1
      #cat("processing month i=", i, "\n")
      tdates <- month_df %>%
        dplyr::filter(month==i)
      
      #subset raster based on relevant indices
      r_yrs <- r_monthly[[tdates$index]]
      r_yrs <- mean(r_yrs, na.rm=TRUE)
      names(r_yrs) <- paste0(sprintf("%02.0f",i))
      
      #append
      r_clm <- c(r_clm, r_yrs)
    }
    #create raster, write
    r_clm <- terra::rast(r_clm)
    terra::writeRaster(r_clm, clm_fname)
  } else {
    r_clm <- terra::rast(clm_fname)
  }
  
  #return object
  return(list(monthly=r_monthly, climatology=r_clm))
}
####

####
#function to interpolate monthly anomalies
intp_anomalies <- function(his_clm, rcp_clm, anom_dir, ref, gcm_name, rcp, varname, period, month) {
  #final output file name
  anom_fname <- paste0(anom_dir, '/', gcm_name, '_', rcp, '_r1i1p1f1_', varname, '_Africa_monthly_intp_anomaly_', min(period), '_', max(period), '_fixed.tif')
  
  #create temporary directory to store monthly output (to reduce impact of server crashes)
  anom_dname <- paste0(anom_dir, '/tmp_', gcm_name, '_', rcp, '_r1i1p1f1_', varname, '_Africa_monthly_intp_anomaly_', min(period), '_', max(period))
  if (!file.exists(anom_dname)) {dir.create(anom_dname)}
  
  if (!file.exists(anom_fname)) {
    
    rcp_clm <- terra::wrap(rcp_clm)
    his_clm <- terra::wrap(his_clm)
    ref <- terra::wrap(ref)
    
    plan(multisession, workers = 12)
    r_anom <- furrr::future_map(.x = 1:12, .f = function(.x){
      
      cat('processing interpolation for month i=', .x, '\n')
      
      rcp_clm_u <- terra::unwrap(rcp_clm)
      his_clm_u <- terra::unwrap(his_clm)
      ref_u <- terra::unwrap(ref); names(ref_u) <- 'mean'
      
      #get climatology rasters
      avg_fut <- rcp_clm_u[[.x]]
      avg_his <- his_clm_u[[.x]]
      
      #calculate anomaly
      cat('calculating anomaly...\n')
      if (varname %in% c('tasmax','tasmin','tas','hurs')) {
        anom <- avg_fut - avg_his
      } else {
        # if precipitation is below or very close to zero make it 0.01
        avg_his <- terra::classify(avg_his, cbind(-Inf, 30, 30)) # if historical has exactly zero make it non-zero, even values close to zero ### avg_his[avg_his[] == 0] <- 0.01
        avg_fut <- terra::classify(avg_fut, cbind(-Inf, 30, 30)) # if future has exactly zero make it non-zero, even values close to zero
        
        #calculate anomaly in fraction
        anom <- (avg_fut - avg_his)/avg_his
        
        # Truncate the top 2% of anomaly values
        thr <- as.numeric(terra::global(x = anom, fun = stats::quantile, probs = 0.98, na.rm = T))
        anom[anom >= thr] <- thr
      }
      names(anom) <- 'mean'
      
      #as data.frame
      if (gcm_name %in% c('CMIP6_EC-Earth3','CMIP6_MPI-ESM1-2-HR')) {
        crds <- terra::as.data.frame(anom, xy = T)
        crds <- terra::spatSample(x = anom, size = ceiling(nrow(crds) * 0.3), method = 'regular', xy = T)
      } else {
        crds <- terra::as.data.frame(anom, xy = T)
      }
      
      #fit tps interpolation model
      cat('fitting thin plate spline\n')
      if (!file.exists(paste0(anom_dname,'/tps_mth_',.x,'_fixed.RData'))) {
        library(fields)
        tps <- fields::Tps(x = crds[,c('x','y')], Y = crds[,'mean'])
        save(tps, file=paste0(anom_dname,'/tps_mth_',.x,'_fixed.RData'))
      } else {
        load(paste0(anom_dname,'/tps_mth_',.x,'_fixed.RData'))
      }
      
      #interpolate
      cat('interpolating onto the reference raster\n')
      if (!file.exists(paste0(anom_dname,'/raster_mth_',.x,'_fixed.tif'))) {
        # intp <- terra::interpolate(object = ref_u, model = tps)
        xy <- terra::as.data.frame(ref_u, xy = T)[,1:2]
        intp_vls <- fields::predict.Krig(tps, xy)
        xyz <- cbind(xy, base::as.data.frame(intp_vls))
        intp <- terra::rast(xyz)
        terra::crs(intp) <- terra::crs(ref_u)
        names(intp) <- paste0(sprintf('%02.0f',.x))
        terra::writeRaster(intp, paste0(anom_dname,'/raster_mth_',.x,'_fixed.tif'), overwrite = T)
      } else {
        intp <- terra::rast(paste0(anom_dname,'/raster_mth_',.x,'_fixed.tif'))
      }
      
      #clean-up
      rm(list=c('tps','anom','avg_fut','avg_his','crds'))
      gc(verbose = F, full = T, reset = T)
      
      #append
      r_anom <- terra::wrap(intp)
      return(r_anom)
      
    }, .progress = T)
    plan(sequential)
    r_anom <- r_anom |> purrr::map(terra::unwrap) |> terra::rast()
    
    #create raster, write
    terra::writeRaster(r_anom, anom_fname)
    
    #final cleanup
    rm(intp)
    gc(verbose = F, full = T, reset = T)
  } else {
    r_anom <- terra::rast(anom_fname)
  }
  
  #delete temporary directory
  if (file.exists(anom_dname)) {system(paste0('rm -rf ',anom_dname))}
  
  #return object
  return(r_anom)
}
####

####
#parameters
####
gcm_i <- 3
varname <- 'pr'
rcp <- ssp
futperiod <- dplyr::case_when(
  prd == '2021_2040' ~ 'near',
  prd == '2041_2060' ~ 'mid',
  prd == '2061_2080' ~ 'far',
  prd == '2081_2100' ~ 'end',
  .default = as.character(prd)
)

#define historical and future periods
his_period <- 1995:2014
if (futperiod == "near") {rcp_period <- 2021:2040}
if (futperiod == "mid") {rcp_period <- 2041:2060}
if (futperiod == "far") {rcp_period <- 2061:2080}
if (futperiod == "end") {rcp_period <- 2081:2100}

#data files
if (varname == "rsds") {
  his_file <- paste0(rsds_dir, "/", gcm_list[gcm_i], "_historical_", varname, "_Africa_daily.tif")
  rcp_file <- paste0(rsds_dir, "/", gcm_list[gcm_i], "_", rcp, "_", varname, "_Africa_daily.tif")
} else if (varname == "hurs") {
  his_file <- paste0(hurs_dir, "/", gcm_list[gcm_i], "_historical_", varname, "_Africa_daily.tif")
  rcp_file <- paste0(hurs_dir, "/", gcm_list[gcm_i], "_", rcp, "_", varname, "_Africa_daily.tif")
} else {
  his_file <- paste0(raw_dir, "/", gcm_list[gcm_i], "_historical_r1i1p1f1_", varname, "_Africa_daily.tif")
  if (futperiod %in% c("near", "mid")) {
    rcp_file <- paste0(raw_dir, "/", gcm_list[gcm_i], "_", rcp, "_r1i1p1f1_", varname, "_Africa_daily.tif")
  } else {
    rcp_file <- paste0(raw_dir, "/", gcm_list[gcm_i], "_", rcp, "_r1i1p1f1_", varname, "_Africa_daily_2061-2100.tif")
  }
}

#historical climatology
his_clm <- calc_climatology(data_file = his_file, 
                            period    = his_period, 
                            sce_lab   = 'historical',
                            gcm_name  = gcm_list[gcm_i],
                            varname   = varname,
                            mth_dir   = mth_dir,
                            clm_dir   = clm_dir)$climatology
his_clm <- terra::classify(his_clm, cbind(-Inf, 0, 0))

#future climatology
rcp_clm <- calc_climatology(data_file = rcp_file, 
                            period    = rcp_period, 
                            sce_lab   = rcp,
                            gcm_name  = gcm_list[gcm_i],
                            varname   = varname,
                            mth_dir   = mth_dir,
                            clm_dir   = clm_dir)$climatology
rcp_clm <- terra::classify(rcp_clm, cbind(-Inf, 0, 0))

#reference CHIRPS/CHIRTS raster
r_ref <- terra::rast('~/common_data/chirts/Tmax/1995/Tmax.1995.01.01.tif')
r_ref <- r_ref %>% terra::crop(terra::ext(his_clm), snap = 'out')
r_ref <- terra::classify(r_ref, rcl = cbind(-Inf, -9990, NA))

#interpolate anomalies
rcp_anom <- intp_anomalies(his_clm  = his_clm,
                           rcp_clm  = rcp_clm,
                           anom_dir = anom_dir,
                           ref      = r_ref,
                           gcm_name = gcm_list[gcm_i],
                           rcp      = rcp,
                           varname  = varname,
                           period   = rcp_period)

#clean-up
rm(list = c('his_clm','rcp_clm','rcp_anom'))
gc(verbose = F, full = T, reset = T)
