#calculate and interpolate anomalies from GCM data
#JRV, Dec 2022

#load libraries
library(tidyverse)
library(terra)

#working directory
wd <- "~/common_data/esfg_cmip6"
raw_dir <- paste0(wd, "/raw")
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
#function to interpolate monthly anomalies
intp_anomalies <- function(his_clm, rcp_clm, anom_dir, ref, gcm_name, rcp, varname, period) {
  anom_fname <- paste0(anom_dir, "/", gcm_name, "_", rcp, "_r1i1p1f1_", varname, "_Africa_monthly_intp_anomaly_", min(period), "_", max(period), ".tif")
  if (!file.exists(anom_fname)) {
    r_anom <- c()
    for (i in 1:12) {
      #i <- 1
      cat("processing interpolation for month i=", i, "\n")
      
      #get climatology rasters
      avg_fut <- rcp_clm[[i]]
      avg_his <- his_clm[[i]]
      
      #calculate anomaly
      cat("calculating anomaly...\n")
      if (varname %in% c('tasmax','tasmin','tas')) {
        anom <- avg_fut - avg_his
      } else {
        #if precip is below zero make it zero
        avg_his[avg_his[] < 0] <- 0
        avg_fut[avg_fut[] < 0] <- 0
        
        #calculate anomaly in fraction
        anom <- (avg_fut - avg_his)/avg_his
        
        # Truncate the top 2% of anomaly values
        thr <- as.numeric(terra::global(x = anom, fun = stats::quantile, probs = 0.98, na.rm = T))
        anom[anom >= thr] <- thr
      }
      names(anom) <- "mean"
      
      #as data.frame
      crds <- anom %>% terra::as.data.frame(xy = T)
      
      #empty raster at GCM resolution
      empty <- anom
      terra::values(empty) <- NA
      
      #anomaly values
      anomalies_values <- unique(crds[,'mean'])
      
      #fit tps interpolation model
      cat("fitting thin plate spline\n")
      tps  <- fields::Tps(x = crds[,c('x','y')], Y = crds[,'mean'])
      
      #interpolate
      cat("interpolating onto the reference raster\n")
      #intp <- raster::interpolate(raster::raster(ref), tps)
      #intp <- terra::rast(intp) %>% terra::mask(mask = ref)
      intp <- terra::interpolate(object=terra::rast(ref), model=tps, fun=predict)
      intp <- intp %>%
        terra::mask(mask = ref)
      names(intp) <- paste0(sprintf("%02.0f",i))
      
      #clean-up
      rm(tps)
      gc(verbose=FALSE, full=TRUE)
      
      #append
      r_anom <- c(r_anom, intp)
    }
    
    #create raster, write
    r_anom <- terra::rast(r_anom)
    terra::writeRaster(r_anom, anom_fname)
  } else {
    r_anom <- terra::rast(anom_fname)
  }
  return(r_anom)
}


####
#loop rcp, variables, and period for given gcm
#add these rcps later "ssp126", "ssp370"
#add this variable later "tas"
#gcm_i <- 5
#rcp <- "ssp245" #"ssp585"
#varname <- "tasmin"#"tasmax", "pr"

for (rcp in c("ssp245", "ssp585")) {
  for (varname in c("tasmin", "tasmax", "pr")) {
    for (futperiod in c("near", "mid")) {
      #rcp <- "ssp585"; varname <- "tasmin"; futperiod <- "mid"
      cat("processing gcm=", gcm_list[gcm_i], "/ rcp=",rcp, "/ variable=", varname, "/ period=", futperiod, "\n")
      
      #define historical and future periods
      his_period <- 1995:2014
      if (futperiod == "near") {rcp_period <- 2021:2040}
      if (futperiod == "mid") {rcp_period <- 2041:2060}
      
      #data files
      his_file <- paste0(raw_dir, "/", gcm_list[gcm_i], "_historical_r1i1p1f1_", varname, "_Africa_daily.tif")
      rcp_file <- paste0(raw_dir, "/", gcm_list[gcm_i], "_", rcp, "_r1i1p1f1_", varname, "_Africa_daily.tif")
      
      #historical climatology
      his_clm <- calc_climatology(data_file=his_file, 
                                  period=his_period, 
                                  sce_lab="historical",
                                  gcm_name=gcm_list[gcm_i],
                                  varname=varname,
                                  mth_dir=mth_dir,
                                  clm_dir=clm_dir)$climatology
      
      #future climatology
      rcp_clm <- calc_climatology(data_file=rcp_file, 
                                  period=rcp_period, 
                                  sce_lab=rcp,
                                  gcm_name=gcm_list[gcm_i],
                                  varname=varname,
                                  mth_dir=mth_dir,
                                  clm_dir=clm_dir)$climatology
      
      #reference CHIRPS/CHIRTS raster
      r_ref <- terra::rast("~/common_data/chirts/Tmax/1995/Tmax.1995.01.01.tif")
      r_ref <- r_ref %>% terra::crop(terra::ext(his_clm), snap = 'out')
      r_ref[r_ref[] < -9990] <- NA
      
      #interpolate anomalies
      rcp_anom <- intp_anomalies(his_clm=his_clm, 
                                 rcp_clm=rcp_clm, 
                                 anom_dir=anom_dir, 
                                 ref=r_ref,
                                 gcm_name=gcm_list[gcm_i], 
                                 rcp=rcp, 
                                 varname=varname, 
                                 period=rcp_period)
    }
  }
}
