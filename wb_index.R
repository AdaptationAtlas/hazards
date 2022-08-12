# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
# .rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation

OSys <<- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/catalogue/WFP_ClimateRiskPr1',
                'Windows' = '//CATALOGUE/WFP_ClimateRiskPr1') #'Windows' = '//CATALOGUE/Workspace14/WFP_ClimateRiskPr'
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, sf, furrr, future, rgdal))


# Function to compute basic Agro-climatic indices
calc_AgrClm <- function(season = season, shp_fl = shp_fl, iso = iso){
  source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')
  outfile <- paste0(root, '/7.Results/Africa/',iso); if(!dir.exists(outfile)){dir.create(outfile, F, T)}
  
  calc_SHIMP <- function(TMAX, RH){SHI <- ifelse(TMAX >= 29 & RH > 50, 1, 0); return(SHI)}; calc_SHIMP <- compiler::cmpfun(calc_SHIMP)
  calc_THIMP <- function(TMAX, RH){THI <- (1.8 * TMAX + 32) - ((0.55 - 0.0055 * RH) * (1.8 * TMAX - 26.8)); THI_n <- ifelse(test = THI >= 79 & THI < 89 | THI >= 89, yes = 1, no = 0); return(THI_n)}; calc_THIMP <- compiler::cmpfun(calc_THIMP)
  calc_HSIMP <- function(TMAX, RH){
    # 2 ~ danger, 3 ~ emergency
    # HSI_n <- ifelse((TMAX <= 26 & RH > 70 |
    #                 TMAX <= 27 & RH >= 40 & RH < 85 |
    #                 TMAX <= 28 & RH < 85 |
    #                 TMAX <= 29 & RH < 60 |
    #                 TMAX <= 30 & RH < 40), 1, 0)
    HSI_n <- ifelse((TMAX <= 27 & RH >= 85 |
                       TMAX <= 28 & RH >= 85 |
                       TMAX <= 29 & RH >= 60 |
                       TMAX <= 30 & RH >= 40 |
                       TMAX > 30), 1, 0)
    return(HSI_n)
  }; calc_HSIMP <- compiler::cmpfun(calc_HSIMP)
  
  
 
  
  ## ROI: regions of interest
  shp <- terra::vect(shp_fl)
  
  ## Daily files
  # Precipitation
  # chr_pth <- '//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/Chirps'
  chr_pth <- paste0(root,'/1.Data/Chirps')
  chr_fls <- gtools::mixedsort(list.files(chr_pth, pattern = '*.tif$', full.names = T))
  chr_dts <- strsplit(x = chr_fls, split = 'chirps-v2.0.', fixed = T) %>% purrr::map(2) %>% unlist()
  chr_dts <- strsplit(x = chr_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
  chr_dts <- as.Date(gsub('.', '-', chr_dts, fixed = T))
  
  # root
  # era5Dir <- '//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/ERA5'
  # chirtsDir <- '//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/Chirts'
  era5Dir <- paste0(root,'/1.Data/ERA5')
  chirtsDir <- paste0(root,'/1.Data/Chirts')
  
  # Tmax
  tmx_pth <- paste0(chirtsDir,'/Tmax')
  tmx_fls <- gtools::mixedsort(list.files(tmx_pth, pattern = '*.tif$', full.names = T))
  tmx_dts <- strsplit(x = tmx_fls, split = 'chirts-v2.0.', fixed = T) %>% purrr::map(2) %>% unlist()
  tmx_dts <- strsplit(x = tmx_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
  tmx_dts <- gsub(pattern = '[.]', replacement = '',tmx_dts) %>% purrr::map(1) %>% unlist()
  tmx_dts <- as.Date(tmx_dts, "%Y%m%d")
  
  # Tmin
  tmn_pth <- paste0(chirtsDir,'/Tmin')
  tmn_fls <- gtools::mixedsort(list.files(tmn_pth, pattern = '*.tif$', full.names = T))
  tmn_dts <- strsplit(x = tmn_fls, split = 'chirts-v2.0.', fixed = T) %>% purrr::map(2) %>% unlist()
  tmn_dts <- strsplit(x = tmn_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
  tmn_dts <- gsub(pattern = '[.]', replacement = '',tmn_dts) %>% purrr::map(1) %>% unlist()
  tmn_dts <- as.Date(tmn_dts, "%Y%m%d")
  
  # Tmean
  tav_pth <- paste0(chirtsDir,'/Tmean')
  tav_fls <- gtools::mixedsort(list.files(tav_pth, pattern = '*.tif$', full.names = T))
  tav_dts <- strsplit(x = tav_fls, split = 'chirts-v2.0.', fixed = T) %>% purrr::map(2) %>% unlist()
  tav_dts <- strsplit(x = tav_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
  tav_dts <- gsub(pattern = '[.]', replacement = '',tav_dts) %>% purrr::map(1) %>% unlist()
  tav_dts <- as.Date(tav_dts, "%Y%m%d")
  
  # Solar radiation
  srd_pth <- paste0(era5Dir,'/solar_radiation_flux')
  srd_fls <- gtools::mixedsort(list.files(srd_pth, pattern = '*.nc$', full.names = T))
  srd_dts <- strsplit(x = srd_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
  srd_dts <- strsplit(x = srd_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
  srd_dts <- as.Date(srd_dts, "%Y%m%d")
  
  # Filtering days within the season
  yrs <- lubridate::year(tmx_dts)
  yrs <- names(table(yrs)[table(yrs) %in% 365:366])
  
  tmx_fls <- tmx_fls[lubridate::year(tmx_dts) %in% yrs]
  tmn_fls <- tmn_fls[lubridate::year(tmn_dts) %in% yrs]
  tav_fls <- tav_fls[lubridate::year(tav_dts) %in% yrs]
  srd_fls <- srd_fls[lubridate::year(srd_dts) %in% yrs]
  
  tmx_dts <- tmx_dts[lubridate::year(tmx_dts) %in% yrs]
  tmn_dts <- tmn_dts[lubridate::year(tmn_dts) %in% yrs]
  tav_dts <- tav_dts[lubridate::year(tav_dts) %in% yrs]
  srd_dts <- srd_dts[lubridate::year(srd_dts) %in% yrs]
  
  if(length(season) < 12){
    cnd <- lubridate::month(tmx_dts) %in% season # Days within the season
    yrs_dts <<- split(tmx_dts[cnd],cumsum(c(1,diff(tmx_dts[cnd])!=1)))
  } else {
    yrs <- lubridate::year(tmx_dts)
    grp <- with(rle(yrs), rep(seq_along(values), lengths))
    yrs_dts <<- split(tmx_dts, grp)
  }
  
  # Raster template
  tmp <- terra::rast('//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif')
  shp <- terra::vect(shp_fl)
  tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
  tmp[!is.na(tmp)] <- 1
  
  cat('..... Computing water balance model.\n')
  WTBL <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmn <- terra::rast(tmn_fls[tmn_dts %in% yrs_dts[[i]]])
      tmn <- tmn %>% terra::crop(terra::ext(tmp)) %>% terra::mask(tmp)
      
      tav <- terra::rast(tav_fls[tav_dts %in% yrs_dts[[i]]])
      tav <- tav %>% terra::crop(terra::ext(tmp)) %>% terra::mask(tmp)
      
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(tmp)) %>% terra::mask(tmp)
      
      
      srd <- terra::rast(srd_fls[srd_dts %in% yrs_dts[[i]]])
      srd <- srd %>% terra::resample(x = ., y = tmn[[1]]) %>% terra::mask(tmn[[1]])
      srd <- srd %>% terra::crop(terra::ext(tmn[[1]])) %>% terra::mask(tmn[[1]])
      srd <- srd/1000000
      
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(tmn[[1]])) %>% terra::mask(tmn[[1]])
      prc[prc == -9999] <- 0
      
      cat('..... Base de datos cargadas .\n')
      
      # Maximum evapotranspiration
      ETMAX <- terra::lapp(x = terra::sds(srd,tmn,tav,tmx), fun = peest)
      
      # Soil data
      scp <- terra::rast("//catalogue/workspace_cluster_14/WFP_ClimateRiskPr/Africa/soil_info/scp.tif")
      sst <- terra::rast("//catalogue/workspace_cluster_14/WFP_ClimateRiskPr/Africa/soil_info/ssat.tif")
      scp <- scp %>% terra::resample(tmp) %>% terra::mask(tmp) # Soil water capacity
      sst <- sst %>% terra::resample(tmp) %>% terra::mask(tmp) # Soil water saturation point
      
      # Compute water balance model
      AVAIL <<- tmp
      AVAIL[!is.na(AVAIL)] <- 0
      system.time( watbal <- 1:terra::nlyr(ETMAX) %>%
                     purrr::map(.f = function(i){
                       water_balance <- eabyep_calc(soilcp  = scp,
                                                    soilsat = sst,
                                                    avail   = AVAIL[[terra::nlyr(AVAIL)]],
                                                    rain    = prc[[i]],
                                                    evap    = ETMAX[[i]])
                       AVAIL <<- water_balance$Availability
                       return(water_balance)
                     })
      )
      cat('..... Compute water balance model ..... listo .\n')
      
      ERATIO  <- watbal %>% purrr::map('Eratio') %>% terra::rast()
      LOGGING <- watbal %>% purrr::map('Logging') %>% terra::rast()
      NDWS    <- terra::app(x = ERATIO, fun = function(ERATIO){ifelse(ERATIO < 0.5, 1, 0)}) %>% sum()
      NWLD50  <- sum(LOGGING > (sst*0.5))
      
      names(NDWS) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      names(NWLD50) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
       
      year <- 1995:2014
      
      terra::writeRaster(NDWS, paste0(outfile,"/NDWS_",year[i],"_", season,".tif" )) 
      terra::writeRaster(NWLD50, paste0(outfile,"/NWLD50_",year[i],"_", season,".tif" )) 
      
      
      cat('..... .........................    .\n')
      
      })
 }


seasons <- list(s1=1, s2=2, s3=3, s4=4, s5=5, s6=6, s7=7, s8=8, s9=9, s10=10, s11=11, s12=12) 

AFR <- raster::shapefile("//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/shps/world_shapefile/Africa/africa.shp")
country <- unique(c(AFR@data$ISO3))

# Loop through seasons
lapply (30:58, function(i){
    cat(paste0("Procesando pais :::: ",country[i], "\n" ))
    iso <- country[i]
    shp_fl <- paste0(root, '/1.Data/shps/',iso,'/',iso,'_GADM1.shp')
    lapply(1: length(seasons), function(j){
      cat(paste0("Procesando season :::: ",seasons[j], "\n" ))
      calc_AgrClm (season = seasons[j], shp_fl = shp_fl, iso = iso)
    })
    
    })
