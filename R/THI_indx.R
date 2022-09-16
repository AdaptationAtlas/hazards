# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation

OSys <<- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/catalogue/WFP_ClimateRiskPr1',
                'Windows' = '//CATALOGUE/WFP_ClimateRiskPr1') #'Windows' = '//CATALOGUE/Workspace14/WFP_ClimateRiskPr'
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, sf, furrr, future))


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
  
  # Relative humidity
  rhy_pth <- paste0(chirtsDir,'/Rh')
  rhy_fls <- gtools::mixedsort(list.files(rhy_pth, pattern = '*.tif$', full.names = T))
  rhy_dts <- strsplit(x = rhy_fls, split = 'RH.', fixed = T) %>% purrr::map(2) %>% unlist()
  rhy_dts <- strsplit(x = rhy_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
  rhy_dts <- gsub(pattern = '[.]', replacement = '',rhy_dts) %>% purrr::map(1) %>% unlist()
  rhy_dts <- as.Date(rhy_dts, "%Y%m%d")
  
  # Filtering days within the season
  yrs <- lubridate::year(tmx_dts)
  yrs <- names(table(yrs)[table(yrs) %in% 365:366])
  
  tmx_fls <- tmx_fls[lubridate::year(tmx_dts) %in% yrs]
  rhy_fls <- rhy_fls[lubridate::year(rhy_dts) %in% yrs]
  
  tmx_dts <- tmx_dts[lubridate::year(tmx_dts) %in% yrs]
  rhy_dts <- rhy_dts[lubridate::year(rhy_dts) %in% yrs]
  
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
  
  IND <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      
      cat('..... Cortar bases de datos.\n')
      
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(tmp)) %>% terra::mask(tmp)
      
      rhy <- terra::rast(rhy_fls[rhy_dts %in% yrs_dts[[i]]])
      rhy <- rhy %>% terra::resample(.,tmx)%>% terra::crop(terra::ext(tmx)) %>% terra::mask(tmx)
      
      
      cat('..... Heat Stress livestock.\n')
      
      THI <- terra::lapp(x = terra::sds(tmx, rhy), fun = calc_THIMP)
      THI <- sum(THI)/terra::nlyr(tmx)
      THI <- THI %>% terra::mask(shp)
      
      year <- 1995:2014
      
      terra::writeRaster(THI, paste0(outfile,"/THI_",year[i],"_", season,".tif" ))
      
      cat('..................................   . \n')
    })
}


seasons <- list(s1=1, s2=2, s3=3, s4=4, s5=5, s6=6, s7=7, s8=8, s9=9, s10=10, s11=11, s12=12) 

AFR <- raster::shapefile("//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/shps/world_shapefile/Africa/africa.shp")
country <- unique(c(AFR@data$ISO3))

# Loop through seasons
lapply (1:5, function(i){
  cat(paste0("Procesando pais :::: ",country[i], "\n" ))
  iso <- country[i]
  shp_fl <- paste0(root, '/1.Data/shps/',iso,'/',iso,'_GADM1.shp')
  lapply(1: length(seasons), function(j){
    cat(paste0("Procesando season :::: ",seasons[j], "\n" ))
    calc_AgrClm (season = seasons[j], shp_fl = shp_fl, iso = iso)
  })
  
})
