## Africa water bodies
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra))

root <- '//catalogue/Workspace14/WFP_ClimateRiskPr' # CIAT route

# CHIRPS template
tmp <- terra::rast('./africa.tif')

# Water bodies to raster
fls <- list.files(paste0(root,'/1.Data/shps/GLWD'), pattern = '_fixed.shp$', full.names = T)
wbd <- fls %>%
  purrr::map(.f = function(fl){
    shp <- terra::vect(fl)
    wbd <- terra::rasterize(x = shp, y = tmp)
    return(wbd)
  })
wbd <- terra::rast(wbd)
wbd <- terra::app(x = wbd, fun = sum, na.rm = T)

terra::writeRaster(wbd, './africa_wbd.tif', overwrite = T) # Manually uploaded to EiA
