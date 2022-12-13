## Africa mask
## By: H. Achicanoy
## December, 2022

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(if(!require(pacman)){install.packages('pacman');library(pacman)} else {library(pacman)})
suppressMessages(pacman::p_load(tidyverse,raster,terra,sp,geodata,rnaturalearthdata,rnaturalearth))

# Root directory
root <- '/home/jovyan/common_data/atlas_hazards'

## Shapefile
# Output file
out  <- paste0(root,'/roi/africa.gpkg')
dir.create(dirname(out), F, T)
if(!file.exists(out)){
  # World shapefile
  wrld <- rnaturalearth::ne_countries(returnclass = 'sf', scale = 50)
  # Africa shapefile
  afrc <- filter(wrld, region_un == 'Africa')
  # Save shapefile in .gpkg file
  afrc <- terra::vect(afrc)
  terra::writeVector(afrc, out)
}

## Raster
tif  <- paste0(root,'/roi/africa.tif')
if(!file.exists(tif)){
  # Africa shapefile
  afrc <- terra::vect(out)
  # CHIRPS template
  ref <- terra::rast('/home/jovyan/common_data/chirps_wrld/chirps-v2.0.1981.01.02.tif')
  ref <- ref %>% terra::crop(terra::ext(afrc))
  ref <- terra::rasterize(afrc, ref)
  terra::writeRaster(ref, tif)
}

# Increase extent in the northern area
r <- terra::rast(tif)
r <- terra::extend(x = r, y = terra::ext(c(-25.3499976955354, 57.8000035434961,
                                            -46.9500014446676, 40)))
terra::writeRaster(r, tif, overwrite = T)
