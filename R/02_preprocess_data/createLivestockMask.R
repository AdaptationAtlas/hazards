#Create Livestock VoP mask at CHIRPS resolution, combining multiple species
#JRV, Dec. 2022

#load packages
library(terra)
library(tidyverse)

#clean-up environment
rm(list=ls())
gc(verbose=FALSE, full=TRUE, reset=TRUE)

#working directory
wd <- "~/common_data/atlas_hazards/livestock_mask"
if (!file.exists(wd)) {dir.create(wd)}

#read Africa shapefile
r_msk <- terra::rast("~/common_data/atlas_hazards/roi/africa.tif")

#first list all VoP individual species files
lstk_dir <- "~/common_data/atlas_livestock/raw"
lstk_files <- list.files(lstk_dir, pattern="\\.tif")

#load them as raster
lstk_rs <- terra::rast(paste0(lstk_dir,"/",lstk_files)) %>%
  terra::crop(., r_msk)
lstk_rs <- terra::app(lstk_rs, fun=sum, na.rm=TRUE)
lstk_rs[lstk_rs[] == 0] <- NA

#resample resulting raster into CHIRPS resolution, use nn
chirps_rs <- terra::rast("~/common_data/chirps_wrld/chirps-v2.0.1995.01.01.tif") %>%
  terra::crop(., lstk_rs)
chirps_rs[chirps_rs[]<0] <- NA
chirps_rs[!is.na(chirps_rs[])] <- 1
lstk_rs <- terra::resample(lstk_rs, chirps_rs, method="bilinear")

#make zeroes as NA
lstk_rs[!is.na(lstk_rs[])] <- 1

#write raster
terra::writeRaster(lstk_rs, paste0(wd, "/livestock_allspecies_mask.tif"), overwrite=TRUE)


