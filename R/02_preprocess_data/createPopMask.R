#Create population mask at CHIRPS resolution, combining urban+rural
#JRV, Dec. 2022

#load packages
library(terra)
library(tidyverse)

#clean-up environment
rm(list=ls())
gc(verbose=FALSE, full=TRUE, reset=TRUE)

#working directory
wd <- "~/common_data/atlas_hazards/population_mask"
if (!file.exists(wd)) {dir.create(wd)}

#read Africa shapefile
r_msk <- terra::rast("~/common_data/atlas_hazards/roi/africa.tif")

#load population raster
pop_rs <- terra::rast("~/common_data/atlas_pop/raw/cell5m_afripop2020_urbanrural_ssa_popheadcount_total.tif") %>%
  terra::crop(., r_msk)
pop_rs[pop_rs[] == 0] <- NA

#resample resulting raster into CHIRPS resolution, use nn
chirps_rs <- terra::rast("~/common_data/chirps_wrld/chirps-v2.0.1995.01.01.tif") %>%
  terra::crop(., r_msk)
chirps_rs[chirps_rs[]<0] <- NA
chirps_rs[!is.na(chirps_rs[])] <- 1
pop_rs <- terra::resample(pop_rs, chirps_rs, method="bilinear")

#make zeroes as NA
pop_rs[!is.na(pop_rs[])] <- 1

#write raster
terra::writeRaster(pop_rs, paste0(wd, "/pop_mask.tif"), overwrite=TRUE)


