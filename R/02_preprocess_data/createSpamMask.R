#Create SPAM mask at CHIRPS resolution, combining multiple crops
#JRV, Dec. 2022

#load packages
library(terra)
library(tidyverse)

#clean-up environment
rm(list=ls())
gc(verbose=FALSE, full=TRUE, reset=TRUE)

#working directory
wd <- "~/common_data/atlas_hazards/mapspam_mask"
if (!file.exists(wd)) {dir.create(wd)}

#read Africa shapefile
r_msk <- terra::rast("~/common_data/atlas_hazards/roi/africa.tif")

#first list all MapSPAM files labeled as "All technologies" (\\_A.tif)
spam_dir <- "~/common_data/mapspam_2017/raw"
spam_files <- list.files(spam_dir, pattern="\\_A.tif")

#filter by _H_ (harvested area)
spam_files <- spam_files[grep("_SSA_H_", spam_files)]

#load them as raster
spam_rs <- terra::rast(paste0(spam_dir,"/",spam_files)) %>%
  terra::crop(., r_msk)
spam_rs <- terra::app(spam_rs, fun=sum, na.rm=TRUE)

#resample resulting raster into CHIRPS resolution, use nn
chirps_rs <- terra::rast("~/common_data/chirps_wrld/chirps-v2.0.1995.01.01.tif") %>%
  terra::crop(., r_msk)
chirps_rs[chirps_rs[]<0] <- NA
chirps_rs[!is.na(chirps_rs[])] <- 1
spam_rs <- terra::resample(spam_rs, chirps_rs, method="bilinear")

#make zeroes as NA
spam_rs[!is.na(spam_rs[])] <- 1

#write raster
terra::writeRaster(spam_rs, paste0(wd, "/spam2017V2r1_allcrop_mask.tif"), overwrite=TRUE)


