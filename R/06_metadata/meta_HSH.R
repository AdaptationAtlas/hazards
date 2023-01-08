#Create Atlas metadata for HSH: Human Heat Stress index
#JRV, Jan 2023

#source metadata function
source("https://raw.githubusercontent.com/AdaptationAtlas/metadata/main/R/metadata.R")

#libraries
library(terra)
library(tidyverse)

#working directory
wd <- "~/common_data/atlas_hazards/cmip6"

#first read datasets
r_data <- terra::rast(paste0(wd, "/indices/..."))

hsh_meta <- atlas_metadata(data=...)