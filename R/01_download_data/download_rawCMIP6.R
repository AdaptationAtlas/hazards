library(devtools)
#install_github(c("SantanderMetGroup/loadeR.java",
#                   "SantanderMetGroup/climate4R.UDG",
#                   "SantanderMetGroup/loadeR",
#                   "SantanderMetGroup/transformeR",
#                   "SantanderMetGroup/visualizeR",
#                   "SantanderMetGroup/convertR",
#                   "SantanderMetGroup/climate4R.indices",
#                   "SantanderMetGroup/downscaleR"))

library(loadeR) 
library(climate4R.UDG) 
library(transformeR)
library(downscaleR)
library(climate4R.indices)
library(visualizeR)
library(rgdal)
library(tidyverse)
library(terra)

#data extraction requirements
lons <- c(-23, 59)  # Africa
lats <- c(-37, 40)   # Africa
season <- 1:12  # All year
years.hist <- 1995:2014
years.rcp <- 2021:2060

#GCMs of interest: MPI-ESM1-2-HR, ACCESS-ESM1-5, EC-Earth3-Veg, INM-CM5-0, MPI-ESM1-2-HR, MRI-ESM2-0

dataset.hist <- "CMIP6_ACCESS-ESM1-5_historical_r1i1p1f1"
dataset.rcp <- "CMIP6_ACCESS-ESM1-5_ssp585_r1i1p1f1"

#function to load data
load.data <- function (dset, years, var) loadGridData(dataset = dset, var = var,
                                                      years = years,
                                                      latLim = lats, lonLim = lons,
                                                      season = season) 
# Loading mean temperature, historical
tas_his <- load.data(dataset.hist, years.hist, "tas")

#base raster using characteristics
r1 <- terra::rast(nrows=length(clm_data$xyCoords$y), 
                  ncols=length(clm_data$xyCoords$x), 
                  nlyrs=1,
                  xmin=min(clm_data$xyCoords$x),
                  xmax=max(clm_data$xyCoords$x),
                  ymin=min(clm_data$xyCoords$y),
                  ymax=max(clm_data$xyCoords$y))

# function to return a SpatRaster object from the 'tas_his' list
makeRaster <- function(.x, clm_data, r1) {
  #cat(.x,"\n")
  datamat <- clm_data$Data[.x,,] #create matrix from array for day .x
  r2 <- terra::rast(datamat) %>%
    terra::flip(., direction='vertical') #create raster from matrix
  rout <- terra::rast(r1) #create output raster
  rout[] <- r2[] #transfer values
  names(rout) <- paste(clm_data$Dates$start[.x]) #name raster using date
  rm(r2) #cleanup
  return(rout)
}

#convert the data object/list into a raster
r_his <- 1:length(tas_his$Dates$start) %>%
  purrr::map(clm_data=tas_his, r1=r1, .f=makeRaster)
r_his <- terra::rast(r_his)
terra::writeRaster(r_his, "~/nfs/CMIP6_ACCESS-ESM1-5_historical_r1i1p1f1_daily.tif")
#plot(r_his)

#rcp data
tas_rcp <- load.data(dataset.rcp, years.rcp, "tas")

#convert the data object/list into a raster
r_rcp <- 1:length(tas_rcp$Dates$start) %>%
  purrr::map(clm_data=tas_rcp, r1=r1, .f=makeRaster)
r_rcp <- terra::rast(r_rcp)
terra::writeRaster(r_rcp, "~/nfs/CMIP6_ACCESS-ESM1-5_ssp585_r1i1p1f1_daily.tif")
#plot(r_rcp)

### below lines just as example
# # Loading minimum temperature
# y.tasmin <- load.data(dataset.obs, years.hist, "tasmin")
# x.tasmin <- load.data(dataset.hist, years.hist, "tasmin")
# newdata.tasmin <- load.data(dataset.rcp, years.rcp, "tasmin")
# 
# # Loading maximum temperature
# y.tasmax <- load.data(dataset.obs, years.hist, "tasmax")
# x.tasmax <- load.data(dataset.hist, years.hist, "tasmax")
# newdata.tasmax <- load.data(dataset.rcp, years.rcp, "tasmax")
# 
# 
