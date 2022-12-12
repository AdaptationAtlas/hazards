#download CMIP6 data from UDG mirror, using UDG tools
#JRV, Dec 2022

#load libraries
library(devtools)
library(loadeR) 
library(climate4R.UDG) 
library(transformeR)
library(downscaleR)
library(visualizeR)
library(rgdal)
library(tidyverse)
library(terra)

#in case the above packages do not exist, install as follows
#install_github(c("SantanderMetGroup/loadeR.java",
#                   "SantanderMetGroup/climate4R.UDG",
#                   "SantanderMetGroup/loadeR",
#                   "SantanderMetGroup/transformeR",
#                   "SantanderMetGroup/visualizeR",
#                   "SantanderMetGroup/convertR",
#                   "SantanderMetGroup/climate4R.indices",
#                   "SantanderMetGroup/downscaleR"))


#working directory
wd <- "~/common_data/esfg_cmip6/raw"
if (!file.exists(wd)) {dir.create(wd, recursive=TRUE)}

#lons <- c(-23, 59)  # Africa
#lats <- c(-37, 40)   # Africa
#years.hist <- 1995:2014
#years.rcp <- 2021:2060
#varname one of "tas","tasmin","tasmax","pr"
#rcp one of "ssp126", "ssp245", "ssp370", "ssp585"

#GCMs of interest: ACCESS-ESM1-5, EC-Earth3-Veg, INM-CM5-0, MPI-ESM1-2-HR, MRI-ESM2-0
dataset_list <- c("CMIP6_ACCESS-ESM1-5_scenario_r1i1p1f1",
                  "CMIP6_MPI-ESM1-2-HR_scenario_r1i1p1f1",
                  "CMIP6_EC-Earth3_scenario_r1i1p1f1",
                  "CMIP6_INM-CM5-0_scenario_r1i1p1f1",
                  "CMIP6_MRI-ESM2-0_scenario_r1i1p1f1")

#function to download CMIP6 data
downloadCMIP6 <- function(ds_name="CMIP6_ACCESS-ESM1-5_scenario_r1i1p1f1", rcp="ssp585", varname="tas", 
                          years.hist=1995:2014, years.rcp=2021:2060, lons=c(-23, 59), lats=c(-37, 40),
                          basedir) {
  #info
  cat("dataset=", ds_name, "/ rcp=", rcp, "/ variable=", varname, "\n")
  
  #names of datasets per CMIP6 Atlas (see https://github.com/SantanderMetGroup/ATLAS/)
  dataset.hist <- gsub(pattern="_scenario_", replacement="_historical_", x=ds_name)
  dataset.rcp <- gsub(pattern="_scenario_", replacement=paste0("_", rcp, "_"), x=ds_name)
  
  #function to load data
  load.data <- function (dset, years, var) loadGridData(dataset = dset, var = var,
                                                        years = years,
                                                        latLim = lats, lonLim = lons,
                                                        season = 1:12) 
  #file name, historical
  fname_his <- paste0(basedir, "/", dataset.hist,"_",varname,"_Africa_daily.tif")
  if (!file.exists(fname_his)) {
    #loading mean temperature, historical
    cat("downloading historical data, please wait...\n")
    data_his <- load.data(dataset.hist, years.hist, var=varname)
    
    #convert 'grid' to sp object
    r_his <- grid2sp(data_his)
    r_his <- terra::rast(r_his)
    
    #write file
    terra::writeRaster(r_his, filename=fname_his)
    
    #clean up
    rm(data_rcp)
    gc(verbose=FALSE, full=TRUE)
  } else {
    cat("historical data already exists, loading\n")
    r_his <- terra::rast(fname_his)
  }
  
  #rcp data, file name
  fname_rcp <- paste0(basedir, "/", dataset.rcp,"_",varname,"_Africa_daily.tif")
  if (!file.exists(fname_rcp)) {
    #loading mean temperature, rcp
    cat("downloading rcp data, please wait...\n")
    data_rcp <- load.data(dataset.rcp, years.rcp, var=varname)
    
    #convert 'grid' to sp object
    r_rcp <- grid2sp(data_rcp)
    r_rcp <- terra::rast(r_rcp)
    
    #write file
    terra::writeRaster(r_rcp, filename=fname_rcp)
    
    #clean up
    rm(data_rcp)
    gc(verbose=FALSE, full=TRUE)
  } else {
    cat("rcp data already exists, loading\n")
    r_rcp <- terra::rast(fname_rcp)
  }
  return(list(his=r_his, rcp=r_rcp))
}

#run function
for (i in 1:length(dataset_list)) {
  for (scenario in c("ssp126", "ssp245", "ssp370", "ssp585")) {
    for (varname in c("tas", "tasmin", "tasmax", "pr")) {
      cmip6_data <- downloadCMIP6(ds_name=dataset_list[i], 
                                  rcp=scenario, 
                                  varname=varname, 
                                  years.hist=1995:2014, 
                                  years.rcp=2021:2060, 
                                  lons=c(-23, 59), 
                                  lats=c(-37, 40), 
                                  basedir=wd)
      rm(cmip6_data)
      gc(verbose=FALSE, full=TRUE)
    }
  }
}
