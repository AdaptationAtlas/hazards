#Download crop calendar
#JRV, Dec. 2022

#load packages
library(terra)
library(tidyverse)
library(geodata)
library(devtools)
#install_github(c("inbo/inborutils"))
library(inborutils)

#clean-up environment
rm(list=ls())
gc(verbose=FALSE, full=TRUE, reset=TRUE)

#working directory
wd <- "~/common_data/atlas_crop_calendar/raw_jagermeyr"
if (!file.exists(wd)) {dir.create(wd, recursive=TRUE)}

#set working directory
setwd(wd)

#download from zenodo, Jagermeyr et al. crop calendar (improved version of Sacks et al.)
download_zenodo(doi="10.5281/zenodo.5062513", path = ".", parallel = FALSE, quiet = FALSE)

####
#crop calendar from GAEZ+ 2015
wd <- "~/common_data/atlas_crop_calendar/raw_gaez2015"
if (!file.exists(wd)) {dir.create(wd, recursive=TRUE)}

#set working directory
setwd(wd)

#url of interest
this_url <- "https://mygeohub.org/publications/60/serve/1?el=1"

#download the data, file name given per download
status <- geodata:::.downloadDirect(url=this_url,
                                    filename="bundle.zip", 
                                    unzip = TRUE, 
                                    quiet = FALSE, 
                                    mode = "wb", 
                                    cacheOK = FALSE)

#unzip file
try(utils::unzip("GAEZ2015_Monthly.zip", exdir = "."), silent = TRUE)

