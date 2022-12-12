#download CHIRTS data
#JRV, 2022

#load needed libraries
library(geodata)
library(tidyverse)

#function to get CHIRTS data (Tmax, Tmin, RHum).
CHIRTS_data <- function(vars=c("RH", "Tmax", "Tmin"), years=1983:2016, basedir="~/common_data/chirts") {
  #create table of files and file names e.g., RH.1983.01.01.tif
  file_tb <- expand.grid(day=1:31, month=1:12, year=years, variable=vars)
  
  #file_tb <- expand.grid(day=1:10, month=1, year=1983, variable=c("RH"))
  filenames <- file_tb %>%
    base::split(.[,c("day","month","year","variable")]) %>%
    purrr::map(.f=function(.x) {
      fname <- paste0(.x$variable, ".", .x$year, ".", sprintf("%02.0f", .x$month), ".", sprintf("%02.0f", .x$day), ".tif")
      return(fname)
    }) %>%
    as.character(unlist(.))
  file_tb <- file_tb %>%
    dplyr::mutate(filename=filenames)
  
  #base url
  base_url <- "http://data.chc.ucsb.edu/products/CHIRTSdaily/v1.0/global_tifs_p05/"
  
  #get CHIRTS function
  .getCHIRTS <- function(.x, baseurl, bdir) {
    if (.x$variable == "RH") {varname <- "RHum"} else {varname <- .x$variable}
    this_url <- paste0(baseurl, varname, "/", .x$year, "/", .x$filename)
    outfolder <- paste0(bdir, "/", varname, "/", .x$year)
    if (!file.exists(outfolder)) {dir.create(outfolder, recursive=TRUE)}
    setwd(outfolder)
    status <- geodata:::.downloadDirect(url=this_url,
                                        filename=.x$filename, 
                                        unzip = FALSE, 
                                        quiet = FALSE, 
                                        mode = "wb", 
                                        cacheOK = FALSE)
    setwd(basedir)
    return(status)
  }
  
  #kick off download
  download_stat <- file_tb %>%
    base::split(.[,c("filename")]) %>%
    purrr::map(baseurl=base_url, bdir=basedir, .f=.getCHIRTS) %>%
    as.character(unlist(.))
  
  #append download status to table
  file_tb <- file_tb %>%
    dplyr::mutate(status=download_stat)
  
  #return values of function
  return(file_tb)
}

#run function
#out_table <- CHIRTS_data(vars=c("RH", "Tmax", "Tmin"), 
#                         years=1983:2016, 
#                         basedir="~/common_data/chirts")


