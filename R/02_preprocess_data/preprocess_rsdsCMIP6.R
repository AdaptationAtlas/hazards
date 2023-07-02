#preprocess raw CMIP6 data for downwards shortwave solar radiation (rsds)
#Julian Ramirez-Villegas
#Jul, 2023

#libraries
library(terra)
library(tidyverse)

#working directory
wd <- "~/common_data/esfg_cmip6/raw_rsds"

#GCMs of interest: ACCESS-ESM1-5, EC-Earth3-Veg, INM-CM5-0, MPI-ESM1-2-HR, MRI-ESM2-0
dataset_list <- c("rsds_day_ACCESS-ESM1-5",
                  "rsds_day_MPI-ESM1-2-HR",
                  "rsds_day_EC-Earth3",
                  "rsds_day_INM-CM5-0",
                  "rsds_day_MRI-ESM2-0")

#function to download CMIP6 data r1i1p1f1
read_rsdsCMIP6 <- function(ds_name="rsds_day_ACCESS-ESM1-5", variant="r1i1p1f1", rcp="ssp126",
                           years.hist=1995:2014, years.rcp=2021:2100, lons=c(-23, 59), lats=c(-37, 40),
                           basedir) {
  
  #info
  cat("dataset=", ds_name, "/ rcp=", rcp, "/ variable=", variant, "\n")
  
  #read data: here i need to build a function that simply lists files and returns the 2 lists
  #this is for ACCESS, but i'll need it for each GCM
  gcm.data <- getGCMFileList(ds_name, rcp, variant)
  hist.data <- gcm.data$his
  rcp.data <- gcm.data$rcp
  
  #load historical data
  fname_his <- paste0(basedir, "/CMIP6_", ds_name, "_historical_rsds_Africa_daily.tif")
  fname_his <- gsub(pattern="_rsds_day", replacement="", x=fname_his)
  if (!file.exists(fname_his)) {
    cat("loading historical data, please wait...\n")
    r_his <- rast(paste0(wd,"/", hist.data))
    years <- lubridate::year(terra::time(r_his)) %in% years.hist
    r_his <- r_his[[years]]
    r_his <- r_his %>%
      terra::rotate(.) %>%
      terra::crop(., terra::ext(lons, lats))
    names(r_his) <- paste(terra::time(r_his))
    
    #unit conversion #w/m2 = J/m2/s / 1000000 * 86400 = MJ/m2/day
    r_his <- r_his * 24 * 3600 / 1000000
    
    #write file, this is in crop model units (MJ/m2/day)
    terra::writeRaster(r_his, filename=fname_his)
    
    #clean up
    rm(r_his)
    gc(verbose=FALSE, full=TRUE)
  } else {
    cat("historical data already exists, loading\n")
    r_his <- terra::rast(fname_his)
  }
  
  #load historical data
  fname_rcp <- paste0(basedir, "/CMIP6_", ds_name, "_", rcp, "_rsds_Africa_daily.tif")
  fname_rcp <- gsub(pattern="_rsds_day", replacement="", x=fname_rcp)
  if (!file.exists(fname_rcp)) {
    cat("loading rcp data, please wait...\n")
    r_rcp <- rast(paste0(wd,"/", rcp.data))
    years <- lubridate::year(terra::time(r_rcp)) %in% years.rcp
    r_rcp <- r_rcp[[years]]
    r_rcp <- r_rcp %>%
      terra::rotate(.) %>%
      terra::crop(., terra::ext(lons, lats))
    
    #unit conversion #w/m2 = J/m2/s / 1000000 * 86400 = MJ/m2/day
    r_rcp <- r_rcp * 24 * 3600 / 1000000
    
    #write file, this is in crop model units (MJ/m2/day)
    terra::writeRaster(r_rcp, filename=fname_rcp)
    
    #clean up
    rm(r_rcp)
    gc(verbose=FALSE, full=TRUE)
  } else {
    cat("rcp data already exists, loading\n")
    r_rcp <- terra::rast(fname_rcp)
  }
  return(list(his=r_his, rcp=r_rcp))
}

###get list of gcm files
getGCMFileList <- function(ds_name="rsds_day_ACCESS-ESM1-5", rcp="ssp585", variant="r1i1p1f1") {
  out.list <- list()
  if (ds_name == "rsds_day_ACCESS-ESM1-5") {
    out.list$his <- c(paste0(ds_name, "_historical_", variant, "_gn_19500101-19991231.nc"),
                   paste0(ds_name, "_historical_", variant, "_gn_20000101-20141231.nc"))
    out.list$rcp <- c(paste0(ds_name, "_", rcp, "_", variant, "_gn_20150101-20641231.nc"),
                  paste0(ds_name, "_", rcp, "_", variant, "_gn_20650101-21001231.nc"))
  } else if (ds_name == "rsds_day_MPI-ESM1-2-HR") {
    out.list$his <- c(paste0(ds_name, "_historical_", variant, "_gn_19950101-19991231.nc"),
                      paste0(ds_name, "_historical_", variant, "_gn_20000101-20041231.nc"),
                      paste0(ds_name, "_historical_", variant, "_gn_20050101-20091231.nc"),
                      paste0(ds_name, "_historical_", variant, "_gn_20100101-20141231.nc"))
    out.list$rcp <- c(paste0(ds_name, "_", rcp, "_", variant, "_gn_20150101-20191231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20200101-20241231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_202500101-20291231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20300101-20341231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20350101-20391231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20400101-20441231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20450101-20491231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20500101-20541231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20550101-20591231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20600101-20641231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20650101-20691231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20700101-20741231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20750101-20791231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20800101-20841231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20850101-20891231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20900101-20941231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20950101-20991231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_21000101-21001231.nc"))
  } else if (ds_name == "rsds_day_EC-Earth3") {
    #historical
    hdates1 <- seq(from=as.Date('1995-01-01'), to=as.Date('2014-01-01'), by='year') %>%
      paste0(.) %>% gsub("-", "", .)
    hdates2 <- seq(from=as.Date('1995-12-31'), to=as.Date('2014-12-31'), by='year') %>%
      paste0(.) %>% gsub("-", "", .)
    hdates <- data.frame(d1=hdates1, d2=hdates2) %>%
      dplyr::mutate(fullname = paste0(ds_name, "_historical_", variant, "_gr_",d1, "-", d2, ".nc"))
    out.list$his <- paste0(hdates$fullname)
    
    #rcp
    rdates1 <- seq(from=as.Date('2015-01-01'), to=as.Date('2100-01-01'), by='year') %>%
      paste0(.) %>% gsub("-", "", .)
    rdates2 <- seq(from=as.Date('2015-12-31'), to=as.Date('2100-12-31'), by='year') %>%
      paste0(.) %>% gsub("-", "", .)
    rdates <- data.frame(d1=rdates1, d2=rdates2) %>%
      dplyr::mutate(fullname = paste0(ds_name, "_", rcp, "_", variant, "_gr_",d1, "-", d2, ".nc"))
    out.list$rcp <- paste0(rdates$fullname)
  } else if (ds_name == "rsds_day_INM-CM5-0") {
    out.list$his <- c(paste0(ds_name, "_historical_", variant, "_gr1_19500101-19991231.nc"),
                      paste0(ds_name, "_historical_", variant, "_gr1_20000101-20141231.nc"))
    out.list$rcp <- c(paste0(ds_name, "_", rcp, "_", variant, "_gr1_20150101-20641231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gr1_20650101-21001231.nc"))
  } else if (ds_name == "rsds_day_MRI-ESM2-0") {
    out.list$his <- c(paste0(ds_name, "_historical_", variant, "_gn_19500101-19991231.nc"),
                      paste0(ds_name, "_historical_", variant, "_gn_20000101-20141231.nc"))
    out.list$rcp <- c(paste0(ds_name, "_", rcp, "_", variant, "_gn_20150101-20641231.nc"),
                      paste0(ds_name, "_", rcp, "_", variant, "_gn_20650101-21001231.nc"))
  }
  return(out.list)
}


#run function
for (i in 1:length(dataset_list)) {
  #i <- 3
  for (scenario in c("ssp126", "ssp245", "ssp370", "ssp585")) {
    rsds_data <- read_rsdsCMIP6(ds_name=dataset_list[i], 
                                variant="r1i1p1f1",
                                rcp=scenario, 
                                years.hist=1995:2014, 
                                years.rcp=2061:2100,
                                lons=c(-23, 59), 
                                lats=c(-37, 40), 
                                basedir=wd)
    rm(rsds_data)
    gc(verbose=FALSE, full=TRUE)
  }
}
