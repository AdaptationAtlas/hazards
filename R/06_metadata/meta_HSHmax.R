#Create Atlas metadata for HSH: Human Heat Stress index (monthly max). This function 
#also creates a single uploadable dataset for the S3/Google buckets.
#JRV, Jan 2023

#clean-up
rm(list = ls()) # Remove objects
g <- gc(reset = T) # Empty garbage collector

#source metadata function
source("https://raw.githubusercontent.com/AdaptationAtlas/metadata/main/R/metadata.R")

#libraries
library(terra)
library(tidyverse)

#working directory
wd <- "~/common_data/atlas_hazards/cmip6"

#gcms and scenarios
gcms <- c("ACCESS-ESM1-5", "MPI-ESM1-2-HR", "EC-Earth3", "INM-CM5-0", "MRI-ESM2-0", "ENSEMBLE")
ssps <- c("ssp245", "ssp585")
periods <- c("2021_2040", "2041_2060")
full_tb <- expand.grid(per=periods, ssp=ssps, gcm=gcms) %>%
  as.data.frame() %>%
  dplyr::mutate(fullname=paste0(.$ssp, "_", .$gcm, "_", .$per))

#directories
his_dir <- paste0(wd, "/indices/historical/HSH")
ssp_dir <- paste0(wd, "/indices/", full_tb$fullname, "/HSH")
meta_dir <- paste0(wd, "/metadata")
if (!file.exists(meta_dir)) {dir.create(meta_dir)}
buck_dir <- paste0(wd, "/bucket_upload")
if (!file.exists(buck_dir)) {dir.create(buck_dir)}

#files to load
fls <- c("max_year_masked.tif", "max_year_masked_categorical.tif", 
         "mean_year_masked.tif", "mean_year_masked_categorical.tif",
         "median_year_masked.tif", "median_year_masked_categorical.tif")

#this expands the list of files
fls <- c(sapply(paste0(his_dir, "/long_term_stats_max/"), FUN=function(i) paste0(i, fls)),
         sapply(paste0(ssp_dir, "/long_term_stats_max/"), FUN=function(i) paste0(i, fls)),
         sapply(paste0(ssp_dir, "_max/"), FUN=function(i) paste0(i, fls)))
fls <- fls[file.exists(fls)]

#first read datasets
r_data <- terra::rast(fls)

#give appropriate names to the layers
fnames <- gsub(pattern=paste0(wd, "/indices/"), replacement="", x=fls)
fnames <- gsub(pattern="/HSH_max", replacement="", x=fnames)
fnames <- gsub(pattern="_year_masked", replacement="", x=fnames)
fnames <- gsub(pattern="/HSH/long_term_stats_max", replacement="", x=fnames)
names(r_data) <- fnames

#write metadata
hsh_meta <- atlas_metadata(data = r_data,
                           folder = meta_dir,
                           dataset.title_short = "Maximum HSH",
                           dataset.title_long = "Human heat stress index, maximum value of the month",
                           dataset.desc = "The heat index (HI) was calculated from daily data and empirical equations based on Steadman (1979a; 1979b) relating the HI to dry-bulb temperature (°C) and relative humidity (%). Using the daily values the maximum of each month is first calculated, and following that we calculate the climatological mean of each month, and then the maximum, mean and median of the year.",
                           dataset.author = "Ramirez-Villegas",
                           dataset.contact = "Julian Ramirez-Villegas",
                           dataset.contact_email = "j.r.villegas@cgiar.org",
                           dataset.pub_doi = list(dataset.pub_doi1="10.1175/1520-0450(1979)018<0861:TAOSPI>2.0.CO;2",
                                                  dataset.pub_doi1="10.1175/1520-0450(1979)018<0874:TAOSPI>2.0.CO;2"),
                           dataset.data_doi = NA,
                           dataset.sourceurl = NA,
                           dataset.projecturl = "http://adaptationatlas.cgiar.org",
                           dataset.citation = "Ramirez-Villegas, J., Achicanoy, H., Thornton, P.K. 2023. CMIP6 climate hazards: human heat stress index. CGIAR. Dataset.",
                           dataset.licence = "CC-BY-4.0",
                           file.filename = "HSH_max.tif",
                           file.format = "GeoTiff",
                           file.data_type = list(file.data_type1="float", 
                                                 file.data_type2="integer"),
                           file.no_value_data = "NA",
                           file.file_naming_convention = "for each layer: [scenario]_[model/ENSEMBLE]_[period]/[statistic]_[categorical].tif",
                           file.flags = NA,
                           variable.theme = "hazards",
                           variable.subtheme = "heat stress",
                           variable.name = "Human heat stress index",
                           variable.subname = "Maximum of month",
                           variable.commodity = NA,
                           variable.type = list(variable.type1="continuous", 
                                                variable.type2="categorical"),
                           variable.statistic = list(variable.statistic1="mean", 
                                                     variable.statistic2="median", 
                                                     variable.statistic3="max"),
                           variable.unit = "dimensionless",
                           method.analysis_type = "empirical equation",
                           method.description = "The equation is as follows, HI = c1 + c2 * T + c3 * R + c4 * T * R + c5 * T^2 + c6 * R^2 + c7 * T^2 * R + c8 * T * R^2 + c9 * T^2 * R^2. T is the dry bulb temperature in ºC (taken to be the daily mean temperature). R is the relative humidity in %, and the values of the constants c1 to c9 are: c1=-8.78469475556; c2=1.61139411; c3=2.33854883889; c4=-0.14611605; c5=-0.012308094; c6=-0.0164248277778; c7=2.211732 x 10^-3; c8=7.2546 x 10^-4; c9=-3.582 x 10^-6.",
                           method.github = "https://github.com/AdaptationAtlas/hazards",
                           method.qual_indicator = NA,
                           method.qual_availability = NA,
                           temporal.resolution = "long-term average",
                           temporal.start_date = list(temporal.start_date1="1995", 
                                                      temporal.start_date2="2021", 
                                                      temporal.start_date3="2041"),
                           temporal.end_date = list(temporal.end_date1="2014", 
                                                    temporal.end_date2="2040", 
                                                    temporal.end_date3="2060"),
                           data.categorical_val = list(data.categorical_val1=1,
                                                       data.categorical_val2=2,
                                                       data.categorical_val3=3,
                                                       data.categorical_val4=4),
                           data.categorical_desc = list(data.categorical_desc1="No significant stress", 
                                                        data.categorical_desc2="Moderate", 
                                                        data.categorical_desc3="Severe", 
                                                        data.categorical_desc4="Extreme"),
                           data.shapefile_field = NA,
                           data.shapefile_type = NA,
                           data.shapefile_unit = NA,
                           data.shapefile_description = NA)

#write a plain text file with the information in the metadata object
out_file <- file(paste0(meta_dir,"/HSH_max.txt"), open="w")
for (i in 1:ncol(hsh_meta)) {
  writeLines(text=names(hsh_meta)[i], con=out_file)
  writeLines(text=paste0(hsh_meta[1,i], "\n"), con=out_file)
}
writeLines(text="Full list of layer names within file:", con=out_file)
writeLines(names(r_data), con=out_file)
close(con=out_file)

#write bucket transfer file
terra::writeRaster(r_data, paste0(buck_dir, "/HSH_max.tiff"), overwrite=TRUE)

