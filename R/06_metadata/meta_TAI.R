#Create Atlas metadata for TAI: Thornthwaite's Aridity Index
#also creates a single uploadable dataset for the S3/Google buckets.
#JRV, Jan 2023

#clean-up
rm(list = ls()) # Remove objects
g <- gc(reset = T); rm(g) # Empty garbage collector

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
his_dir <- paste0(wd, "/indices/historical/TAI")
ssp_dir <- paste0(wd, "/indices/", full_tb$fullname, "/TAI")
meta_dir <- paste0(wd, "/metadata")
if (!file.exists(meta_dir)) {dir.create(meta_dir)}
buck_dir <- paste0(wd, "/bucket_upload")
if (!file.exists(buck_dir)) {dir.create(buck_dir)}

#files to load
fls <- c("mean_year_masked.tif", "mean_year_masked_categorical.tif",
         "median_year_masked.tif", "median_year_masked_categorical.tif")

#this expands the list of files
fls <- c(sapply(paste0(his_dir, "/long_term_stats/"), FUN=function(i) paste0(i, fls)),
         sapply(paste0(ssp_dir, "/long_term_stats/"), FUN=function(i) paste0(i, fls)),
         sapply(paste0(ssp_dir, "/"), FUN=function(i) paste0(i, fls)))
fls <- fls[file.exists(fls)]

#first read datasets
r_data <- terra::rast(fls)

#give appropriate names to the layers
fnames <- gsub(pattern=paste0(wd, "/indices/"), replacement="", x=fls)
fnames <- gsub(pattern="_year_masked", replacement="", x=fnames)
fnames <- gsub(pattern="/TAI/long_term_stats", replacement="", x=fnames)
fnames <- gsub(pattern="/TAI", replacement="", x=fnames)
names(r_data) <- fnames

#write metadata
tai_meta <- atlas_metadata(data = r_data,
                           folder = meta_dir,
                           dataset.title_short = "TAI",
                           dataset.title_long = "Thornthwaite's Aridity Index",
                           dataset.desc = "Thornwaite aridity index = 100 * (d / n), where d = sum of monthly differences between prec and PET for months where precip < PET, and n = sum of monthly PET for those months. Note this index is calculated using the ENVIREM R package directly. ENVIREM takes a global solar radiation estimate (available at https://www.dropbox.com/sh/e5is592zafvovwf/AAAijCvHNiE4mYvYqWDpeJ3Ga/Global%20PET%20and%20Aridity%20Index?dl=0&lst=) to compute PET, and then uses that and precipitation to compute the Thornthwaite Aridity Index.",
                           dataset.author = "Ramirez-Villegas",
                           dataset.contact = "Julian Ramirez-Villegas",
                           dataset.contact_email = "j.r.villegas@cgiar.org",
                           dataset.pub_doi = list(dataset.pub_doi1="10.1175/JCLI-D-18-0698.1",
                                                  dataset.pub_doi2="10.1038/sdata.2015.66",
                                                  dataset.pub_doi3="10.1111/ecog.02880"),
                           dataset.data_doi = NA,
                           dataset.sourceurl = NA,
                           dataset.projecturl = "http://adaptationatlas.cgiar.org",
                           dataset.citation = "Ramirez-Villegas, J., Achicanoy, H. 2023. CMIP6 climate hazards: Thornthwaite's Aridity Index. CGIAR. Dataset.",
                           dataset.licence = "CC-BY-4.0",
                           file.filename = "TAI.tif",
                           file.format = "GeoTiff",
                           file.data_type = list(file.data_type1="float", 
                                                 file.data_type2="integer"),
                           file.no_value_data = "NA",
                           file.file_naming_convention = "for each layer: [scenario]_[model/ENSEMBLE]_[period]/[statistic]_[categorical].tif",
                           file.flags = NA,
                           variable.theme = "hazards",
                           variable.subtheme = "drought stress",
                           variable.name = "Thornthwaite's Aridity Index",
                           variable.subname = NA,
                           variable.commodity = NA,
                           variable.type = list(variable.type1="continuous", 
                                                variable.type2="categorical"),
                           variable.statistic = list(variable.statistic1="mean", 
                                                     variable.statistic2="median"),
                           variable.unit = "%",
                           method.analysis_type = "simple water balance, empirical equations",
                           method.description = "Thornwaite aridity index = 100 * (d / n), where d = sum of monthly differences between prec and PET for months where precip < PET, and n = sum of monthly PET for those months. Note this index is calculated using the ENVIREM R package directly. ENVIREM takes a global solar radiation estimate (available at https://www.dropbox.com/sh/e5is592zafvovwf/AAAijCvHNiE4mYvYqWDpeJ3Ga/Global%20PET%20and%20Aridity%20Index?dl=0&lst=) to compute PET, and then uses that and precipitation to compute the Thornthwaite Aridity Index.",
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
out_file <- file(paste0(meta_dir,"/TAI.txt"), open="w")
for (i in 1:ncol(tai_meta)) {
  writeLines(text=names(tai_meta)[i], con=out_file)
  writeLines(text=paste0(tai_meta[1,i], "\n"), con=out_file)
}
writeLines(text="Full list of layer names within file:", con=out_file)
writeLines(names(r_data), con=out_file)
close(con=out_file)

#write bucket transfer file
terra::writeRaster(r_data, paste0(buck_dir, "/TAI.tiff"), overwrite=TRUE)

