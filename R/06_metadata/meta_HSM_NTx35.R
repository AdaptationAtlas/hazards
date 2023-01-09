#Create Atlas metadata for NDD: total number of dry days per month
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
his_dir <- paste0(wd, "/indices/historical/HSM_NTx35")
ssp_dir <- paste0(wd, "/indices/", full_tb$fullname, "/HSM_NTx35")
meta_dir <- paste0(wd, "/metadata")
if (!file.exists(meta_dir)) {dir.create(meta_dir)}
buck_dir <- paste0(wd, "/bucket_upload")
if (!file.exists(buck_dir)) {dir.create(buck_dir)}

#files to load
fls <- c("max_year_masked.tif", "max_year_masked_categorical.tif", 
         "mean_year_masked.tif", "mean_year_masked_categorical.tif",
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
fnames <- gsub(pattern="/HSM_NTx35/long_term_stats", replacement="", x=fnames)
fnames <- gsub(pattern="/HSM_NTx35", replacement="", x=fnames)
names(r_data) <- fnames

#write metadata
hsm_meta <- atlas_metadata(data = r_data,
                           folder = meta_dir,
                           dataset.title_short = "HSM_NTx35",
                           dataset.title_long = "Total number of heat stress days per month for maize",
                           dataset.desc = "The number of days with daily maximum temperatures above a given threshold during the growing season. To compute this, we use a maize crop calendar (and focus on the main growing season), and assume the flowering period occurs in the middle month of the season. Cropping calendars are taken from Jagermeyr et al. (2022), which is a modified version of Sacks et al. (2010). For maize, we use a temperature threshold of 35ºC, and assume that at this temperature heat stress starts affecting maize plants.",
                           dataset.author = "Ramirez-Villegas",
                           dataset.contact = "Julian Ramirez-Villegas",
                           dataset.contact_email = "j.r.villegas@cgiar.org",
                           dataset.pub_doi = list(dataset.pub_doi1="10.1038/s43016-021-00400-y",
                                                  dataset.pub_doi2="10.1175/JCLI-D-18-0698.1",
                                                  dataset.pub_doi3="10.1111/j.1466-8238.2010.00551.x"),
                           dataset.data_doi = NA,
                           dataset.sourceurl = NA,
                           dataset.projecturl = "http://adaptationatlas.cgiar.org",
                           dataset.citation = "Ramirez-Villegas, J., Achicanoy, H. 2023. CMIP6 climate hazards: number of dry days per month. CGIAR. Dataset.",
                           dataset.licence = "CC-BY-4.0",
                           file.filename = "HSM_NTx35.tif",
                           file.format = "GeoTiff",
                           file.data_type = list(file.data_type1="float", 
                                                 file.data_type2="integer"),
                           file.no_value_data = "NA",
                           file.file_naming_convention = "for each layer: [scenario]_[model/ENSEMBLE]_[period]/[statistic]_[categorical].tif",
                           file.flags = NA,
                           variable.theme = "hazards",
                           variable.subtheme = "heat stress",
                           variable.name = "number of heat stress days for maize",
                           variable.subname = "count of month",
                           variable.commodity = NA,
                           variable.type = list(variable.type1="continuous", 
                                                variable.type2="categorical"),
                           variable.statistic = list(variable.statistic1="mean", 
                                                     variable.statistic2="median", 
                                                     variable.statistic3="max"),
                           variable.unit = "dimensionless",
                           method.analysis_type = "count of days",
                           method.description = "The number of days with daily maximum temperatures above 35 Celsius degrees during the growing season. To compute this, we use a maize crop calendar (and focus on the main rainy season), and assume the flowering period occurs in the middle month of the season.",
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
out_file <- file(paste0(meta_dir,"/HSM_NTx35.txt"), open="w")
for (i in 1:ncol(hsm_meta)) {
  writeLines(text=names(hsm_meta)[i], con=out_file)
  writeLines(text=paste0(hsm_meta[1,i], "\n"), con=out_file)
}
writeLines(text="Full list of layer names within file:", con=out_file)
writeLines(names(r_data), con=out_file)
close(con=out_file)

#write bucket transfer file
terra::writeRaster(r_data, paste0(buck_dir, "/HSM_NTx35.tiff"), overwrite=TRUE)

