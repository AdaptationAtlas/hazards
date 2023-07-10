#Create Atlas metadata for PTOT: total monthly precipitation
#also creates a single uploadable dataset for the S3/Google buckets.
#JRV, July 2023

# Clean-up
rm(list = ls()) # Remove objects
gc(reset = T) # Empty garbage collector

# Libraries
library(terra)
library(tidyverse)

# General config
skip_preparation <- FALSE

# ------------------------------------------------------------
# Run calc_LongTermStats: run this section of the code if the
# long-term statistics (long-term mean, etc) need to be rerun
if (!skip_preparation) {
  source("~/Repositories/hazards/R/05_final_maps/calc_LongTermStats.R")
  stp <- expand.grid(sce = sce_list, prd = period_list) %>% 
    as.data.frame() %>%
    dplyr::mutate(sce_prd = paste0(sce, "-", prd)) %>%
    dplyr::filter(!sce_prd %in% c("historical-near", "historical-mid", "ssp245-hist", "ssp585-hist")) %>%
    dplyr::select(-sce_prd)
  1:nrow(stp) %>%
    purrr::map(.f=function(i) {
      continuous_map(index="PTOT", HS.stat=NULL, period=stp$prd[i], scenario=stp$sce[i], domean=TRUE, domedian=TRUE, domax=TRUE, doensemble=TRUE, omitcalendar=FALSE)
    })
}

# ------------------------------------------------------------
# Run calc_discreteMaps: run this section of the code if the
# discrete maps (categorical maps) need to be rerun.
if (!skip_preparation) {
  rm(list = ls()) # Remove objects
  gc(reset = T) # Empty garbage collector
  source("~/Repositories/hazards/R/05_final_maps/calc_discreteMaps.R")
  stp <- expand.grid(sce=sce_list[2:3], prd = period_list[2:3], gcm = gcm_list, stat=stat_list) %>%
    as.data.frame() %>%
    rbind(data.frame(sce="historical", prd="hist", gcm=NA, stat=stat_list), .)
  1:nrow(stp) %>%
    purrr::map(.f=function(i) {
      category_map(index="PTOT", HS.stat=NULL, period=stp$prd[i], scenario=stp$sce[i], gcm=stp$gcm[i], stat=stp$stat[i])
    })
}

# ------------------------------------------------------------
# Compile meta-data
# Clean-up
rm(list = ls()) # Remove objects
gc(reset = T) # Empty garbage collector

# Source metadata and other functions
source("https://raw.githubusercontent.com/AdaptationAtlas/metadata/main/R/metadata.R")

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
his_dir <- paste0(wd, "/indices/historical/PTOT")
ssp_dir <- paste0(wd, "/indices/", full_tb$fullname, "/PTOT")
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
fnames <- gsub(pattern="/PTOT/long_term_stats", replacement="", x=fnames)
fnames <- gsub(pattern="/PTOT", replacement="", x=fnames)
names(r_data) <- fnames

#write metadata
ptot_meta <- atlas_metadata(data = r_data,
                           folder = meta_dir,
                           dataset.title_short = "PTOT",
                           dataset.title_long = "Total precipitation per month (mm)",
                           dataset.desc = "Sum of total precipitation month and per year for historical and future periods, calculated using CHIRPS daily precipitation data.",
                           dataset.author = "Ramirez-Villegas",
                           dataset.contact = "Julian Ramirez-Villegas",
                           dataset.contact_email = "j.r.villegas@cgiar.org",
                           dataset.pub_doi = "10.1038/sdata.2015.66",
                           dataset.data_doi = NA,
                           dataset.sourceurl = NA,
                           dataset.projecturl = "http://adaptationatlas.cgiar.org",
                           dataset.citation = "Ramirez-Villegas, J., Achicanoy, H. 2023. CMIP6 climate hazards: total precipitation per month. CGIAR. Dataset.",
                           dataset.licence = "CC-BY-4.0",
                           file.filename = "PTOT.tif",
                           file.format = "GeoTiff",
                           file.data_type = list(file.data_type1="float", 
                                                 file.data_type2="integer"),
                           file.no_value_data = "NA",
                           file.file_naming_convention = "for each layer: [scenario]_[model/ENSEMBLE]_[period]/[statistic]_[categorical].tif",
                           file.flags = NA,
                           variable.theme = "hazards",
                           variable.subtheme = "drought and excess water stress",
                           variable.name = "total precipitation",
                           variable.subname = "total precipitation amount",
                           variable.commodity = NA,
                           variable.type = list(variable.type1="continuous", 
                                                variable.type2="categorical"),
                           variable.statistic = list(variable.statistic1="mean", 
                                                     variable.statistic2="median", 
                                                     variable.statistic3="max"),
                           variable.unit = "mm",
                           method.analysis_type = "sum of daily precipitation",
                           method.description = "Sum daily rainfall (i.e., accumulated daily precipitation) per month.",
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
out_file <- file(paste0(meta_dir,"/PTOT.txt"), open="w")
for (i in 1:ncol(ptot_meta)) {
  writeLines(text=names(ptot_meta)[i], con=out_file)
  writeLines(text=paste0(ptot_meta[1,i], "\n"), con=out_file)
}
writeLines(text="Full list of layer names within file:", con=out_file)
writeLines(names(r_data), con=out_file)
close(con=out_file)

#write bucket transfer file
terra::writeRaster(r_data, paste0(buck_dir, "/PTOT.tiff"), overwrite=TRUE)

