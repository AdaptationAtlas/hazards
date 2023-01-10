#Upload data files into Google Cloud bucket
#Note that for this script to work the machine where it is being run from
#will have to have access to the bucket. This can be done in the linux/macOS terminal
#by typing: gcloud auth login

#clean-up
rm(list = ls()) # Remove objects
g <- gc(reset = T); rm(g) # Empty garbage collector

#load libraries
library(tidyverse)

#bucket name
bk_name <- "gs://adaptation-atlas/"

#working directory
wd <- "~/common_data/atlas_hazards/cmip6"
meta_dir <- paste0(wd, "/metadata")
buck_dir <- paste0(wd, "/bucket_upload")

#list directories/files in bucket
bk_fls <- system(paste0("gsutil ls ", bk_name), intern=TRUE)

#list objects in bucket_upload (datasets) folder
fls <- list.files(paste0(wd, "/bucket_upload"), full.names = TRUE)
1:length(fls) %>% 
  purrr::map(.f = function(i) {
    objname <- paste0(bk_name, "cmip6_hazards/datasets/", basename(fls[i]))
    cp_file <- system(paste0("gsutil cp ", fls[i], " ", objname), intern=TRUE)
  })


#list objects in metadata folder
fls <- list.files(paste0(wd, "/metadata"), full.names = TRUE)
1:length(fls) %>% 
  purrr::map(.f = function(i) {
    objname <- paste0(bk_name, "cmip6_hazards/metadata/", basename(fls[i]))
    cp_file <- system(paste0("gsutil cp ", fls[i], " ", objname), intern=TRUE)
  })


