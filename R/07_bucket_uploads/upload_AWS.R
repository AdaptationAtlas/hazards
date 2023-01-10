#Upload data files into AWS S3 bucket
#Note that "Access Key ID" and "Secret access key" are both in .Renviron

#clean-up
rm(list = ls()) # Remove objects
g <- gc(reset = T); rm(g) # Empty garbage collector

#load libraries
library(aws.s3)
library(tidyverse)

#bucket name
bk_name <- "s3://digital-atlas/"

#working directory
wd <- "~/common_data/atlas_hazards/cmip6"
meta_dir <- paste0(wd, "/metadata")
buck_dir <- paste0(wd, "/bucket_upload")

#list buckets
buck_list <- bucketlist()

#check if atlas bucket exists
atlas_bk <- bucket_exists(bucket=bk_name)

#list objects in bucket
bucket_df <- get_bucket_df(bucket=bk_name,
                           prefix="Updates_for_MVP_Release/1_hazards",
                           max = 20000) %>%
             as.data.frame()

#create metadata and datasets folders. Note that S3 buckets are flat (i.e., with no
#folder structure). Yet Amazon allows creating pseudo-folders to ease (visual) access
dset_ok <- put_folder(folder = "Updates_for_MVP_Release/1_hazards/datasets",
                      bucket = bk_name)
mdat_ok <- put_folder(folder = "Updates_for_MVP_Release/1_hazards/metadata",
                      bucket = bk_name)

#list objects in bucket_upload (datasets) folder
fls <- list.files(paste0(wd, "/bucket_upload"), full.names = TRUE)
1:length(fls) %>% 
  purrr::map(.f = function(i) {
      objname <- paste0("Updates_for_MVP_Release/1_hazards/datasets/",basename(fls[i]))
      fok <- object_exists(object = objname, bucket = bk_name, quiet = TRUE)
      if (!fok) {put_ok <- put_object(file = fls[i], object = objname, bucket = bk_name, verbose = FALSE, multipart=TRUE)}
    })


#list objects in metadata folder
fls <- list.files(paste0(wd, "/metadata"), full.names = TRUE)
1:length(fls) %>% 
  purrr::map(.f = function(i) {
    objname <- paste0("Updates_for_MVP_Release/1_hazards/metadata/",basename(fls[i]))
    fok <- object_exists(object = objname, bucket = bk_name, quiet = TRUE)
    if (!fok) {put_ok <- put_object(file = fls[i], object = objname, bucket = bk_name, verbose = FALSE)}
  })


