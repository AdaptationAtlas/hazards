# Download CMIP6 downscaled+bias corrected data from Nex-GDDP-CMIP6
# By: H. Achicanoy
# Alliance Bioversity International & CIAT, 2025

# R options
options(warn = -1, scipen = 999)

# Load libraries and functions
if(!require(pacman)){install.packages('pacman');library(pacman)} else {library(pacman)}
pacman::p_load(purrr,furrr,future,dplyr,httr)
grep2 <- Vectorize(grep, 'pattern')

urlFileExist <- function(url){
  HTTP_STATUS_OK <- 200
  hd <- httr::HEAD(url)
  status <- hd$all_headers[[1]]$status
  list(exists = status == HTTP_STATUS_OK)
}

scenario <- 'historical'

# # Available files to download
# fls <- readLines('https://nex-gddp-cmip6.s3-us-west-2.amazonaws.com/index_v1.1_md5.txt')
# nms <- strsplit(fls, split = '/') |> purrr::map(6) |> unlist()

# Root URL for downloads
root <- 'https://nex-gddp-cmip6.s3.us-west-2.amazonaws.com/NEX-GDDP-CMIP6'

# Filters to apply
gcms <- c('ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CMCC-ESM2','EC-Earth3','EC-Earth3-Veg-LR','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','TaiESM1')
vars <- c('pr','tasmax','tasmin','hurs','rsds')
if (scenario == 'future') {
  ssps <- c('ssp126','ssp245','ssp370','ssp585')
  yrs <- 2021:2100
} else {
  if (scenario == 'historical') {
    ssps <- 'historical'
    yrs <- 1981:1994
  }
}

# Setup table
stp <- base::expand.grid(gcm = gcms, ssp = ssps, var = vars, yr = yrs, stringsAsFactors = F) |>
  base::as.data.frame(); rm(gcms, ssps, vars, yrs)
# Available files to download
wd <- '/home/jovyan/common_data/nex-gddp-cmip6_raw'
dir.create(wd, F, T)
if (scenario == 'future') {
  outfile <- file.path(wd,'cmip6_future_files_to_download.csv')
} else {
  if (scenario == 'historical') {
    outfile <- file.path(wd,'cmip6_short-baseline_files_to_download.csv')
  }
}
if(!file.exists(outfile)) {
  plan(multisession, workers = 30)
  available_files <-  1:nrow(stp) |>
    furrr::future_map(.f = function(i) {
      input_dir <- paste0(root,'/',stp$gcm[i],'/',stp$ssp[i],'/r1i1p1f1/',stp$var[i])
      input_file <- c(paste0(stp$var[i],'_day_',stp$gcm[i],'_',stp$ssp[i],'_r1i1p1f1_gn_',stp$yr[i],'_v1.2.nc'),
                      paste0(stp$var[i],'_day_',stp$gcm[i],'_',stp$ssp[i],'_r1i1p1f1_gr_',stp$yr[i],'_v1.2.nc'),
                      paste0(stp$var[i],'_day_',stp$gcm[i],'_',stp$ssp[i],'_r1i1p1f1_gr1_',stp$yr[i],'_v1.2.nc'),
                      paste0(stp$var[i],'_day_',stp$gcm[i],'_',stp$ssp[i],'_r1i1p1f1_gn_',stp$yr[i],'_v1.1.nc'),
                      paste0(stp$var[i],'_day_',stp$gcm[i],'_',stp$ssp[i],'_r1i1p1f1_gr_',stp$yr[i],'_v1.1.nc'),
                      paste0(stp$var[i],'_day_',stp$gcm[i],'_',stp$ssp[i],'_r1i1p1f1_gr1_',stp$yr[i],'_v1.1.nc'),
                      paste0(stp$var[i],'_day_',stp$gcm[i],'_',stp$ssp[i],'_r1i1p1f1_gn_',stp$yr[i],'.nc'),
                      paste0(stp$var[i],'_day_',stp$gcm[i],'_',stp$ssp[i],'_r1i1p1f1_gr_',stp$yr[i],'.nc'),
                      paste0(stp$var[i],'_day_',stp$gcm[i],'_',stp$ssp[i],'_r1i1p1f1_gr1_',stp$yr[i],'.nc'))
      
      #available files to download
      input_file_avl <- input_file[which(unlist(lapply(file.path(input_dir, input_file), urlFileExist)))]
      cndt2 <- grep(pattern = 'v1.2', x = input_file_avl)
      cndt1 <- grep(pattern = 'v1.1', x = input_file_avl)
      if (length(cndt2) > 0) {
        input_file_dwn <- input_file_avl[cndt2]
      } else {
        if (length(cndt1) > 0) {
          input_file_dwn <- input_file_avl[cndt1]
        } else {
          input_file_dwn <- input_file_avl
        }
      }
      res <- data.frame(pth_dir = input_dir, file = input_file_dwn)
      return(res)
    }, .progress = T) |> dplyr::bind_rows()
  plan(sequential)
  gc(F,T,T)
  
  stp <- cbind(stp, available_files); rm(available_files)
  utils::write.csv(x = stp, file = outfile, row.names = F)
  
} else {
  stp <- utils::read.csv(outfile)
}

# Download files (serial)
1:nrow(stp) |>
  purrr::map(.f = function(i) {
    
    outd <- paste0(wd,'/',stp$var[i],'/',stp$ssp[i],'/',stp$gcm[i])
    dir.create(outd,F,T)
    outfile <- file.path(outd,stp$file[i])
    
    if (!file.exists(outfile) || file.size(outfile) < 1e8) {
      download.file(url = file.path(stp$pth_dir[i], stp$file[i]), destfile = outfile, method = 'curl')
    } else {
      cat(stp$file[i],'downloaded!\n')
    }
    return('Done.\n')
  })
gc(F, T, T)

# # Download files (parallel)
# plan(multisession, workers = 30)
# 1:nrow(stp) |>
#   furrr::future_map(.f = function(i) {
#     
#     outd <- paste0(wd,'/',stp$var[i],'/',stp$ssp[i],'/',stp$gcm[i])
#     dir.create(outd,F,T)
#     outfile <- file.path(outd,stp$file[i])
#     
#     if (!file.exists(outfile) || file.size(outfile) < 1e8) {
#       download.file(url = file.path(stp$pth_dir[i], stp$file[i]), destfile = outfile, method = 'curl')
#     } else {
#       cat(stp$file[i],'downloaded!\n')
#     }
#     return('Done.\n')
#   }, .progress = T)
# plan(sequential)
# gc(F, T, T)