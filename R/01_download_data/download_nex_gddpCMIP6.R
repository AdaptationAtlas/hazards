#download CMIP6 downscaled+bias corrected data from Nex-GDDP-CMIP6
#HAE, Dec 2024

#load libraries
library(pacman)
pacman::p_load(devtools)
# devtools::install_github("hllauca/RClimChange")
pacman::p_load(RClimChange,furrr,future,dplyr)

#define working directory
wd <- '~/common_data/nex-gddp-cmip6'; dir.create(wd, F, T)

#setup table
vars <- c('pr','tasmax','tasmin')
gcms <- c('ACCESS-ESM1-5','MPI-ESM1-2-HR','EC-Earth3','INM-CM5-0','MRI-ESM2-0')
ssps <- c('ssp126','ssp245','ssp370','ssp585')
prds <- c('2021_2040','2041_2060','2061_2080','2081_2100')
stp <- expand.grid(var = vars, gcm = gcms, ssp = ssps, period = prds,
                   stringsAsFactors = F) |>
  base::as.data.frame()
stp <- stp |> dplyr::arrange(var, gcm, ssp, period)
rm(vars, gcms, ssps, prds)

# # Download daily precipitation flux from BCC-CSM2-MR model
# # historical period
# nex_download(location = wd,
#              model    = gcm,
#              scenario = 'historical',
#              variable = 'pr',
#              years    = 1995:2014,
#              version  = 'v1.2',
#              roi      = c(-26,58,-47,38),
#              method   = 'curl')

# #download per row from setup table
# for (i in 1:nrow(stp)) {
#   
#   #initial and ending years
#   ini <- strsplit(stp$period[i], split = '_')[[1]][1] |> as.numeric()
#   end <- strsplit(stp$period[i], split = '_')[[1]][2] |> as.numeric()
#   # download
#   RClimChange::nex_download(location = wd,
#                             model    = stp$gcm[i],
#                             scenario = stp$ssp[i],
#                             variable = stp$var[i],
#                             years    = ini:end,
#                             version  = 'v1.1',
#                             roi      = c(-26,58,-47,38), # xmin,xmax,ymin,ymax
#                             method   = 'curl')
#   
# }

1:nrow(stp) |>
  purrr::map(.f = function(.x){
  
  #initial and ending years
  ini <- strsplit(stp$period[.x], split = '_')[[1]][1] |> as.numeric()
  end <- strsplit(stp$period[.x], split = '_')[[1]][2] |> as.numeric()
  # download
  RClimChange::nex_download(location = wd,
                            model    = stp$gcm[.x],
                            scenario = stp$ssp[.x],
                            variable = stp$var[.x],
                            years    = ini:end,
                            version  = 'v1.1',
                            roi      = NULL, # c(-26,58,-47,38), # xmin,xmax,ymin,ymax
                            method   = 'curl')
  #clean-up
  gc(verbose = F, full = T, reset = T)
  cat(paste0(stp$var[.x],'/',stp$gcm[.x],'/',stp$ssp[.x],'/',stp$period[.x],' ready.\n'))
})

plan(multisession, workers = 20)
furrr::future_map(.x = 1:nrow(stp), .f = function(.x){
  
  #initial and ending years
  ini <- strsplit(stp$period[.x], split = '_')[[1]][1] |> as.numeric()
  end <- strsplit(stp$period[.x], split = '_')[[1]][2] |> as.numeric()
  # download
  RClimChange::nex_download(location = wd,
                            model    = stp$gcm[.x],
                            scenario = stp$ssp[.x],
                            variable = stp$var[.x],
                            years    = ini:end,
                            version  = 'v1.1',
                            roi      = NULL, # c(-26,58,-47,38), # xmin,xmax,ymin,ymax
                            method   = 'curl')
  #clean-up
  gc(verbose = F, full = T, reset = T)
  
}, .progress = T)
plan(sequential)
gc(verbose = F, full = T, reset = T)
