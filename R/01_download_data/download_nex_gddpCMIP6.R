#download CMIP6 downscaled+bias corrected data from Nex-GDDP-CMIP6
#HAE, Dec 2024

#load libraries
library(pacman)
pacman::p_load(devtools)
# devtools::install_github("hllauca/RClimChange")
pacman::p_load(RClimChange)

#define working directory
wd <- '~/common_data/nex-gddp-cmip6'

#setup table
vars <- c('pr','tasmax','tasmin')
gcms <- c('ACCESS-ESM1-5','MPI-ESM1-2-HR','EC-Earth3','INM-CM5-0','MRI-ESM2-0')
ssps <- c('ssp126','ssp245','ssp370','ssp585')
prds <- c('2021_2040','2041_2060','2061_2080','2081_2100')
stp <- expand.grid(var = vars, gcm = gcms, ssp = ssps, period = prds,
                   stringsAsFactors = F) |>
  base::as.data.frame()
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

#download per row from setup table
for (i in 1:nrow(stp)) {
  
  #initial and ending years
  ini <- strsplit(stp$period[i], split = '_')[[1]][1] |> as.numeric()
  end <- strsplit(stp$period[i], split = '_')[[1]][2] |> as.numeric()
  # download
  RClimChange::nex_download(location = wd,
                            model    = stp$gcm[i],
                            scenario = stp$ssp[i],
                            variable = stp$var[i],
                            years    = ini:end,
                            version  = 'v1.1',
                            roi      = c(-26,58,-47,38), # xmin,xmax,ymin,ymax
                            method   = 'curl')
  
}
