## ------------------------------------------ ##
## Download daily global CHIRPS
## By: H. Achicanoy
## By: Harold Achicanoy
## WUR & ABC
## Feb 2025
## ------------------------------------------ ##

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,lubridate,R.utils))

# Time frame
ini <- as.Date('2024-10-01')
end <- as.Date('2025-01-31')
dts <- seq(from = ini, to = end, by = 'day'); rm(ini, end)

# Output directory
Out  <- '//catalogue/WFP_ClimateRiskPr1/1.Data/Chirps'
dir.create(Out,F,T)

# Main function
getChirps <- function(date = dts[1]){
  # CHIRPS base URL
  chrps <- 'https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/tifs/p05'
  # Get day and year
  Day  <- date
  Year <- lubridate::year(Day)
  # Target file
  tfile <- paste0(chrps,'/',Year,'/chirps-v2.0.',gsub('-','.',Day,fixed=T),'.tif.gz')
  # Destination file
  dfile <- paste0(Out,'/',basename(tfile))
  # Raster file
  rfile <- gsub('.gz','',dfile,fixed = T)
  
  if(!file.exists(rfile)){
    
    # Downloading
    if(!file.exists(dfile)){
      tryCatch(expr = {
        utils::download.file(url = tfile, destfile = dfile)
      },
      error = function(e){
        cat(paste0(basename(tfile),' failed.\n'))
      })
    }
    
    # Unzip
    R.utils::gunzip(dfile)
    return(cat(paste0('Image ',basename(rfile),' processed correctly!!!\n')))
  } else {
    return(cat(paste0('Image ',basename(rfile),' already exists!\n')))
  }
  
}

# Loop through the dates
dts |> purrr::map(.f = getChirps)
