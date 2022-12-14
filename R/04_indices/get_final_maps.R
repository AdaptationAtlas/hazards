## Final maps to share
# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra))

OSys <<- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/catalogue/WFP_ClimateRiskPr1',
                'Windows' = '//CATALOGUE/WFP_ClimateRiskPr1')

gcms <- c('MPI-ESM1-2-HR','MRI-ESM2-0','EC-Earth3-Veg') # 'ACCESS-ESM1-5','INM-CM5-0'
prds <- c('2030s','2050s')
ssps <- c('ssp245','ssp585')
indx <- c('HSH','THI') # 'NDD','NTx40','NDWS','NWLD50','TAI'

# Africa shapefile
shp <- terra::vect(paste0(root,'/1.Data/shps/world_shapefile/Africa/africa.shp'))
ref <- terra::rast(paste0(root,'/10.New_Results/tmp.tif'))

# Template with water bodies
wbd <- terra::rast(paste0(root,'/10.New_Results/water_bodies.tif'))
wbd <- terra::crop(x = wbd, y = terra::ext(shp))

## Generate final numerical maps
numeric_map <- function(ind, prd, ssp){
  
  outfile <- paste0(root,'/10.New_Results/toShare/',ind,'/numeric/',ind,'_',prd,'_',ssp,'.tif')
  dir.create(dirname(outfile),F,T)
  
  # if(!file.exists(outfile)){
    
    in_dir <- paste0(root,'/10.New_Results/continental/numeric/',ind)
    fls <- list.files(path = in_dir, pattern = paste0('^',ind,'_',prd,'_',ssp), full.names = T)
    fls <- fls[grep(pattern = '.tif$', x = fls)]
    
    r <- fls %>%
      purrr::map(.f = function(fl){
        r <- terra::rast(fl)
        r <- terra::resample(r, ref)
      }) %>% terra::rast()  # Read GCM files
    r <- mean(r, na.rm = T) # Compute ensemble
    r[r == 0] <- NA
    r <- terra::focal(r, w = 9, fun = mean, na.policy = 'only', na.rm = T) # Fill gaps if exists
    r[is.na(r)] <- 0
    r <- terra::mask(r, shp) # Mask by continent
    wbd <- terra::resample(x = wbd, y = r) # Adjusting the extent
    r <- terra::mask(r, wbd, inverse = T)  # Mask removing water bodies
    
    if(ind %in% c('HSM','NTx30')){
      crp_msk <- terra::rast("//catalogue/WFP_ClimateRiskPr1/Atlas_1/1_hazards/hot_days_maize_s1/hts_maize_s1_hist.tif")
      crp_msk <- terra::resample(x = crp_msk, r)
      r <- terra::mask(r, crp_msk)
    }
    
    terra::writeRaster(r, outfile, overwrite = T)
    
  # }
  
  return('Done\n')
  
}

for(i in 1:length(prds)){
  for(j in 1:length(ssps)){
    for(l in 1:length(indx)){
      
      numeric_map(ind = indx[l], prd = prds[i], ssp = ssps[j])
      
    }
  }
}

## Generate final categorical maps
category_map <- function(ind, prd, ssp){
  
  infile <- paste0(root,'/10.New_Results/toShare/',ind,'/numeric/',ind,'_',prd,'_',ssp,'.tif')
  outfile <- paste0(root,'/10.New_Results/toShare/',ind,'/categorical/',ind,'_',prd,'_',ssp,'.tif')
  dir.create(dirname(outfile),F,T)
  
  # if(!file.exists(outfile)){
    
    r <- terra::rast(infile)
    
    if(ind %in% c('NDD','NDWS')){
      cat <- c(-Inf,  15,  1,
               15,  20,  2,
               20,  25,  3,
               25, Inf,  4)
      cat_m <- matrix(cat, ncol = 3, byrow = T)
      rc <- terra::classify(x = r, rcl = cat_m, right = F)
      cls <- data.frame(id = 1:4, class = c('No significant stress','Moderate','Severe','Extreme'))
      levels(rc) <- cls
    }
    if(ind == 'NTx40'){
      cat <- c(-Inf,   1,  1,
               1,   5,  2,
               5,  10,  3,
               10, Inf,  4)
      cat_m <- matrix(cat, ncol = 3, byrow = T)
      rc <- terra::classify(x = r, rcl = cat_m, right = F)
      cls <- data.frame(id = 1:4, class = c('No significant stress','Moderate','Severe','Extreme'))
      levels(rc) <- cls
    }
    if(ind == 'NWLD50'){
      cat <- c(-Inf,   2,  1,
               2,   5,  2,
               5,   8,  3,
               8, Inf,  4)
      cat_m <- matrix(cat, ncol = 3, byrow = T)
      rc <- terra::classify(x = r, rcl = cat_m, right = F)
      cls <- data.frame(id = 1:4, class = c('No significant stress','Moderate','Severe','Extreme'))
      levels(rc) <- cls
    }
    if(ind == 'TAI'){
      cat <- c(-Inf,  40,  1,
               40,  60,  2,
               60,  80,  3,
               80, Inf,  4)
      cat_m <- matrix(cat, ncol = 3, byrow = T)
      rc <- terra::classify(x = r, rcl = cat_m, right = F)
      cls <- data.frame(id = 1:4, class = c('No significant stress','Moderate','Severe','Extreme'))
      levels(rc) <- cls
    }
    if(ind == 'HSH'){
      # cat <- c(-Inf,  25,  1,
      #          25,  30,  2,
      #          30,  35,  3,
      #          35, Inf,  4)
      cat <- c(-Inf,  20,  1,
               20,  25,  2,
               25,  30,  3,
               30, Inf,  4)
      cat_m <- matrix(cat, ncol = 3, byrow = T)
      rc <- terra::classify(x = r, rcl = cat_m, right = F)
      cls <- data.frame(id = 1:4, class = c('No significant stress','Moderate','Severe','Extreme'))
      levels(rc) <- cls
    }
    if(ind == 'THI'){
      cat <- c(-Inf,   72,  1,
               72,   78,  2,
               78,  89,  3,
               89, Inf,  4)
      cat_m <- matrix(cat, ncol = 3, byrow = T)
      rc <- terra::classify(x = r, rcl = cat_m, right = F)
      cls <- data.frame(id = 1:4, class = c('No significant stress','Moderate','Severe','Extreme'))
      levels(rc) <- cls
    }
    if(ind %in% c('HSM','NTx30')){
      cat <- c(-Inf,  10,  1,
               10,  20,  2,
               20,  25,  3,
               25, Inf,  4)
      cat_m <- matrix(cat, ncol = 3, byrow = T)
      rc <- terra::classify(x = r, rcl = cat_m, right = F)
      cls <- data.frame(id = 1:4, class = c('No significant stress','Moderate','Severe','Extreme'))
      levels(rc) <- cls
    }
    
    terra::writeRaster(rc, outfile, overwrite = T)
    
  # }
  return('Done\n')
}

for(i in 1:length(prds)){
  for(j in 1:length(ssps)){
    for(l in 1:length(indx)){
      
      category_map(ind = indx[l], prd = prds[i], ssp = ssps[j])
      
    }
  }
}
