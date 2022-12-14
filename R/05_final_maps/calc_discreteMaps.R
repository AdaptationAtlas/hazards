
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
