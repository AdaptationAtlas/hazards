## Check precipitation/temperature daily files
## By: H. Achicanoy
## January, 2023

# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse,terra,gtools,lubridate,future,furrr,future.apply,parallel))
grep2 <- Vectorize(grep, vectorize.args = 'pattern')

root <- '/home/jovyan/common_data'

# Setup
gcms <- c('ACCESS-ESM1-5','MPI-ESM1-2-HR','EC-Earth3','INM-CM5-0','MRI-ESM2-0')
ssps <- c('ssp245','ssp585')
prds <- c('2021_2040','2041_2060')
stp <- base::expand.grid(gcms,ssps,prds) %>% base::as.data.frame()
names(stp) <- c('gcm','ssp','prd'); rm(gcms, ssps, prds)
stp <- stp %>%
  dplyr::arrange(gcm,ssp,prd) %>%
  base::as.data.frame()
stp[] <- lapply(stp, as.character) %>% base::as.data.frame()
stp$corrupted <- NA
stp <- tibble::tibble(stp)

# Looping the process
for(i in 1:nrow(stp)){
  
  cat(paste0('>>> Processing: ',stp$gcm[i],', ',stp$ssp[i],', ',stp$prd[i],'\n'))
  
  # Identify index files with issues
  indices <- c('NDWS') # 'NDWL50'
  idx <- indices[1]
  prd <- stp$prd[i] # prd <- '2021_2040'
  ssp <- stp$ssp[i] # ssp <- 'ssp245'
  gcm <- stp$gcm[i] # gcm <- 'ACCESS-ESM1-5'
  ifl <- list.files(path = paste0(root,'/atlas_hazards/cmip6/indices/',ssp,'_',gcm,'_',prd,'/',idx), pattern = '.tif$', full.names = T)
  ifl <- ifl[grep(pattern = 'AVAIL', x = ifl, invert = T)]
  
  # Detect corrupted index files
  future::plan(multicore, workers = parallel::detectCores()-1)
  vrf <- 1:length(ifl) %>%
    furrr::future_map(.f = function(j){
      r <- terra::rast(ifl[j])
      v <- sum(r[!is.na(r)|r != 0])
      return(v)
    }) %>% base::unlist()
  future:::ClusterRegistry("stop")
  gc(reset = T)
  future::plan(sequential)
  
  toCheck <- ifl[vrf == 0]
  if(length(toCheck) > 0){
    date <- basename(toCheck) %>% gsub('NDWS-', '', ., fixed = T) %>% gsub('.tif', '', ., fixed = T)
    
    # Setup daily files per variable
    vrss <- c('pr','tasmax','tasmin')
    stp_dly <- merge(x = data.frame(gcm, ssp, prd),
                     y = base::data.frame(vrss)); rm(vrss)
    names(stp_dly)[ncol(stp_dly)] <- 'var'
    stp_dly[] <- lapply(stp_dly, as.character) %>% base::as.data.frame()
    
    # Corrupted files identification
    crrptd <- 1:nrow(stp_dly) %>%
      purrr::map(.f = function(k){
        
        var <- stp_dly$var[k]
        
        # Define daily file's path
        if(var == 'pr'){
          pth <- paste0(root,'/chirps_cmip6_africa/Prec_',gcm,'_',ssp,'_',prd)
        } else {
          if(var %in% c('tasmax','tasmin')){
            if(var == 'tasmax'){vr <- 'Tmax'}
            if(var == 'tasmin'){vr <- 'Tmin'}
            if(prd == '2021_2040'){yrs <- 2021:2040}
            if(prd == '2041_2060'){yrs <- 2041:2060}
            pth <- paste0(root,'/chirts_cmip6_africa/',vr,'_',gcm,'_',ssp,'_',prd,'/',yrs)
          }
        }
        
        # Filter years with index's problems (only applies for temperature files)
        if(var %in% c('tasmax','tasmin')){
          pth <- pth[grep2(strsplit(date,'-') %>% purrr::map(1) %>% base::unlist() %>% unique(),pth) %>% base::unlist()]
        }
        
        # Corrupted files
        fls2chck <- 1:length(pth) %>%
          purrr::map(.f = function(j){
            fls <- list.files(path = pth[j], pattern = '.tif$', full.names = T) # List daily files
            szs <- file.size(fls) * 1e-6    # File size in MB
            avg <- mean(szs)                # Average files size
            std <- sd(szs)                  # Standard deviation files size
            vrf <- fls[szs < (avg - 5*std)] # Identify files which are 5 std below average (corrupted files)
            return(vrf)
          }) %>% base::unlist() %>% unique()
        
        return(fls2chck) # Return corrupted files per combination
        
      }) %>% base::unlist()
    
    if(length(crrptd) > 0){
      stp$corrupted[i] <- crrptd
    } else {
      stp$corrupted[i] <- NA
    }
    
  } else {
    stp$corrupted[i] <- NA
  }
  
}

# saveRDS(object = stp, file = paste0())
