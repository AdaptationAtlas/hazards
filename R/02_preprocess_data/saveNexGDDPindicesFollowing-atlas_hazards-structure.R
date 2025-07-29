# Generate atlas_hazards structure for Nex-GDDP-CMIP6 indices
# Target directory: ~/common_data/atlas_nex-gddp_hazards/cmip6/indices
# By: H. Achicanoy
# Alliance of Bioversity-CIAT, 2025

# R options
options(warn = -1, scipen = 999)
suppressMessages(library(pacman))
pacman::p_load(purrr)
list.files2 <- Vectorize(FUN = list.files, vectorize.args = 'pattern')

# Setup table
ssps <- c('ssp126','ssp245','ssp370','ssp585')
gcms <- c('ACCESS-CM2','ACCESS-ESM1-5','CanESM5','CMCC-ESM2','EC-Earth3','EC-Earth3-Veg-LR','GFDL-ESM4','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NorESM2-LM','NorESM2-MM','TaiESM1') # c('CanESM5','GFDL-ESM4','INM-CM4-8','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','TaiESM1') # c('ACCESS-CM2','ACCESS-ESM1-5','CMCC-ESM2','EC-Earth3','EC-Earth3-Veg-LR','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0')
prds <- c('2021_2040','2041_2060','2061_2080','2081_2100')
stp_tbl <- expand.grid(ssp = ssps, gcm = gcms, prd = prds, stringsAsFactors = F) |>
  base::as.data.frame() |>
  dplyr::arrange(gcm, ssp, prd)
rm(ssps, gcms, prds)
stp_tbl$folder <- paste0(stp_tbl$ssp,'_',stp_tbl$gcm,'_',stp_tbl$prd)
stp_tbl$ini_year <- strsplit(stp_tbl$prd, '_') |> purrr::map(1) |> unlist() |> as.numeric()
stp_tbl$end_year <- strsplit(stp_tbl$prd, '_') |> purrr::map(2) |> unlist() |> as.numeric()

root <- '/home/jovyan/common_data'

# Available indices
indices <- c('TAVG') # c('NDD','NDWL0','NDWL50','NDWS','PTOT','TAI','TMAX','TMIN','NTX30','NTX35')

for (index in indices) {
  
  cat('.... Copying files for index:',index,'\n')
  
  for (i in 1:nrow(stp_tbl)) {
    
    cat('Over the scenario:',stp_tbl$folder[i],'\n')
    
    # Origin
    org_pth <- paste0(root,'/nex-gddp-cmip6_indices/',stp_tbl$ssp[i],'_',stp_tbl$gcm[i],'/',index)
    fls2copy <- list.files2(path = org_pth, pattern = paste0(index,'-',stp_tbl$ini_year[i]:stp_tbl$end_year[i]), full.names = T) |> as.character()
    
    if (index %in% c('NTX30','NTX35')) {
      # Origin
      trg_index <- gsub('X','x',index)
      # Target
      trg_pth <- paste0(root,'/atlas_nex-gddp_hazards/cmip6/indices/',stp_tbl$ssp[i],'_',stp_tbl$gcm[i],'_',stp_tbl$prd[i],'/',trg_index)
      dir.create(trg_pth, F, T)
      file.copy(from = fls2copy, to = file.path(trg_pth,gsub('X','x',basename(fls2copy))))
    } else {
      # Target
      trg_pth <- paste0(root,'/atlas_nex-gddp_hazards/cmip6/indices/',stp_tbl$ssp[i],'_',stp_tbl$gcm[i],'_',stp_tbl$prd[i],'/',index)
      dir.create(trg_pth, F, T)
      file.copy(from = fls2copy, to = file.path(trg_pth,basename(fls2copy)))
    }
    
  }
  
}
