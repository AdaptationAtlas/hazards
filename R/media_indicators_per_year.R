## --------------------------------------------------------------------------------- ##
## --------------------------------------------------------------------------------- ##
# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
# .rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, lubridate,foreach,fields,gtools))
root <- "//catalogue/WFP_ClimateRiskPr1/7.Results/"

results <- function(iso,ind, year){ 
  root <- "//catalogue/WFP_ClimateRiskPr1/7.Results/"
  # path <- paste0(root, "Africa/",iso,"/",ind,"_",year)
  path <- paste0(root, "Africa/",iso,"/")
  fls <- list.files(path = path, pattern = ind)
  fls1 <- fls[grep(year, fls)]
  
  fls_f <- paste0(path, fls1) 
  ras <- raster::stack(fls_f)  
  res <- mean(ras, na.rm = T)
  
  outfile <- paste0("//catalogue/WFP_ClimateRiskPr1/8.Final_results_current/country/",iso,"/"); if(!dir.exists(outfile)){dir.create(outfile, F, T)}
  writeRaster(res, paste0(outfile,ind,"_",year,".tif"), overwrite = T)  
  
}

AFR <- raster::shapefile("//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/shps/world_shapefile/Africa/africa.shp")
country <- unique(c(AFR@data$ISO3))
year <- 1995:2014


# Loop through seasons
lapply (1:length(country), function(i){
  cat(paste0("Procesando pais :::: ",country[i], "\n" ))
  iso <- country[i]
  lapply(1: length(year), function(j){
    cat(paste0("Procesando year :::: ",year[j], "\n" ))
    results(year = year[j], ind = "HSH", iso = iso)
  })
  
})
