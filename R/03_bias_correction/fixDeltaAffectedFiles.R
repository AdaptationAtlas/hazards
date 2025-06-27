# Define parameters
stp <- data.frame(var = 'pr',
                  ssp = c('ssp245',rep('ssp585',5),'ssp370'),
                  gcm = 'EC-Earth3',
                  prd = c('2021_2040',rep('2081_2100',3),rep('2061_2080',2),'2081_2100'),
                  mn  = c('02','08','09','10','08','09','09'))

# var = 'pr'; ssp = 'ssp245'; gcm = 'EC-Earth3'; prd = '2021_2040'; mn = '02'
# var = 'pr'; ssp = 'ssp585'; gcm = 'EC-Earth3'; prd = '2081_2100'; mn = '08'
# var = 'pr'; ssp = 'ssp585'; gcm = 'EC-Earth3'; prd = '2081_2100'; mn = '09'
# var = 'pr'; ssp = 'ssp585'; gcm = 'EC-Earth3'; prd = '2081_2100'; mn = '10'
# var = 'pr'; ssp = 'ssp585'; gcm = 'EC-Earth3'; prd = '2061_2080'; mn = '08'
# var = 'pr'; ssp = 'ssp585'; gcm = 'EC-Earth3'; prd = '2061_2080'; mn = '09'
var = 'pr'; ssp = 'ssp370'; gcm = 'EC-Earth3'; prd = '2081_2100'; mn = '09'
afr <- terra::vect('~/common_data/atlas_hazards/roi/africa.gpkg')

# Load anomalies interpolation function (modified)
source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/bc_calc_anomalies_fixed.R'); gc()

# Load get daily future data function (modified)
source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/getDailyFutureData_fixed.R'); gc()

# Load the monthly total precipitation function (modified)
source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/calc_PTOT_fixed.R'); gc()

px_cnt <- utils::read.csv('C:/Users/haachicanoy/Downloads/pixel_counts_thr_15.csv')
px_cnt$month <- base::sprintf('%02.0f', px_cnt$month)
px_cnt$id <- paste0(px_cnt$month,'--',px_cnt$folder)
px_tfx <- px_cnt[px_cnt$pixel_value != 0,]
px_tfx <- px_tfx |> dplyr::arrange(-count)
rownames(px_tfx) <- 1:nrow(px_tfx)
px_tfx[,c('month','ssp','gcm','prd','folder')] |> unique() |> dim()

dfm <- utils::read.csv('~/common_data/affected_geographies/pr_deltas_summaries.csv')
cond1 <- which(dfm$min < -2)
cond2 <- which(dfm$max > 2)

mtch <- sort(unique(union(cond1, cond2))); rm(cond1, cond2)

dfm_ok <- dfm[setdiff(1:nrow(dfm),mtch),]
sort(table(dfm_ok$file), decreasing = T)

# Code to fix all precipitation files
ssps <- c('ssp126','ssp245','ssp370','ssp585')
gcms <- c('ACCESS-ESM1-5','EC-Earth3','INM-CM5-0','MPI-ESM1-2-HR','MRI-ESM2-0')
prds <- c('2021_2040','2041_2060','2061_2080','2081_2100')
mns  <- sprintf('%02.0f',1:12)
stp_tbl <- expand.grid(ssp = ssps, gcm = gcms, prd = prds, mn = mns, stringsAsFactors = F) |>
  base::as.data.frame() |>
  dplyr::arrange(gcm, ssp, prd)

for (i in 1:nrow(stp_tbl)) {
  
  lgs_file <- '~/common_data/affected_geographies/pr_files_fixed.csv'
  if (!file.exists(lgs_file)) {
    
    ssp <- stp_tbl$ssp[i]
    gcm <- stp_tbl$gcm[i]
    prd <- stp_tbl$prd[i]
    mn  <- stp_tbl$mn[i]
    var <- 'pr'
    
    cat(paste0('>>> Fix anomalies for: ',ssp,', ',gcm,', and, ',prd,'\n'))
    # Load anomalies interpolation function (modified)
    source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/bc_calc_anomalies_fixed.R'); gc()
    cat(paste0('... Fix daily rasters for month: ',mn,'\n'))
    # Load get daily future data function (modified)
    source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/getDailyFutureData_fixed.R'); gc()
    cat(paste0('... Fix monthly rasters for month: ',mn,'\n'))
    # Load the monthly total precipitation function (modified)
    source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/calc_PTOT_fixed.R'); gc()
    
    logs <- stp_tbl[i,]
    utils::write.csv(logs, file = lgs_file, row.names = F)
  } else {
    logs <- utils::read.csv(lgs_file)
    
    i <- (nrow(logs) + 1)
    ssp <- stp_tbl$ssp[i]
    gcm <- stp_tbl$gcm[i]
    prd <- stp_tbl$prd[i]
    mn  <- stp_tbl$mn[i]
    var <- 'pr'
    
    cat(paste0('>>> Fix anomalies for: ',ssp,', ',gcm,', and, ',prd,'\n'))
    # Load anomalies interpolation function (modified)
    source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/bc_calc_anomalies_fixed.R'); gc()
    cat(paste0('... Fix daily rasters for month: ',mn,'\n'))
    # Load get daily future data function (modified)
    source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/getDailyFutureData_fixed.R'); gc()
    cat(paste0('... Fix monthly rasters for month: ',mn,'\n'))
    # Load the monthly total precipitation function (modified)
    source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/calc_PTOT_fixed.R'); gc()
    
    logs <- rbind(logs,
                  stp_tbl[i,])
    utils::write.csv(logs, file = lgs_file, row.names = F)
  }
  cat('\n\n')
}
