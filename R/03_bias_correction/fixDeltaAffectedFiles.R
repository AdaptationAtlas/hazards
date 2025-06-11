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
source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/bc_calc_anomalies_fixed.R')

# Load get daily future data function (modified)
source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/getDailyFutureData_fixed.R')

# Load the monthly total precipitation function (modified)
source('https://raw.githubusercontent.com/AdaptationAtlas/hazards/refs/heads/main/R/03_bias_correction/calc_PTOT_fixed.R')
