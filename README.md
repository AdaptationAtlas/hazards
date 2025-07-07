# Adaptation atlas hazards data generation scripts

Scripts for processing of historical and future hazards data. For more information see the [wiki](https://github.com/AdaptationAtlas/hazards/wiki).

**Updates: July, 2025**

Currently, we have hazard data from two sources. The first dataset uses the bias correction delta method and has a spatial resolution of 0.05° (~5 km) at the equator (```main``` GitHub branch). For this dataset, we used five GCMs across the four SSPs (SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5) and the 2021-2100 period:

```ACCESS-ESM1-5, EC-Earth3, INM-CM5-0, MPI-ESM1-2-HR, MRI-ESM2-0```

The second dataset corresponds to [Nex-GDDP-CMIP6](https://developers.google.com/earth-engine/datasets/catalog/NASA_GDDP-CMIP6). This dataset has a spatial resolution of 0.25° (```nexgddp``` GitHub branch). Eighteen GCMs were used across the four SSPs and for the 2021-2100 period:

```ACCESS-CM2, ACCESS-ESM1-5, CMCC-ESM2, CanESM5, EC-Earth3, EC-Earth3-Veg-LR, GFDL-ESM4, INM-CM4-8, INM-CM5-0, IPSL-CM6A-LR, KACE-1-0-G, MIROC6, MPI-ESM1-2-HR, MPI-ESM1-2-LR, MRI-ESM2-0, NorESM2-LM, NorESM2-MM, TaiESM1```

Both datasets cover the five required daily variables (minimum and maximum temperatures, precipitation, solar radiation, and relative humidity) under the same initial simulation conditions (r1i1p1f1).

## Hazards included

As of January 2022, we have calculated a total of 10 hazard indices. The descriptions of these are available in the wiki, under [hazard definitions](https://github.com/AdaptationAtlas/hazards/wiki/Hazards-definitions). The indices are as follows,

1. **Drought stress**
  - NDD: number of dry days
  - NDWS: number of soil moisture stress days
  - TAI: Thornthwaite's aridity index

2. **Heat stress**
  - NTx35: number of heat stress (Tmax \> 35ºC) days for crops
  - NTx40: number of extreme heat stress (Tmax \> 40ºC) days for crops
  - HSM_NTx35: number of heat stress days for maize during the crop growing season
  - HSH: human heat stress index
  - THI: cattle thermal humidity index

3. **Waterlogging and flooding**
  - NDWL0: number of days with soil waterlogging at moisture content at start of saturation or above.
  - NDWL50: number of days with soil waterlogging at moisture content 50% between field capacity and saturation.

These data are available for four SSPs (SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5), four periods (2021-2040, 2041-2060, 2061-2080, 2081-2100), and five GCMs (ACCESS-ESM1-5, MPI-ESM1-2-HR, EC-Earth3, INM-CM5-0, MRI-ESM2-0) in addition to the multi-model ensemble mean.

## Where is the data?
While the Atlas is under development, the data is internally available to Atlas scientists in the EiA server under `~/common_data/atlas_hazards/cmip6`. The data are also available in the [AWS S3](s3://digital-atlas/Updates_for_MVP_Release/1_hazards/) and [Google Cloud](gs://adaptation-atlas/cmip6_hazards/) buckets.

## What code is available?
This repository is divided into 7 basic processes, each available under specific folders containing [the code](https://github.com/AdaptationAtlas/hazards/tree/main/R) for processing all Atlas hazard layers. Each folder contains various functions designed to work as generically as possible so the workflow can be replicated. The below provides an overview of each process, and some of the functions therein.

To get access to the code it is best to first clone the repository locally, as most of the scripts require small modifications (e.g., local folders). To clone the repository, use

```bash
git clone https://github.com/AdaptationAtlas/hazards.git
```

### Data downloads
For hazards calculations, right now we require five datasets, namely, CHIRPS/CHIRTS (from the Climate Hazards Group), AgERA5 (from the EU Copernicus service), CMIP6 projections, and crop calendars. For each dataset, we have built a script / function that helps download the dataset based on specified dates and/or spatial domain. For instance, for CHIRPS, the following code will download a global CHIRPS data layer (in GeoTiff format)

```r
getChirps(date = as.Date("1981-01-01""))
```

An interesting function we built is the `downloadCMIP6` function, which leverages work done by the [SantanderMetGroup](https://github.com/SantanderMetGroup/ATLAS) on the CMIP6 Atlas. This function allows downloading raw CMIP6 data in an automated way, without having to pass through the ESFG CMIP6 portal.

```r
cmip6_data <- downloadCMIP6(ds_name    = "CMIP6_ACCESS-ESM1-5_scenario_r1i1p1f1", 
                            rcp        = "ssp246", 
                            varname    = "tasmax", 
                            years.hist = 1995:2014, 
                            years.rcp  = 2021:2060, 
                            lons       = c(-23, 59), 
                            lats       = c(-37, 40), 
                            basedir    = my_local_dir)
```

### Data pre-processing
Contains functions to create various masks for Africa that are needed for the post-processing of the hazard layers. There is also a function to pre-process crop calendar data and make it compatible (resolution and extent), with the rest of the data layers.

### Bias-correction
The two main scripts here are `bc_calc_anomalies.R` and `getDailyFutureData.R`. These scripts help first calculating monthly climate anomalies, and then apply these to the historical climate data. An example is provided here,

```r
#historical climatology
his_clm <- calc_climatology(data_file = "CMIP6_ACCESS-ESM1-5_historical_r1i1p1f1__Africa_daily.tif", 
                            period   = 1995:2014, 
                            sce_lab  = "historical",
                            gcm_name = "CMIP6_ACCESS-ESM1-5",
                            varname  = "tasmax",
                            mth_dir  = "my_folder_with_monthly_annual_files",
                            clm_dir  = "my_folder_with_climatology_files")
      
#future climatology
rcp_clm <- calc_climatology(data_file = "CMIP6_ACCESS-ESM1-5_ssp245_r1i1p1f1__Africa_daily.tif", 
                            period    = 2021:2040, 
                            sce_lab   = "ssp245",
                            gcm_name  = "CMIP6_ACCESS-ESM1-5",
                            varname   = "tasmax",
                            mth_dir   = "my_folder_with_monthly_annual_files",
                            clm_dir   = "my_folder_with_climatology_files")
      
#interpolate anomalies
rcp_anom <- intp_anomalies(his_clm  = "my_folder_with_climatology_files_hist", 
                           rcp_clm  = "my_folder_with_climatology_files_rcp", 
                           anom_dir = "my_output_folder_for_anomalies", 
                           ref      = reference_raster_object,
                           gcm_name = "CMIP6_ACCESS-ESM1-5",, 
                           rcp      = "ssp245", 
                           varname  = "tasmax", 
                           period   = 2021:2040)
```

Now produce the daily future data:

```r
get_daily_future_data(gcm = "ACCESS-ESM1-5",
                      ssp = "ssp245",
                      var = "tasmax",
                      prd = "2021_2040")
```

**Updates: July, 2025**

The daily future precipitation files were corrected to preserve a reasonable range of variation. Due to mathematical indeterminacies in the delta formula [Navarro-Racines et al., (2020)](https://www.nature.com/articles/s41597-019-0343-8), extreme delta values were constrained from -1.5 to 1.5; this means maximum decreases and increases of 150%. Additionally, the minimum monthly precipitation values used to calculate anomalies were capped at 1mm/month in both the historical and future GCM datasets.

### Hazard indices
This set of functions calculates each index in spatial form, taking as input the downscaled / bias-corrected CMIP6 daily data. A script is provided for each hazard index. Note that these scripts contain both paths, execution loops (or `purrr::map()`), and functions. So some modifications will be in order for them to work in a new computing environment than ours. Most functions take the year and month as input, and in some cases other additional parameters. For example, for NTx40, we built the function `calc_ntx()`, which is used as follows:

```r
gcm <- "ACCESS-ESM1-5"
ssp <- "ssp245"
prd <- "2021_2040"
calc_ntx(yr=2021, mn="01", thr=40)
```

### Final maps, and meta-data
Each of the indices has a meta-data function, which uses the Atlas [metadata](https://github.com/AdaptationAtlas/metadata) functions. This function will calculate long-term statistics (using the function `continuous_map()` in script `calc_LongTermStats.R`). It will also generate discrete maps of hazard severity levels (function `category_map()` in script `calc_discreteMaps.R`), and create the actual meta-data and bucket uploadable data files. These scripts need to be sourced (i.e., they are not functions).

### Bucket uploads
These scripts will automatically upload the processed output into the AWS S3 and Google Buckets. They use the Linux system library `gsutil` and the R package `aws.s3`.


