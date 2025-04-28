screen -dmS ndwl50_earth Rscript fast_calc_NDWL50.R 'ssp370' 'EC-Earth3' '2041_2060'
screen -dmS ndwl50_mpi Rscript fast_calc_NDWL50.R 'ssp370' 'MPI-ESM1-2-HR' '2021_2040'
screen -dmS ndwl50_mpi2 Rscript fast_calc_NDWL50.R 'ssp370' 'MPI-ESM1-2-HR' '2041_2060'


