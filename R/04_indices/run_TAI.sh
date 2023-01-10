#!/bin/bash

for GCM in ACCESS-ESM1-5 MPI-ESM1-2-HR EC-Earth3 INM-CM5-0 MRI-ESM2-0
do
    for SCENARIO in ssp245 ssp585
    do
        for PERIOD in 2021_2040 2041_2060
        do
            Rscript --vanilla ~/Repositories/hazards/R/04_indices/calc_TAI_sh.R ${GCM} ${SCENARIO} ${PERIOD} 1 4
            Rscript --vanilla ~/Repositories/hazards/R/04_indices/calc_TAI_sh.R ${GCM} ${SCENARIO} ${PERIOD} 5 8
            Rscript --vanilla ~/Repositories/hazards/R/04_indices/calc_TAI_sh.R ${GCM} ${SCENARIO} ${PERIOD} 9 12
            Rscript --vanilla ~/Repositories/hazards/R/04_indices/calc_TAI_sh.R ${GCM} ${SCENARIO} ${PERIOD} 13 17
            Rscript --vanilla ~/Repositories/hazards/R/04_indices/calc_TAI_sh.R ${GCM} ${SCENARIO} ${PERIOD} 18 20
        done
    done
done

