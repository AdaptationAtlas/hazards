#!/bin/bash

#ACCESS-ESM1-5 MPI-ESM1-2-HR EC-Earth3 INM-CM5-0 MRI-ESM2-0

for GCM in INM-CM5-0
do
    for SCENARIO in ssp245 ssp585
    do
        for PERIOD in 2021_2040 2041_2060
        do
            for YEAR in {1..20}
            do
                echo ----------------------------------------------------------------------
                echo ---- processing ${GCM} - ${SCENARIO} - ${PERIOD} - ${YEAR} -----------
                echo ----------------------------------------------------------------------
                
                Rscript --vanilla ~/Repositories/hazards/R/04_indices/calc_NDWL50_sh.R ${GCM} ${SCENARIO} ${PERIOD} ${YEAR}
            done
        done
    done
done

