# Free up space deleting Nex-GDDP-CMIP6 raw files
# By: H. Achicanoy
# Alliance Bioversity-International & CIAT, 2025

# The purpose of this script is to free up space in our
# storage system so that we can process the baselines.
# However, it has one disadvantage: if any pre-processed
# file became corrupted, we must download the raw file
# to fix it again. BE CAREFUL WITH ITS USE.

vrs <- c('pr','tasmax','tasmin','hurs','rsds')

for (vr in vrs) {
  
  pth <- file.path('~/common_data/nex-gddp-cmip6_raw',vr)
  drs <- list.dirs(path = pth, recursive = F)
  drs <- drs[grep('ssp', drs)]
  
  for (dr in drs) {
    system(paste0('rm -r ',dr))
  }
  
}
