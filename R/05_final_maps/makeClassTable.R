#Calculate discrete maps
#HA/JRV, Dec 2022

#function to create the categories
make_class_tb <- function() {
  #NDD: number of dry days per month
  cat_df <- c(-Inf,  15,  1,
              15,  20,  2,
              20,  25,  3,
              25, Inf,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  full_tb <- cat_df %>%
    dplyr::mutate(index_name="NDD", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  
  #NDWS: number of water stress days per month
  cat_df <- c(-Inf,  15,  1,
              15,  20,  2,
              20,  25,  3,
              25, Inf,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="NDWS", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #NTx40: number of heat stress for crops using a generic threshold of 40 degC (Tmax >= 40)
  cat_df <- c(-Inf,   1,  1,
              1,   5,  2,
              5,  10,  3,
              10, Inf,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="NTx40", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #NDWL50: number of waterlogging days per month at 50% saturation
  cat_df <- c(-Inf,   2,  1,
              2,   5,  2,
              5,   8,  3,
              8, Inf,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="NDWL50", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #NDWL0: number of waterlogging days per month at 0% saturation
  cat_df <- c(-Inf,   2,  1,
              2,   5,  2,
              5,   8,  3,
              8, Inf,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="NDWL0", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #TAI: Thornthwaite's aridity index
  cat_df <- c(-Inf,  40,  1,
              40,  60,  2,
              60,  80,  3,
              80, Inf,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="TAI", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #HSH: human heat stress index as calculated by Philip Thornton (https://en.wikipedia.org/wiki/Heat_index)
  cat_df <- c(-Inf,  27,  1,
              27,  32,  2,
              32,  41,  3,
              41, Inf,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="HSH", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #HSH_TW: human heat stress index equivalent to wet bulb temperature (Stull, 2011)
  cat_df <- c(-Inf,  25,  1,
              25,  30,  2,
              30,  35,  3,
              35, Inf,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="HSH_TW", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #THI: thermal humidity index for livestock
  cat_df <- c(-Inf,   72,  1,
              72,   78,  2,
              78,  89,  3,
              89, Inf,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="THI", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #THIMS: THI multiple species expressed as probability (frequency) of occurrence of extreme heat
  cat_df <- c(0.00, 0.25,  1,
              0.25, 0.50,  2,
              0.50, 0.75,  3,
              0.75, 1.00,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="THIMS", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #HSM_NTx35: number of heat stress days during the growing season of maize (Tmax >= 35)
  cat_df <- c(-Inf,   7,  1,
              7,   14,  2,
              14,  21,  3,
              21, Inf,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="HSM_NTx35", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #NTx35: number of heat stress days in the entire year in maize areas (Tmax >= 35)
  cat_df <- c(-Inf,   7,  1,
              7,   14,  2,
              14,  21,  3,
              21, Inf,  4)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="NTx35", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #PTOT: total monthly precipitation
  cat_df <- c(-Inf,   50,  4,
              50,   100,  3,
              100,  150,  2,
              150, Inf,  1)
  cat_df <- as.data.frame(matrix(cat_df, ncol = 3, byrow = T))
  names(cat_df) <- c("lower_lim", "upper_lim", "class")
  cat_df <- cat_df %>%
    dplyr::mutate(index_name="PTOT", version="mvp") %>%
    dplyr::select(index_name, version, lower_lim:class)
  full_tb <- full_tb %>%
    dplyr::bind_rows(cat_df)
  
  #add class names
  cls <- data.frame(class = 1:4, description = c('No significant stress','Moderate','Severe','Extreme'))
  full_tb <- full_tb %>%
    dplyr::full_join(., cls, by="class")
  
  #return object
  return(full_tb)
}
