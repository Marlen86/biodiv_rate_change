# ------------------------------------------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------------------------------------------
require(tidyverse)
require(reshape2)
require(broom)
require(terra) 
require(sf)
# ------------------------------------------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------------------------------------------
# To process the climatic variables to get the mean for each year
var_calc <- function(directory1, directory2, select_files2, period1, period2, season, var_name){
  # 1938-1980 => select rasters for this period
  setwd(directory1)
  files1 <- list.files( , pattern = "bil.bil")[c(TRUE, FALSE)]
  # Keeping selected months
  files1      <- files1[select_files2]
  files1_vars <- files1[period1] 
  st_80_vars  <- rast(files1_vars)
  # 1981-2019 => select rasters for this period
  setwd(directory2)
  files2 <- list.files( , pattern="bil.bil")[c(TRUE, FALSE)]
  # Keeping selected months
  files2       <- files2[select_files2]
  files2_vars  <- files2[period2]
  st_8119_vars <- rast(files2_vars) 
  # Join rasters for both periods
  st_vars      <- c(st_80_vars, st_8119_vars)
  # Change resolution to 10km (10000m) with aggregate and function mean (default function)
  aggr_raster_vars <- aggregate(st_vars, 2, na.rm = T)
  # Create a vector of indexes to use in the function tapp
  indices_vars <- sort(rep(1:82, 3)) # 1:82 -> number of years, 3 -> number of layers per year
  # Calculating the mean value of the variable per season per year from 1938 to 2019
  vars_all         <- tapp(aggr_raster_vars, indices_vars, fun = mean) 
  # indices_vars vector indicates to perform the function using only the 3 layers with the same number
  names(vars_all)  <- paste(var_name, season, c(1938:2019), sep = '_')
  # Change the projection of climatic layers to the projection of the LULC layers 
  # (so all environmental layers can be stacked together and worked with)
  # LULC projection crs: 
  # +proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs 
  vars_all1 <- project(vars_all, template) # bi-linear interpolation (default) => appropriate for continuous variables
  setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/Climatic_vars_1938_1979_2009_2018_10km/may2023")
  writeRaster(vars_all1, paste(var_name, "_1938_2019_10km_july2024.tif", sep = ''), overwrite = TRUE)
  rm(files1, st_80_vars, st_8119_vars, st_vars, indices_vars, vars_all, aggr_raster_vars, vars_all1)
} 

# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
# LULC data
# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
# Historical Backcasting
setwd("/Volumes/G2T/memory/Thesis/CONUS_LandUse_1938_2100_250m/ModeledHistoric")
st_3891 <- rast(list.files( , pattern = 'tif')[1:54])
# Historical LULC change in the conterminous United States from 1992 to 2005. 
setwd("/Volumes/G2T/memory/Thesis/CONUS_LandUse_1938_2100_250m/CONUS_Landcover_Historical")
st_9205 <- rast(list.files( , pattern = 'tif')[1:14])
# Projected LULC change in the conterminous United States from 2006 to 2018. 
# Annual LULC maps were produced at 250-m resolution, with 14 LULC classes.
# For description of the four scenarios available for CONUS, read the 
# pdf "sres_Future_Scenarios.pdf" available online.
# A1B scenario
setwd("/Volumes/G2T/memory/Thesis/CONUS_LandUse_1938_2100_250m/CONUS_Landcover_A1B")
st_0619_A1B <- rast(list.files( , pattern = 'tif')[1:14])
# A2 scenario
setwd("/Volumes/G2T/memory/Thesis/CONUS_LandUse_1938_2100_250m/CONUS_Landcover_A2")
st_0619_A2  <- rast(list.files( , pattern = 'tif')[1:14])
# B1 scenario
setwd("/Volumes/G2T/memory/Thesis/CONUS_LandUse_1938_2100_250m/CONUS_Landcover_B1")
st_0619_B1  <- rast(list.files( , pattern = 'tif')[1:14])
# B2 scenario
setwd("/Volumes/G2T/memory/Thesis/CONUS_LandUse_1938_2100_250m/CONUS_Landcover_B2")
st_0619_B2  <- rast(list.files( , pattern = 'tif')[1:14])
# Get all rasters together
st_3819     <- rast(list(st_3891, st_9205, st_0619_A1B, st_0619_A2,st_0619_B1, st_0619_B2)) 
# Clean environment a little
rm(st_0619_A1B,st_0619_A2, st_0619_B1, st_0619_B2, st_3891,st_9205)
# ------------------------------------------------------------------------------------------------------------------
# Layerizing each year's map into one map for each land cover type for each year
# ------------------------------------------------------------------------------------------------------------------
# Loop for layerizing each year's map into one map for each land cover type for 
# each year and combining layers into six main LULC types
# Layers:
# 1: 0, 2: water, 3: developed, 4: mining,  5: barren, 6: deciduous forest, 
# 7: evergreen forest; 8: mixed forest, 9: grassland, 10: shrubland, 11 cropland, 
# 12: hay/pasture, 13: herbaceous wetland, 14: woody wetland, 15: ice
indices <- c(1,1,2,1,7,3,3,3,4,4,5,5,6,6,1)
# Line above was used to combine individual layers into broader LULC groups (Others, 
# Barren, Urban, Forests, Grassland_Shrubland, Crops_Pasture, Wetlands)
for (i in 1 : nlyr(st_3819)) {
  # Get year and LULC scenario (if applicable)
  if(strsplit(names(st_3819[[i]]),"_")[[1]][2]%in%c('Backcasting','Historical')){
    year <- strsplit(names(st_3819[[i]])[1],"_")[[1]][3] # For years 1938-2005
  }else{
    year <- paste(strsplit(names(st_3819[[i]]),"_")[[1]][3], strsplit(names(st_3819[[i]]),"_")[[1]][2], sep = '_') # For years 2006-2019
  }
  print(year)
  # Convert each land cover into an individual layer
  layers <- segregate(st_3819 [[i]])
  # Change resolution to 10km (10000m) with aggregate and function mean
  aggr_raster <- aggregate(layers, 40, na.rm = T)
  rm(layers)
  # Select the raster to mask
  raster_to_maks <- aggr_raster[[1]]
  raster_to_maks[raster_to_maks > 0] <- NA
  # Mask all the rasters
  aggr_raster <- rast(lapply(aggr_raster, function(x) {x[is.na(raster_to_maks[])] <- NA; return(x)}))
  if(length(names(aggr_raster)) != 15){
    aggr_raster <- aggr_raster[[c(1:3,7:18)]]
  }
  # Write rasters individually
  setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/CONUS_LandUse_1938_1979_2009_2018_10km/per_year/Individual_lulc")
  writeRaster(aggr_raster,
              paste(year,"_CONUS_LandUse_10km.tif", sep = ''),
              filetype = 'GTiff',
              overwrite = TRUE)
  # To combine individual layers into broader LULC groups 
  # (others, Barren, Urban, Forests ,Grassland_Shrubland, Crops_Pasture, Wetlands)
  setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/CONUS_LandUse_1938_1979_2009_2018_10km/per_year/Grouped_lulc")
  datasum <- tapp(aggr_raster, indices, fun = sum)
  names(datasum) <- c(paste(year,'_Others', sep = ''),
                      paste(year,'_Urban_Developed', sep=''),
                      paste(year,'_Barren', sep=''),
                      paste(year,'_Forests', sep=''),
                      paste(year,'_Grassland_Shrubland', sep=''), 
                      paste(year,'_Crops_Pasture', sep=''),
                      paste(year,'_Wetlands', sep=''))
  writeRaster(datasum,
              paste(year,"_CONUS_LandUse_Grouped_10km.tif", sep = ''),
              filetype = 'GTiff', 
              names = names(datasum),
              overwrite = TRUE)
  gc()
}
# ------------------------------------------------------------------------------------------------------------------
# Making raster bricks of the same land use each for each period
# ------------------------------------------------------------------------------------------------------------------
# For the combined layers
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/CONUS_LandUse_1938_1979_2009_2018_10km/per_year/Grouped_lulc")
bricks_3819 <- lapply(list.files(), rast) 
dt <- c("y1990_2019_A1B_","y1990_2019_A2_","y1990_2019_B1_","y1990_2019_B2_")
files_j <- list(c(53:68, seq(69,124,4)), # A1B
                c(53:68, seq(70,124,4)), # A2
                c(53:68, seq(71,124,4)), # B1
                c(53:68, seq(72,124,4))) # B2

for(j in 1: length(dt)){ # Select files for different periods
  for(i in 1:length(names(bricks_3819 [[1]]))) { # Combine layers of same LULC for the years in each period
    same_LULC_type <- lapply(bricks_3819[files_j[[j]]], function (x) {layer <- x[[i]]; return(layer)})
    same_LULC_type <- rast(same_LULC_type)
    name_i <- gsub("y1990_","" ,names(same_LULC_type)[1])
    print(name_i)
    setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/CONUS_LandUse_1938_1979_2009_2018_10km/may2023/per_lulc_1990_2019")
    writeRaster(same_LULC_type,
                paste(dt[j],name_i,"_CONUS_LandUse_Grouped_10km.tif", sep=''),
                filetype = 'GTiff',
                overwrite = TRUE)
  }
}
#--------------------------------------------------------------------------------------------------------------------------------------
# Extracting the mean value of lulc variable for each route using mean function
#--------------------------------------------------------------------------------------------------------------------------------------
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/CONUS_LandUse_1938_1979_2009_2018_10km/may2023/per_lulc_1990_2019")
template <- rast(list.files( , pattern = '.tif')[1])[[1]]
files_j  <- list.files( , pattern = '.tif')
# Template for North America
setwd("~/Desktop")
world <- read_sf(dsn = ".", layer = "countriesWGS84")
america <- world [c(252), ] 
america <- st_transform(america, crs = crs(template))
ext <- ext(-2357953, 2292047, 235292.6, 3175293)
america <- st_crop(america, ext)
# BBS routes layer
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity")
# Sampling units
data        <- read.csv('taxonomic_bdiv_ch_1990_2019_nomore1missing_data_July2023.csv', h = T)
sampl_units <- data[ ,1:3][!duplicated(data$CountryStateRoute),]
# Transforming the matrix into spatial points data frame, so it can be used to extract env. variable values
coord_units <- st_as_sf(sampl_units, 
                        coords = c("Longitude","Latitude"), 
                        crs = crs("+proj=longlat +datum=WGS84 +no_defs"))
# Transforming the coordinates of my point to match the environmental layers
transf_coord_units <-  st_transform(coord_units, crs(template))
rm(world,ext,template,data,coord_units)
# Set the buffer resolutions to be used in the loops below
resol     <- c("12_5","25","50")
buff_size <-list(st_buffer(transf_coord_units,12500),
                 st_buffer(transf_coord_units,25000),
                 st_buffer(transf_coord_units,50000))
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/CONUS_LandUse_1938_1979_2009_2018_10km/may2023/per_lulc_1990_2019")
for (j in 1 : length(files_j)){ # Calculations by layers
  name_layers <- gsub("_CONUS_LandUse_Grouped_10km", "",         # Obtain the name without the _CONUS_LandUse_Grouped_10km
                      gsub("y1990_2019",  "",                    # Obtain the name without the y1990_2019
                           strsplit(files_j[j],'.tif')[[1]][1])) # Obtain the name without the .tif
  print(name_layers)
  layer_j <- rast(files_j[j])
  layer_j <- mask (layer_j, america)
  names(layer_j) <- paste('y',c(1990:2019), name_layers, sep = '')
  env_dt <- data.frame()
  for (i in 1: length(resol)){
    resol_i     <- resol[i]
    buff_size_i <- buff_size[[i]]
    env_var_dtf <-  data.frame(extract(layer_j,
                                       buff_size_i, 
                                       method = 'simple', 
                                       fun    = mean, 
                                       na.rm  = TRUE,
                                       bind   = T,
                                       xy     = T))  
  
    env_var_dtf$resoln <- resol_i
    env_var_dtf <- env_var_dtf[ ,c(1,32,2:31)]
    env_dt <- rbind(env_dt, env_var_dtf)
  }
  sampl_units <- full_join(sampl_units, env_dt)
}
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity")
write.csv(sampl_units, 'lulc_vars_mean_before_temp_trend.csv', row.names = F)
# ------------------------------------------------------------------------------------------------------------------
# Univariate rate of change of LULC
# ------------------------------------------------------------------------------------------------------------------
# Temporal trends
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity")
lulc_vars <- read.csv('lulc_vars_mean_before_temp_trend.csv')
long <- lulc_vars%>% 
  pivot_longer(
    cols = c("y1990_A1B_Barren":"y2019_B2_Wetlands"), 
    names_to = "Year",
    values_to = "value"
  )
long$lulc_var <- sub("^[^_]*_", "", long$Year)
long$Year<- as.integer(gsub("y", "", gsub("\\_.*","",long$Year)))
temp_trend <- data.frame(long %>%
                           group_by( Longitude, Latitude,CountryStateRoute, resoln, lulc_var) %>%
                           nest() %>%
                           mutate(
                             fit = map(data, ~ lm(value ~ Year, data = .x)),
                             tidied = map(fit, tidy)
                           ) %>%
                           unnest(tidied, .drop = TRUE))[ ,-c(6,7)] %>%
  filter(term =='Year')
wide <- arrange(dcast(temp_trend[ ,-6], Longitude + Latitude + CountryStateRoute ~ lulc_var+resoln, value.var = 'estimate'), CountryStateRoute)
write.csv(wide,'lulc_vars_lm_tempTrend_1990_2019.csv', row.names = F)
# ------------------------------------------------------------------------------------------------------------------
# Extracting percent of LULC per ecoregion
# ------------------------------------------------------------------------------------------------------------------
# LULC data
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/CONUS_LandUse_1938_1979_2009_2018_10km/may2023/per_lulc_1990_2019")
template <- rast(list.files( , pattern = '.tif')[1])[[1]]
files_j  <- list.files( , pattern = '.tif')
# Template for North America
setwd("~/Desktop")
world <- read_sf(dsn = ".", layer = "countriesWGS84")
america <- world [c(252), ] 
america <- st_transform(america, crs = crs(template))
ext <- ext(-2357953, 2292047, 235292.6, 3175293)
america <- st_crop(america, ext)
# Ecoregions
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/NA_Terrestrial_Ecoregions_v2_Level_I_Shapefile/NA_TerrestrialEcoregions_LI/data")
ecoregions <- read_sf(dsn = ".", layer = "NA_Terrestrial_Ecoregions_v2_level1")
ecoregions <- st_transform(ecoregions, crs = st_crs(template))
# Extracting the mean value of lulc variable for each route using mean function
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/CONUS_LandUse_1938_1979_2009_2018_10km/may2023/per_lulc_1990_2019")
perc_dt <- c("Eastern Temperate Forests", "Great Plains", "North American Deserts",
                                    "Northern Forests", "Northwestern Forested Mountains")
all_data <- data.frame()
for (i in 1:length(perc_dt)){
  ecor_i <- ecoregions%>% filter(NameL1_En == perc_dt[i])
  ecor_i_sv <- vect(ecor_i)
for (j in 1 : length(files_j)){ # Calculations by layers
  name_layers <- gsub("_CONUS_LandUse_Grouped_10km", "",         # Obtain the name without the _CONUS_LandUse_Grouped_10km
                      gsub("y1990_2019",  "",                    # Obtain the name without the y1990_2019
                           strsplit(files_j[j],'.tif')[[1]][1])) # Obtain the name without the .tif
  print(name_layers)
  layer_j <- rast(files_j[j])
  layer_j <- mask (layer_j, america)
  names(layer_j) <- paste('y',c(1990:2019), name_layers, sep = '')
  # Extract values from the raster for each county #
  lulc_prop_dtf <-  data.frame(terra::extract(layer_j ,
                                              ecor_i_sv , 
                                              method = 'simple', 
                                              na.rm  = TRUE,
                                              bind   = T,
                                              xy     = T))  
 summ_lulc <- lulc_prop_dtf %>% 
    summarise(across(where(is.numeric), 
                     list(mean = mean), 
                     na.rm = T))
  summ_lulc$mean_perc <- rowMeans(summ_lulc[,-c(1,32,33) ])
  summ_lulc$ lulc <- name_layers
  summ_lulc$ NameL1_En <- perc_dt[i]
  all_data <- rbind(all_data , summ_lulc[ ,c(34:36)])
}
  }

setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity")
write.csv(all_data, 'lulc_vars_mean_prop_ecoregions_1990_2019.csv', row.names = F)

# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
# Climatic Variables
# ------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------
# Defining the variables of the function to loop through
directories1 <- c("/Volumes/G2T/memory/Thesis/PRISM_1895_2020/PRISM_tmax_stable_4kmM3_189501_198012_bil",
                  "/Volumes/G2T/memory/Thesis/PRISM_1895_2020/PRISM_tmin_stable_4kmM3_189501_198012_bil",
                  "/Volumes/G2T/memory/Thesis/PRISM_1895_2020/PRISM_ppt_stable_4kmM2_189501_198012_bil")
directories2 <- c("/Volumes/G2T/memory/Thesis/PRISM_1895_2020/PRISM_tmax_stable_4kmM3_198101_202008_bil",
                  "/Volumes/G2T/memory/Thesis/PRISM_1895_2020/PRISM_tmin_stable_4kmM3_198101_202008_bil",
                  "/Volumes/G2T/memory/Thesis/PRISM_1895_2020/PRISM_ppt_stable_4kmM3_198101_202008_bil")
select_files2s <- list(c(rep(FALSE,4), rep(TRUE,3),rep(FALSE,5)),
                       c(rep(FALSE,4), rep(TRUE,3),rep(FALSE,5)),# Select summer months
                       c(rep(FALSE,4), rep(TRUE,3),rep(FALSE,5)))
period1s  <- list(c(130:258),c(130:258),c(130:258)) # Select years before 1981
period2s  <- list(c(1:117),c(1:117),c(1:117)) # Select years after 1980
seasons   <- c("summ","summ","summ")
var_names <- c("tmax_may_jul","tmin_may_jul","pptmean_may_jul")
# Template to change the projection of climatic layers to the projection of the LULC layers
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/CONUS_LandUse_1938_1979_2009_2018_10km/may2023")
template <- rast(list.files( ,pattern = '.tif')[1])[[1]] # Map created for LULC
# Loop to pass the var_calc function
# i = 1 -> Maximum Temperature for summer; 
# i = 2 -> Minimum Temperature for summer; 
# i = 3 -> Mean Precipitation for summer
for(i in 1:3){
  directory1    <- directories1[i]
  directory2    <- directories2[i]
  select_files2 <- select_files2s[[i]]
  period1       <- period1s[[i]]
  period2       <- period2s[[i]]
  season        <- seasons[i]
  var_name      <- var_names[i]
  print (var_name)
  var_calc(directory1, directory2, select_files2, period1, period2, season, var_name)
}
#--------------------------------------------------------------------------------------------------------------------------------------
# Extracting the mean value of climatic variable for each route using mean function
#--------------------------------------------------------------------------------------------------------------------------------------
# Period 1990 - 2019
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/Climatic_vars_1938_1979_2009_2018_10km/may2023")
files_j  <- rast(list.files( , pattern = '_july2024.tif')) # [j]
pptmean_may_jul <- files_j [[1:82]][[53:82]] # 1990-2019
tmax_may_jul    <- files_j [[83:164]][[53:82]] # 1990 -2019
tmin_may_jul    <- files_j [[165:246]][[53:82]] # 1990 -2019
layers <- list(pptmean_may_jul,tmax_may_jul,tmin_may_jul)
name_layers <- c('pptmean_may_jul','tmax_may_jul','tmin_may_jul')
# Template for North America
setwd("~/Desktop")
world <- read_sf(dsn = ".", layer = "countriesWGS84")
# Template to use for transforming the CRS and extent of US shapefile to the CRS 
# of climatic and LULC variables:
america <- world [c(252), ] 
america <- st_transform(america, crs = crs(layers[[1]])) 
ext <- ext(-2357953, 2292047, 235292.6, 3175293)
america <- st_crop(america, ext)
# BBS routes layer
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity")
# Sampling units
data        <- read.csv('taxonomic_bdiv_ch_1990_2019_nomore1missing_data_July2024.csv', h = T)
sampl_units <- data[ ,c(1:3)][!duplicated(data$CountryStateRoute),]
# Transforming the matrix into spatial points data frame, so it can be used to 
# extract env. variable values
coord_units <- st_as_sf(sampl_units, 
                        coords = c("Longitude","Latitude"), 
                        crs = crs("+proj=longlat +datum=WGS84 +no_defs"))
# Transforming the coordinates of my point to match the environmental layers
transf_coord_units <-  st_transform(coord_units, crs(layers[[1]]))
rm(files_j,pptmean_may_jul,tmax_may_jul,tmin_may_jul,world,ext,data,coord_units)
# Set the buffer resolutions to be used in the loops below
resol     <- c("12_5","25","50")
buff_size <-list(st_buffer(transf_coord_units,12500),
                 st_buffer(transf_coord_units,25000),
                 st_buffer(transf_coord_units,50000))
for (j in 1: length(layers)){
  names(layers[[j]]) <- paste(name_layers[j], c(1990:2019), sep = '_')
  layer_j <- mask (layers[[j]], america)
  env_dt <- data.frame()
  for (i in 1: length(resol)){
    print(resol[i])
    resol_i     <- resol[i]
    buff_size_i <- buff_size[[i]]
    env_var_dtf <-  data.frame(extract(layer_j,
                                       buff_size_i, 
                                       method = 'simple', 
                                       fun    = mean, 
                                       na.rm  = TRUE,
                                       bind   = T,
                                       xy     = T))                
    env_var_dtf$resoln <- resol_i
    env_var_dtf <- env_var_dtf[ ,c(1,32,2:31)]
    env_dt <- rbind(env_dt, env_var_dtf)
  }
  sampl_units <- full_join(sampl_units, env_dt)
}
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity")
write.csv(sampl_units, 'clim_vars_mean_before_temp_trend_july2024.csv', row.names = F)
# ------------------------------------------------------------------------------------------------------------------
# Univariate rate of change in climate
# ------------------------------------------------------------------------------------------------------------------
# Temporal trends
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity")
clim_vars <- read.csv('clim_vars_mean_before_temp_trend_july2024.csv')
long <- clim_vars%>% 
  pivot_longer(
    cols = c("pptmean_may_jul_1990":"pptmean_may_jul_2019", 
             "tmax_may_jul_1990" :"tmax_may_jul_2019", 
             "tmin_may_jul_1990" :"tmin_may_jul_2019"), 
    names_to = "Year",
    values_to = "value"
  )
long$clim_var <- sub("(_[0-9]+).*", "\\2", long$Year)
long$Year<- as.integer(gsub('^(?:[^_]*_){3}','',long$Year))
temp_trend <- data.frame(long %>%
                           group_by(Longitude, Latitude,CountryStateRoute, resoln, clim_var) %>%
                           nest() %>%
                           mutate(fit = map(data, ~ lm(value ~ Year, data = .x)),
                                  tidied = map(fit, tidy)) %>%
                           unnest(tidied, .drop = TRUE))[ ,-c(6,7)] %>%
  filter(term =='Year')
wide <- arrange(dcast(temp_trend[ ,-6], 
                      Longitude + Latitude + CountryStateRoute ~ clim_var + resoln, value.var = 'estimate'), 
                CountryStateRoute)
write.csv(wide,'clim_vars_lm_tempTrend_1990_2019_july2024.csv', row.names = F)
# ------------------------------------------------------------------------------------------------------------------
# Correlation between Tmax and Tmin
# ------------------------------------------------------------------------------------------------------------------
long <- clim_vars [ ,c(1:4,35:94)]%>% 
  pivot_longer(
    cols = c("tmax_may_jul_1990" :"tmax_may_jul_2019", "tmin_may_jul_1990" :"tmin_may_jul_2019"), 
    names_to = "Year",
    values_to = "value"
  )
long$clim_var <- sub("(_[0-9]+).*", "\\2", long$Year)
long$Year<- as.integer(gsub('^(?:[^_]*_){3}','',long$Year))

wide <- arrange(dcast(long, Longitude + Latitude + CountryStateRoute + resoln + Year ~ clim_var, value.var = 'value'), CountryStateRoute)
cor_test <- data.frame(wide  %>%
                         nest(- c(Longitude, Latitude,CountryStateRoute, resoln)) %>%
                         mutate(test = map(data, ~ cor.test(.x$tmax_may_jul, .x$tmin_may_jul)),
                                tidied = map(test, tidy))%>% 
                         unnest(tidied, .drop = TRUE))