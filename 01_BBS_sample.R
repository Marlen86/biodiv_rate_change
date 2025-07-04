#--------------------------------------------------------------------------------------------------------------------------------------
# PACKAGES
#--------------------------------------------------------------------------------------------------------------------------------------
require(tidyverse)
require(reshape2)
require(spThin)
require(stringr)
require(readxl)
require(terra)
require(sf)
#--------------------------------------------------------------------------------------------------------------------------------------
# FUNCTIONS:
#--------------------------------------------------------------------------------------------------------------------------------------
# Function to create the community data of 0 and 1
comm_data <- function(x){
  # Select only useful columns
  x <- x[ ,c(2,7:ncol(x))]
  # Turn data to wide format, so species are columns
  x_wide <- reshape2::recast(x,
                             Longitude + Latitude + CountryStateRoute + BCR + na_elevation_m + Ecor_Lev1_code + Ecor_Lev1_name ~ variable + English_Common_Name,
                             id.var = c("English_Common_Name", "Longitude", "Latitude", "CountryStateRoute", "BCR", "na_elevation_m","Ecor_Lev1_code", "Ecor_Lev1_name"))
  # Order rows
  x_wide <- arrange(x_wide, CountryStateRoute)
  # NA values to 0
  x_wide[is.na(x_wide)] <- 0 
  # Eliminate prefix from species names
  names(x_wide) <- sub("SpeciesTotal_", "", names(x_wide))
  return(x_wide)
  }
# Function to complete data frames with the species (columns) and the routes 
# (rows) that were not recorded for each year 
complete_dtf <- function(x){
  # Complete columns with species names
  list.colname <- lapply(bbs_list, function(x) {data.frame(names(x))})
  list.colname <- bind_rows(list.colname, .id = "column_label")
  list.colname <- unique(list.colname$names.x.)
  list.colname <- c(list.colname[1:7], sort(list.colname[8:length(list.colname)]))
  x[setdiff(list.colname, names(x))] <- 0 
  x <- x[list.colname] 
  return(x)
}
#--------------------------------------------------------------------------------------------------------------------------------------
# BBS data
#--------------------------------------------------------------------------------------------------------------------------------------
# Reading BBS original files 
setwd("/Volumes/G2T/memory/Thesis/2022_BBS_data/States")
datf1_bbs <- unique(bind_rows(lapply(list.files( , pattern = '.csv'), read.csv)))
datf1_bbs$CountryStateRoute <- paste(datf1_bbs$CountryNum,'_', 
                                     datf1_bbs$StateNum, '_',
                                     datf1_bbs$Route, sep = '')
#--------------------------------------------------------------------------------------------------------------------------------------
# Selecting routes present in more than 24 years
#--------------------------------------------------------------------------------------------------------------------------------------
routes_sampled <- unique(datf1_bbs[ ,c(1:6,15)]) %>% 
  filter(Year > 1989,
         Year < 2020, 
         RPID == 101, 
         CountryNum == 840,
         StateNum != 3 # Exclude Alaska
         )
routes_sampled_freq <- data.frame(routes_sampled %>% 
                                    group_by(CountryStateRoute) %>% 
                                    count() %>% 
                                    filter(n > 24))
routes_sampled <- routes_sampled %>% 
  filter(CountryStateRoute %in% routes_sampled_freq$CountryStateRoute)
# To eliminate those routes that were not sampled in 1 or more consecutive years
# Filling with NA years with no data
routes <- data.frame(unique(arrange(routes_sampled %>%
                                      complete(CountryStateRoute, 
                                               nesting(Year = c(min(routes_sampled$Year): max(routes_sampled$Year)))),
                                    CountryStateRoute, Year)))
# To eliminate routes not sampled in 1990
routes90 <- routes %>% 
  filter(Year == 1990 & is.na(CountryNum))
# Check runs of NA 
r_routes <- rle(is.na(routes$Route))
# Select routes with > one years of no sampling or not sampled in 1990 to eliminate them
routes <- routes [rep(r_routes$values & r_routes$lengths > 1, r_routes$lengths), ]
routes <- rbind(routes, routes90)
routes <- data.frame(routes_elim = unique(routes$CountryStateRoute))
# Eliminating from the sample the BBS routes not sampled in more than 1 consecutive years
routes_selected <- routes_sampled %>% 
  filter(!(CountryStateRoute %in% routes$routes_elim))
routes_year9019 <- routes_selected %>% 
  group_by(CountryStateRoute) %>% 
  count() %>% 
  filter(n > 24)
rm(r_routes, routes, routes_sampled, routes_sampled_freq, routes90, routes_selected)
# Thinning routes so we have only one route per map cell
# geographic coordinates of routes
setwd("/Volumes/G2T/memory/Thesis/2022_BBS_data")
routes <- read.csv("routes.csv", h = T)
routes <- routes %>% 
  filter(StateNum != 3) # Exclude Alaska
routes$CountryStateRoute <- paste(routes$CountryNum,'_',routes$StateNum,'_',routes$Route, sep = '')
routes <- unique(routes[ ,c(6,7,9,12)]) %>% 
  filter(CountryStateRoute %in% routes_year9019$CountryStateRoute)
routes <- na.omit(unique(right_join(routes, routes_year9019)))
routes <- unique(routes[,c(1:4)])
routes$Route <-"BBS" 
set.seed(1)
routes_thinned <- thin(routes,
                       lat.col = "Latitude", 
                       long.col = "Longitude", 
                       spec.col = "Route",
                       thin.par = 15, 
                       reps = 100,
                       locs.thinned.list.return = T,
                       write.files = F, 
                       max.files = 1,
                       out.dir, 
                       out.base = "thinned_data", 
                       write.log.file = F, 
                       log.file = "spatial_thin_log.txt",
                       verbose = TRUE)
routes_thinned <- right_join(routes,routes_thinned[[1]])
#--------------------------------------------------------------------------------------------------------------------------------------
# Creating bird assemblage data
#--------------------------------------------------------------------------------------------------------------------------------------
# Species data
setwd("/Volumes/G2T/memory/Thesis/2022_BBS_data")
sp_list_bbs0 <- read.csv("SpeciesList.csv", h = T)[ ,-c(1,4:7)] # This data frame was produced from the original 'sp_list.txt'
# Eliminating unidentified sp, subspecies and hybrids
sp_list_bbs0 <- dplyr::filter(sp_list_bbs0, !grepl("unid", English_Common_Name))
sp_list_bbs0 <- dplyr::filter(sp_list_bbs0, !grepl("Unid", English_Common_Name))
sp_list_bbs0 <- dplyr::filter(sp_list_bbs0, !grepl("hybrid", English_Common_Name))
sp_list_bbs0 <- unique(dplyr::filter(sp_list_bbs0, !grepl("Hybrid", English_Common_Name)))
# Getting rid of subspecies names
sp_list_bbs0$English_Common_Name <- gsub("\\s*\\([^\\)]+\\)","", as.character(sp_list_bbs0$English_Common_Name))
sp_list_bbs0$English_Common_Name <- trimws(sp_list_bbs0$English_Common_Name, "l")
sp_list_bbs0$Species <- gsub("\\s*\\([^\\)]+\\)","", as.character(sp_list_bbs0$Species))
sp_list_bbs0$Species <- sub(" .*", "", sp_list_bbs0$Species)
# Changing the string after the '-' to lowercase
sp_list_bbs0$English_Common_Name <- gsub("(-.)","\\L\\1", sp_list_bbs0$English_Common_Name, perl = TRUE)
sp_list_bbs0$SCI_NAME <- paste(sp_list_bbs0$Genus, ' ', sp_list_bbs0$Species, sep = '')
sp_list_bbs0 <- unique(sp_list_bbs0)
# Joining BBS data with species names
datf1_bbs <- right_join(sp_list_bbs0, datf1_bbs)
# Eliminating data from unidentified sp and hybrids  
datf1_bbs <- na.omit(datf1_bbs)
datf1_bbs$English_Common_Name <- as.factor(datf1_bbs$English_Common_Name)
rm(sp_list_bbs0)
# Selecting the BBS data that have been sampled in 25 years or more
bbs <- datf1_bbs[ , c(1:5,11,18:19)] %>% 
  filter (CountryStateRoute %in% routes_thinned$CountryStateRoute,
          Year > 1989, 
          Year < 2020)
bbs <- left_join(bbs, routes_thinned[ ,-5])
# Selecting species sampled in more than 4 routes each year
bbs_sp_years <- unique(bbs[ ,c(2,6,8)]) %>%
  group_by(English_Common_Name, Year) %>%
  summarise(numb_routes = n()) %>%
  filter(numb_routes < 4) 
# Exclude species sampled in less than 5 routes any year
bbs <- bbs %>% 
  filter(!(English_Common_Name %in% bbs_sp_years$English_Common_Name))
rm(datf1_bbs,routes, routes_year9019, routes_thinned)

# Eliminating water/sea/shore birds using BBS website and Peterjohn and Sauer 1993
# Getting land bird names from BBS website
setwd("/Volumes/G2T/memory/Thesis")
bbs1 <- read.csv("sp_breeding_habitat_BBS.csv")                                                                
bbs2 <- read.csv("sp_migratory_category_BBS.csv")                                                              
bbs3 <- read.csv("sp_nest_location_BBS.csv")                                                                   
bbs4 <- read.csv("sp_nest_type_BBS.csv")        

bbs_sp <- full_join(bbs1,bbs2)
bbs_sp <- full_join(bbs_sp,bbs3)
bbs_sp <- full_join(bbs_sp,bbs4)
rm(bbs1,bbs2,bbs3,bbs4)
# Changing the string after the '-' to lowercase
bbs_sp$English_Common_Name <- gsub("(-.)","\\L\\1", bbs_sp$English_Common_Name, perl = TRUE)
# Changing abbreviations to full words
bbs_sp$English_Common_Name <- gsub("(crn.)","crowned ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(&)","and", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(Am. )","American ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(bell. )","bellied ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(crest. )","crested ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(col. )","collared ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(back. )","backed ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(bkd. )","backed ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(tail. )","tailed ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(tld. )","tailed ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(Gol.-)","Golden-", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(th. )","throated ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(thr. )","throated ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(chin. )","chinned ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(Grt.) ","Great ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(head. )","headed ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(Pac.-sl.-co. )","Pacific-slope ", bbs_sp$English_Common_Name, perl = TRUE)
bbs_sp$English_Common_Name <- gsub("(N. )","Northern ", bbs_sp$English_Common_Name, perl = TRUE)
# Eliminating white space at the begining of names
bbs_sp$English_Common_Name <- str_trim(bbs_sp$English_Common_Name, "left") 
# Uppercase word after white space
bbs_sp$English_Common_Name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",    
                                    bbs_sp$English_Common_Name ,
                                    perl = TRUE)
bbs_sp <- bbs_sp %>% 
  filter(Breeding_Habitat %in% c(NA, 'Grassland', 'Successional-scrub', 'Urban',  'Woodland species'))
# Getting land bird names from Peterjohn and Sauer 1993
setwd("/Volumes/G2T/memory/Thesis/CH1_range_correlation/bibliography/Migratory categ_Peterjohn_1993")
residents               <- read_excel("permanent_residents.xlsx")
colnames(residents)[4]  <- "English_Common_Name"
residents$migrat_categ  <- 'resident'
short_dist              <- read_excel("short_distance_migrants.xlsx")
colnames(short_dist)[1] <- "English_Common_Name"
short_dist$migrat_categ <- 'short_dist'
neot                    <- read_excel("neotropical_migrants.xlsx")
colnames(neot)[3]       <- "English_Common_Name"
neot$migrat_categ       <- 'neot'
peterjohn <- rbind(residents[ ,c(4,6)], short_dist[ ,c(1,7)], neot[ ,c(3,9)])
rm(residents, short_dist,neot)
# Modifying names to match the bird surveys
# Changing the string after the '-' to lowercase
peterjohn$English_Common_Name <- gsub("(-.)","\\L\\1", peterjohn$English_Common_Name, perl = TRUE)
# Changing abbreviations to full words
peterjohn$English_Common_Name <- gsub("(bell. )","bellied ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(crest. )","crested ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(col. )","collared ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(back. )","backed ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(bkd. )","backed ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(tail. )","tailed ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(tld. )","tailed ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(Gol.-)","Golden-", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(th. )","throated ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(thr. )","throated ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(chin. )","chinned ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(Grt.) ","Great ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(head. )","headed ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(Pac.-sl.-co. )","Pacific-slope ", peterjohn$English_Common_Name, perl = TRUE)
peterjohn$English_Common_Name <- gsub("(N. )","Northern ", peterjohn$English_Common_Name, perl = TRUE)
# Standardizing some species names
peterjohn$English_Common_Name[peterjohn$English_Common_Name == "Le Conte's Sparrow"] <- "LeConte's Sparrow"
peterjohn$English_Common_Name[peterjohn$English_Common_Name == "Le Conte's Thrasher"] <- "LeConte's Thrasher"
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Gray Jay'] <- 'Canada Jay'
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Sky Lark'] <- "Eurasian Skylark"
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Western Scrub-jay'] <- "California Scrub-jay"
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Common Ground-dove'] <- "Common Ground Dove"
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Ruddy Ground-dove'] <- "Ruddy Ground Dove"
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Magnificent Hummingbird'] <- "Rivoli's Hummingbird"
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Blue-throated Hummingbird'] <- "Blue-throated Mountain-gem"
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Western Swamphen'] <- "Purple Swamphen"
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Three-toed Woodpecker'] <- "American Three-toed Woodpecker"
peterjohn$English_Common_Name[peterjohn$English_Common_Name == "Nelson's sharp-tailed Sparrow"] <- "Nelson's Sparrow"
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Rock Dove'] <- 'Rock Pigeon'
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Sage Grouse'] <- 'Greater Sage-grouse'
peterjohn$English_Common_Name[peterjohn$English_Common_Name == 'Whip-poor-will'] <- 'Eastern Whip-poor-will'

# Merging the two lists
bbs_sp <- arrange(full_join(bbs_sp, peterjohn), English_Common_Name)
rm(peterjohn)
# Eliminate other waterbirds
bbs_sp <- bbs_sp %>% filter(!(English_Common_Name %in% c("Common Merganser", "Killdeer", "Long-billed Curlew")))
# "Upland Sandpiper" was not eliminated since is a praire obligated species
data_names <- arrange(unique(data.frame(English_Common_Name = bbs$English_Common_Name)), English_Common_Name)
data_names1 <- inner_join(data_names,bbs_sp)
setdiff(data_names$English_Common_Name,data_names1$English_Common_Name)
# Add three land birds that are not in this lists but are in the BBS
to_add <- data.frame(English_Common_Name = c("Cordilleran Flycatcher", "Sagebrush Sparrow", "Woodhouse's Scrub-jay"),
                     Breeding_Habitat = c(NA,NA,NA),
                     Migratory_Category = c(NA,NA,NA),
                     Nest_Location = c(NA,NA,NA),
                     Nest_Type = c(NA,NA,NA), 
                     migrat_categ = c(NA,NA,NA))
data_names1 <- arrange(rbind(data_names1, to_add), English_Common_Name)
data_names1 <- data.frame(English_Common_Name = unique(data_names1$English_Common_Name))

bbs <- bbs %>% 
  filter(English_Common_Name %in% data_names1$English_Common_Name)
rm(bbs_sp, data_names, data_names1, to_add, bbs_sp_years)
# Convert species counts into 0/1 presence absence data
bbs$SpeciesTotal[bbs$SpeciesTotal > 0] <- 1 
#---------------------------------------------------------------------------------------------------------------------
# Include ecoregion and elevation
#---------------------------------------------------------------------------------------------------------------------
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/NA_Terrestrial_Ecoregions_v2_Level_I_Shapefile/NA_TerrestrialEcoregions_LI/data")
ecoregions <- read_sf(dsn = ".", layer = "NA_Terrestrial_Ecoregions_v2_level1")
ecoregions$LEVEL1 <- as.integer(ecoregions$LEVEL1)
# Get the coordinates for BBS routes
sampl_units <- bbs[ ,8:10][!duplicated(bbs$CountryStateRoute),]
# Transforming the matrix into spatial points data frame, so it can be used to 
# extract env. variable values
coord_units <- st_as_sf(sampl_units, coords = c("Longitude","Latitude"), 
                        crs = crs("+proj=longlat +datum=WGS84 +no_defs"))
transf_coord_units <-  st_transform(coord_units, crs(ecoregions))
ecor_dt <- data.frame(st_intersection(ecoregions, transf_coord_units))
# Elevation data
setwd("/Volumes/G2T/memory/Thesis/Elevation_TIF/NA_Elevation/data")
# Distance in meters, < 0 below sea level, > 0 above sea level
elev_lyr <- rast("na_elevation.tif")
# Extract data for the BBS routes 
elev_var_dtf <-  data.frame(extract(elev_lyr,
                                    transf_coord_units, 
                                    method = 'simple', 
                                    na.rm  = TRUE,
                                    raw    = F, 
                                    bind   = T,
                                    xy     = T))  

names(elev_var_dtf) <- c("CountryStateRoute", "na_elevation_m", "Longitude", "Latitude")
elev_var_dtf_na <- elev_var_dtf[is.na(elev_var_dtf$na_elevation_m),]
# Merge BBS data with elevation data
bbs <- na.omit(left_join(bbs, elev_var_dtf[ ,c(1,2)])) # eliminate routes with no elevation data 
# Merge BBS data with ecoregion data
bbs <- left_join(bbs, ecor_dt[ ,c(1,2,6)])
names(bbs)[names(bbs) == 'LEVEL1'] <- 'Ecor_Lev1_code'
names(bbs)[names(bbs) == 'NameL1_En'] <- 'Ecor_Lev1_name'
# Keep only ecoregions with more than 50 BBS routes
table(ecor_dt$NameL1_En)
bbs <- bbs %>% 
  filter(Ecor_Lev1_name %in% c('Eastern Temperate Forests', 'Great Plains', 'North American Deserts',
                               'Northern Forests', 'Northwestern Forested Mountains'))
rm(coord_units,ecor_dt, ecoregions, elev_var_dtf,elev_lyr, elev_var_dtf_na,sampl_units, transf_coord_units) 
#---------------------------------------------------------------------------------------------------------------------
# Split the data set to a list in which each year is a data frame
#---------------------------------------------------------------------------------------------------------------------
bbs_list <- split(bbs, bbs$Year)
# Create community data for each year
bbs_list <- lapply(bbs_list, comm_data)
bbs_list <- sapply(bbs_list, complete_dtf, simplify = FALSE)
# how many routes in single year
range(lapply(bbs_list, nrow))
# Save community data
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/Results_2023May/BBS_community_data_per_year_1991_2019_nomore1missing")
for(i in 1:length(bbs_list)){
  write.csv(bbs_list[[i]], paste0("y",names(bbs_list)[i],"_BBS_species_presence_July2023.csv"), row.names = F)
}