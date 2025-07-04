#--------------------------------------------------------------------------------------------------------------------------------------
# PACKAGES
#--------------------------------------------------------------------------------------------------------------------------------------
require(betapart)
require(tidyverse)
require(magrittr)
require(data.table)
require(broom) 
#--------------------------------------------------------------------------------------------------------------------------------------
# FUNCTIONS:
#--------------------------------------------------------------------------------------------------------------------------------------
# Function to transform the beta.part output into a data frame format
bpair_data <- function(x){
  # Turn the .dist format to matrix
  x <- as.matrix(x)
  # Change diagonal values to NA
  diag(x) <- NA
  # Convert to data frame
  x <- data.frame(x)
  # Set the row names as column to keep track of the cell number
  setDT(x, keep.rownames = TRUE)[]
  # Convert to data frame again
  x <- data.frame(x)
  return(x)
}
# Function to eliminate rows and columns fully filled with NA, and to calculate 
# the mean index value per cell
beta_mean <- function(x) {  
  x$index_mean <- rowMeans(x[,2:ncol(x)], na.rm = T)
  return(x)
}
# Function to calculate the temporal trend of beta-diversity and its components
beta_temp_trend <-function(j){
  beta_lm <- data.frame(j %>%
                          nest(-c(CountryStateRoute, Longitude, Latitude, Ecor_Lev1_name))%>%
                          mutate(
                            test = map(data, ~ lm(index ~ Year, data = .x)),
                            tidied = map(test, tidy)
                          ) %>%
                          unnest(tidied, .drop = TRUE) %>% filter(term == 'Year'))
  beta_lm <- beta_lm %<>%
    mutate(signif_lm = as.factor(case_when(
      p.value < 0.05 & estimate > 0 ~ "signif increasing",
      p.value < 0.05 & estimate < 0 ~ "signif decreasing",
      p.value >= 0.05 ~ "non_signif"
    )))
  return(beta_lm)
}
#--------------------------------------------------------------------------------------------------------------------------------------
# Calculating mean pairwise beta diversity per cell per year  
#--------------------------------------------------------------------------------------------------------------------------------------
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/Results_2023May/BBS_community_data_per_year_1991_2019_nomore1missing")
files <- sapply(list.files( , pattern = 'July2024'), read.csv, simplify = FALSE) 
years <- c(1990:2019)
# Mean beta diversity per cell per year 
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/BBS_taxonomic_bpair_per_year_1991_2019_nomore1missing")
# Beta div. per ecoregion
  # Ecoregion codes:
    # 5 -> Northern Forests  
    # 6 -> Northwestern Forested Mountains 
    # 8 -> Eastern Temperate Forests 
    # 9 -> Great Plains  
    # 10 -> North American Deserts 
dt <- unique(files[[1]][ ,c(6,7)] %>% 
               filter(Ecor_Lev1_code %in% c(5,6,8,9,10)))
for (i in 1: length(files)){
  # Year
  name_i <- years[i]
  print(name_i)
  # Select the file corresponding to the year i
  year_i_0 <- files[[i]]
  for (j in 1:nrow(dt)) { # Calculations per ecoregions
    print(dt[j,2])
    # Eliminate the columns that are not presence absence, so the beta.pair function works
    year_j <- year_i_0 %>% filter(Ecor_Lev1_code == dt[j,1])
    # Set the row names as column to keep track of the cell number
    year_j  <- setDT(year_j , keep.rownames = TRUE)[]
    year_i <- year_j[ ,-c(1:8)]
    year_i[ year_i > 0] <- 1 
    # Calculating pairwise beta-diversity and its partitions using the betapart functions from Baselga 2012
    year_i_bpair_sor <- beta.pair(year_i, index.family = "sorensen")
    # Transforming the beta.part output into a data frame format
    year_i_bpair_sor <- sapply(year_i_bpair_sor, bpair_data, simplify = FALSE)
    # Eliminating rows and columns fully filled with NA, and to calculating the mean index value per cell
    bpair <- lapply(year_i_bpair_sor, beta_mean)
    # Change the name of last column by the name of data frame in list
    bpair <- bpair %>%
      imap(~ {nm1 <- .y
      .x %>%
        rename_at(ncol(bpair[[1]]), ~ nm1) })
    # Selecting only the row name column and the index average column for each dataframe
    bpair <- lapply(bpair, function(x) x[ ,c(1, ncol(x))])
    # Joining data frames in list into a single data frame with the 3 index averages
    bpair <-  bpair %>% reduce(full_join, by = "rn")
    # Including cell number and cell coordinates
    bpair <- full_join(year_j[ ,1:8], bpair)
    # Including year
    bpair$Year <- name_i
    # Save
    write.csv(bpair, paste(name_i,"_",dt[j,2], "_taxonomic_bdiversity_pair_July2024.csv", sep = ''), row.names = F)
  }
}
#--------------------------------------------------------------------------------------------------------------------------------------
# Temporal trends of total assemblage dissimilarity and its components per BBS route
#--------------------------------------------------------------------------------------------------------------------------------------
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/BBS_taxonomic_bpair_per_year_1991_2019_nomore1missing")
files <- sapply(list.files( , pattern = 'July2024'), read.csv, simplify = FALSE)
ecor_names <- c("Eastern_Temperate_Forests","Great_Plains","North_American_Deserts",
                "Northern_Forests","Northwestern_Forested_Mountains")
mean_bdiv <- data.frame()
bdiv_dt   <- data.frame()
for (i in 1:length(ecor_names)){
  print(ecor_names[i])
  files_i <- files[seq(from = i, to = length(files), by = 5)]
  files_i <- arrange(bind_rows(files_i, .id = "column_label"), CountryStateRoute, Year)[ ,-c(1,2)]
  mean_bdiv_year <- files_i %>% 
    group_by(Year) %>%
    summarise(mean_beta.sim = mean (beta.sim),
              sd_beta.sim = sd(beta.sim),
              mean_beta.sne = mean (beta.sne), 
              sd_beta.sne = sd (beta.sne),
              mean_beta.sor = mean (beta.sor),
              sd_beta.sor = sd (beta.sor))
  mean_bdiv_year$Ecor_Lev1_name <- ecor_names[i]
  mean_bdiv <- rbind(mean_bdiv, mean_bdiv_year)
  
  files_i <- list('beta.sim' = files_i[ ,c(1:8,11)], 'beta.sne' = files_i[ ,c(1:7,9,11)], 'beta.sor' = files_i[ ,c(1:7,10:11)])
  # Change the name of the fifth column to "index" so the function beta_temp_trend works properly
  files_i <- lapply(files_i, setNames, c("Longitude", "Latitude", "CountryStateRoute", "BCR", "na_elevation_m" ,"Ecor_Lev1_code", "Ecor_Lev1_name", "index", "Year"))
  # Calculate the temporal trend of beta-diversity and its components using function beta_temp_trend
  beta.pair <- lapply(files_i, beta_temp_trend)
  # Create a column in each data frame with the name of the corresponding metric 
  # as value (the name of the metric is the name of the data frame)
  beta.pair <- Map(cbind, beta.pair, beta.index = names(beta.pair))
  beta.pair <- arrange(bind_rows(beta.pair), CountryStateRoute,beta.index)
  beta.pair$Ecor_Lev1_name <- ecor_names[i]
  beta.pair <- left_join(beta.pair, unique(files_i[[1]][ ,c(3,5)])) 
  bdiv_dt   <- rbind(bdiv_dt, beta.pair)
}
bdiv_dt <- bdiv_dt[ ,c(1:3,14,4,8:13)]
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity")
write.csv(bdiv_dt, 'taxonomic_bdiv_ch_1990_2019_nomore1missing_data_July2024.csv', row.names = F)
write.csv(mean_bdiv,'mean_taxonomic_bdiv_per_year_1990_2019_nomore1missing_data_July2024.csv', row.names = F)
#--------------------------------------------------------------------------------------------------------------------------------------
# Temporal trends of species richness per BBS route
#--------------------------------------------------------------------------------------------------------------------------------------
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity/Results_2023May/BBS_community_data_per_year_1991_2019_nomore1missing")
files <- sapply(list.files( , pattern = 'July2024'), read.csv,simplify = FALSE)
alpha.temp <- data.frame(rn = as.character(c(1:nrow(files[[1]]))), files[[1]][ ,1:7])
for (i in 1:length(files)) {
  t_i <-  files[[i]]
  z <- strsplit(names(files[i]), '_')[[1]][1]
  sp_rich <- t_i %>% 
    mutate(sp_rich = rowSums(t_i[ ,8:ncol(t_i)], na.rm = TRUE))
  sp_rich <- sp_rich[ ,c(1:7, ncol(sp_rich))]
  names(sp_rich) <- c("Longitude","Latitude","CountryStateRoute","BCR", "na_elevation_m" ,"Ecor_Lev1_code", "Ecor_Lev1_name", 
                      paste("sp_rich_", z, sep = ''))
  alpha.temp <- full_join(alpha.temp, sp_rich)
}
alpha.temp1 <- arrange(melt(alpha.temp, id = 1:8, 
                            variable.name = "Year",
                            value.name = "index"), 
                       CountryStateRoute)
alpha.temp1$Year <- as.numeric(sapply(strsplit(as.character(alpha.temp1$Year), "_y"), "[", 2))
# Temporal trends
alpha.temp_lm <- beta_temp_trend(alpha.temp1)
alpha.temp_lm <- full_join(alpha.temp_lm, alpha.temp.sd)
table(alpha.temp_lm$signif_lm)
alpha.temp_lm <- alpha.temp_lm [ ,-c(5:7)]
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity")
write.csv(alpha.temp_lm, 'taxonomic_alpha_itself_ch_1990_2019_nomore1missing_July2024.csv', row.names = F)
write.csv(alpha.temp1, 'taxonomic_alphadiv_itselt_per_year_1990_2019_nomore1missing_data_July2024.csv', row.names = F)