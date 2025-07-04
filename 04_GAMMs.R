# ------------------------------------------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------------------------------------------
require(tidyverse)
require(magrittr)
require(broom)
require(data.table)
require(mgcv) # GAMMs
# ------------------------------------------------------------------------------------------------------------------
# FUNCTION:
# ------------------------------------------------------------------------------------------------------------------
# Function for the gamm of biodiversity as predicted by environmental rate of change
env_gamm <- function(data, use){ 
  # GAMM set up
    # To be able to exclude the correlated variables and/or the LULC variable with 
    # the lowest % of coverage for each region we utilize the variable named "use" 
    # and we create a model formula as below in the "form" object:
    form <- as.formula(paste0("log(estimate+1)~", paste0("s(",use,")", collapse = " + ")))
    gamm_mod <- mgcv::gamm(formula = form,
                           data = data, family = gaussian,
                           correlation = corGaus(value = 0.5, 
                                                 form = ~ Latitude + Longitude, 
                                                 fixed = F))
    # If needed, we used control = list(MaxIter = 100) to handle model convergence. 

  # Summary table for the model:
    # Parametric portion 
    model_summary_lme <- setDT(data.frame(summary(gamm_mod$lme)$tTable), keep.rownames = TRUE)[] 
    model_summary_lme <- model_summary_lme[ ,c(1,2,5,6)]
    names(model_summary_lme)    <- c("predictor", "value","t value","p_value_parametric")
    model_summary_lme$predictor <- gsub("X","",model_summary_lme$predictor)
    model_summary_lme$predictor <- gsub("Fx1","",model_summary_lme$predictor)
    # Smooth parameters portion
    model_summary_gam           <- setDT(data.frame(summary(gamm_mod$gam)$s.table), keep.rownames = TRUE)[] 
    model_summary_gam           <- model_summary_gam[ ,c(1,2,4,5)]
    names(model_summary_gam)    <- c("predictor", "edf","F value","p_value_smooth_terms")
    model_summary <- full_join(model_summary_lme, model_summary_gam)
    # Add other evaluation statistics to table
    model_eval    <- data.frame(predictor = c('AIC', 'cor_str_param', 'r_squared', 'N'),
                                value     = c(round(AIC(gamm_mod$lme), digits = 2), round(summary(gamm_mod$lme)$modelStruct$corStruct, digits = 2), round(summary(gamm_mod$gam)$r.sq, digits = 4), round(summary(gamm_mod$lme)$dims$N[1], digits = 0)))
    model_summary <- full_join(model_summary, model_eval)
    # Add identifiers
    model_summary$index     <- unique(data[,11])
    model_summary$lulc_scen <- lulc_scen[i]
    model_summary$resoln    <- resoln[j]
    model_summary$ecoregion <- ecoregion[k]
    return(model_summary)
  # GAMM model check
    type <- "deviance"
    resid <- residuals(gamm_mod$gam, type = type)
    linpred <- napredict(gamm_mod$gam$na.action, gamm_mod$gam$linear.predictors)
    observed.y <- napredict(gamm_mod$gam$na.action, gamm_mod$gam$y)
    
    tiff(paste(unique(data[,11]), ecoregion[k],lulc_scen[i],resoln[j],
               'km_gam_check.tiff', sep = '_'), 
         units = "mm", width = 200, height = 250, res = 300)
    par(mfrow=c(2,2))
    # QQ plot
    print(qq.gam(gamm_mod$gam, pch=19,cex=.3, ))
    # Residuals vs linear predictor
    print(plot(linpred, resid, main = "Resids vs. linear pred.", 
                        xlab = "linear predictor", ylab = "residuals"))
    # Histogram of residuals
    print(hist(resid, xlab = "Residuals", main = "Histogram of residuals"))
    # Observed vs fitted values
    print(plot(fitted(gamm_mod$gam), observed.y, xlab = "Fitted Values", 
         ylab = "Response", main = "Response vs. Fitted Values"))
    dev.off()
    
    tiff(paste(unique(data[,11]), ecoregion[k],lulc_scen[i],resoln[j],
               'km_gamm_results.tiff', sep = '_'), 
         units = "mm", width = 200, height = 250, res = 300)
    print(plot(gamm_mod$gam,pages=1))
    dev.off()
}
# ------------------------------------------------------------------------------------------------------------------
# DATA
# ------------------------------------------------------------------------------------------------------------------
# Bio data
setwd("/Volumes/G2T/memory/Thesis/CH2_funct_diversity")
  # beta.space --> data frame containing the slopes of temporal trends for the 
  # assemblage dissimilarity metrics per BBS route and ecoregion: 
    # - Total assemblage dissimilarity: beta.sor
    # - Species replacement: beta.sim
    # - Species loss: beta.sne
  beta.space <- read.csv("taxonomic_bdiv_ch_1990_2019_nomore1missing_data_July2024.csv")
# Environmental predictors data:
  # elev --> elevation data per BBS route and ecoregion
  elev <- unique(beta.space[, 1:5]) # get elevation data
  # lulc_lm --> data frame containing the slopes of temporal trends for the 
  # different LULC variables per ecoregion at resolutions 12.5 km, 25 km and 50 km.
  lulc_lm <- read.csv('lulc_vars_lm_tempTrend_1990_2019.csv')
  # clim_lm --> data frame containing the slopes of temporal trends for the 
  # different climatic variables per ecoregion at resolutions 12.5 km, 25 km and 50 km.
  clim_lm <- read.csv('clim_vars_lm_tempTrend_1990_2019_July2024.csv')

  # Join environmental data together
  env_lm  <- full_join(clim_lm, lulc_lm)
  env_lm  <- full_join(elev, env_lm)
  rm(clim_lm, lulc_lm, elev)
  
  # Mean proportion of LULC types per ecoregion between 1990-2019
  perc_ecor <- read.csv('lulc_vars_mean_prop_ecoregions_1990_2019.csv')
  perc_ecor <- perc_ecor %>% filter(! grepl('Others',lulc))
  names(perc_ecor)[names(perc_ecor) == 'NameL1_En'] <- 'Ecor_Lev1_name'
  perc_ecor$Ecor_Lev1_name <- gsub(" ", "_", perc_ecor$Ecor_Lev1_name)
# Split bio data by index and ecoregion into list of data frames
beta.space <- split(beta.space, list(beta.space$beta.index, beta.space$Ecor_Lev1_name)) 
# ------------------------------------------------------------------------------------------------------------------
# GAMMs
# ------------------------------------------------------------------------------------------------------------------
lulc_scen <- "A1B"
resoln <- c("12_5", "25", "50")
ecoregion <- c(unique(perc_ecor$Ecor_Lev1_name))
cor_data <- data.frame()
variables <- data.frame (predictor = c('Intercept','PPTmean May-July','Tmax May-July','Tmin May-July','Barren',
              'Crops Pasture','Forests','Grassland Shrubland','Urban Developed','Wetlands', 
              'Elevation (m)','AIC','cor str param','r squared','N','numb outl rm'))

i = 1 # For LULC scenario A1B
  print(lulc_scen[i])
  # Select the env variables for the i lulc scenario
  env_lm_i <- cbind(env_lm[ ,1:14],env_lm[,grepl(lulc_scen[i], colnames(env_lm))])
  # Select the % of lulc for the i lulc scenario
  perc_ecor_i <- perc_ecor %>% filter(grepl(lulc_scen[i],lulc))
  for (j in 1:length(resoln)) { # Resolution
    print(resoln[j])
    # Select the env variables for the k resolution
    env_lm_j <- cbind(env_lm_i[ ,1:5],env_lm_i[,grepl(resoln[j], colnames(env_lm_i))])
    for (k in 1:length(ecoregion)) { # Ecoregion
      print(ecoregion[k])
      # Select the data frames from the k ecoregion
      data_k <- beta.space[grepl(ecoregion[k], names(beta.space))]
      env_lm_k <- env_lm_j %>% filter(Ecor_Lev1_name == ecoregion[k])
      # Bind env data with each of the biodiversity data frames in list
      data_k <- lapply(data_k, function(x){x <- full_join(x,env_lm_k)})
      # lulc with the lowest % of coverage in ecoregion
      perc_ecor_j <- perc_ecor_i %>% filter(grepl(ecoregion[k],Ecor_Lev1_name))
      perc_ecor_j <- slice_min( perc_ecor_j, order_by =  perc_ecor_j[ ,1], n = 1, with_ties = F)
      perc_ecor_j <- gsub(paste("_", lulc_scen[i], "_", sep = ''), "", perc_ecor_j$lulc)
      # Create the list of predictors to include in the model
      use <- names(env_lm_k)[c(4,6:12,14:length(names(env_lm_k)),4)]
      # Eliminate lulc with the lowest % of coverage in ecoregion  
      use <- unique(as.character(na.omit(gsub(paste('\\S+',perc_ecor_j ,'\\S+', sep = ''), NA, use)))) 
      if(use[1]== 'na_elevation_m'){
       use <- use[c(2:length(use),1)]
      }
      use <- data.frame(use = use)
      # Eliminate highly correlated variables from list of predictors to use in model
      correl_vars <- env_lm_k [,c(4,6:15)]%>% as.matrix %>%
        cor %>% as.data.frame %>% rownames_to_column(var = 'var1') %>% gather(var2, value, -var1)
      correl_vars <- correl_vars[as.character(correl_vars$var1) != as.character(correl_vars$var2),]
      correl_vars <- correl_vars[!duplicated(apply(correl_vars,1,function(x) paste(sort(x),collapse=''))),]
      correl_vars <- correl_vars %>% filter (abs(value) > 0.6 )
      # Decide which one of the highly correlated variables will be kept in model 
      # and which to exclude. Suggestion: use if(){}else{} statements
      correl_vars$lulc_scen <- lulc_scen[i]
      correl_vars$resoln    <- resoln[j]
      correl_vars$ecoregion <- ecoregion[k]
      cor_data <- unique(rbind(cor_data,correl_vars))
      use <- use$use
      setwd("~/Desktop/Plos_ONE/Revisions_R1/R_Scrpits_to_share")
      # Run the GAMMs
      gamm_results    <- lapply(data_k, env_gamm, use)
      # Save results
      lapply(names(gamm_results), 
             function(nm){write.csv(gamm_results[[nm]], 
                                    paste("gamm_results",nm,ecoregion[k], lulc_scen[i], outlier0 [m],resoln[j],'km_july2024.csv', sep = '_'))})
    }
  }
# Save predictors correlation data   
write.csv(cor_data, 'beta_div_env_preds_correlation_vals_july2024.csv', row.names = F)

# ------------------------------------------------------------------------------------------------------------------
# Variations of the script to run the GAMM for species richness (alpha diversity)
# ------------------------------------------------------------------------------------------------------------------
# In line 18, use 
form <- as.formula(paste0("estimate~", paste0("s(",use,")", collapse = " + ")))

# In line 43, use
model_summary$index     <- "alpha.div"

# In line 86, use
alpha.time <- read.csv("taxonomic_alpha_itself_ch_1990_2019_nomore1missing_July2024.csv")
# alpha.time --> data frame containing the slopes of temporal trends for the 
# assemblage species richness per BBS route and ecoregion:

# Do not run line 108, there is only one metric for species richness

# In line 134, use
data_k <- alpha.time %>% filter(Ecor_Lev1_name == ecoregion[k])

# In line 137, use
data_k <- full_join(data_k,env_lm_k)

# In line 165, use
gamm_results    <- env_gamm(data_k, use)

# In lines 169 and 173, change the name of files to reflect these are results
# for species richness models
