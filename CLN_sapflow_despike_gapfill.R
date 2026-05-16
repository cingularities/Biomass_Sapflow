#######################
# TREE & SPECIES GAP FILL
#######################
#Written by Cindy Norton 2026
#Script segments a point cloud into clusters using watershed segmentation on rasterized point cloud
#install.packages('VoxR')
#remotes::install_github('bi0m3trics/spanner')
#install.packages("remotes")
#install.packages(c('spanner','readxl',"dplR","dplyr","tidyr", "geosphere", "lidR", "raster", "TreeLS", "rgdal", "rgeos", "sp","sf","tibble","ggplot2","tidyverse"),dependencies=TRUE)
#remotes::install_github("Jean-Romain/lidR", dependencies=TRUE)
#remotes::install_github('tiagodc/TreeLS',dependencies=TRUE)
#remotes::install_github('cszang/treeclim',dependencies=TRUE)
#install.packages("janitor")
gc()
#Packages <- c("dplyr", "ggplot2",  "httr", "rjson", "splitstackshape", "jsonlite", "curl","purrr","reshape2","tidyr","stringr","broom","modelr","lubridate")
#lapply(Packages, library, character.only = TRUE)
#library("pracma")
Packages <- c('cowplot','matrixStats','patchwork','treeclim','remotes','readxl',"dplR","dplyr","tidyr", "geosphere", "lidR", 'terra',"sf","tibble","ggplot2","tidyverse","lubridate","janitor")
lapply(Packages, library, character.only = TRUE)

gc()
setwd("//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/")

options(scipen=999)
#Written by Cindy Norton 2022
#install.packages('car')
library('ggplot2')
library(dplyr)
library(car)
library(reshape2)
library(gridExtra)
library(tidyverse)
library(mgcv)   # for GAMs
setwd('//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/outputs/')


thermal_data_dry <- readxl::read_xlsx('//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/Thermal/Results/Figures_Metadata_sapwood_sampling_080222.xlsx', sheet = 1) %>% 
  filter(weather == 'dry') %>%
  mutate(
    # Ensure key columns are numeric
    across(c(Core_radius_cm, Sapwood_width_visual_cm, Sapwood_width_thermal_cm), as.numeric))%>%
  mutate(
    Heartwood_radius_thermal = as.numeric(Core_radius_cm) - as.numeric(Sapwood_width_thermal_cm),
    Heartwood_radius_visual = as.numeric(Core_radius_cm) - as.numeric(Sapwood_width_visual_cm),
    SH_ratio_thermal = as.numeric(Sapwood_width_thermal_cm) / Heartwood_radius_thermal,
    SH_ratio_visual = as.numeric(Sapwood_width_visual_cm) / Heartwood_radius_visual,
    Sapwood_area_thermal = pi * as.numeric(Core_radius_cm)^2 - pi * Heartwood_radius_thermal^2,
    Sapwood_area_visual = pi * as.numeric(Core_radius_cm)^2 - pi * Heartwood_radius_visual^2,
    Heartwood_area_thermal = pi * Heartwood_radius_thermal^2,
    Heartwood_area_visual = pi * Heartwood_radius_visual^2
  )%>%
  rename(DBH = DBH_cm)%>%
  group_by(Species, Tree_ID) %>%
  filter(!is.na(DBH) & !is.na(Sapwood_area_thermal) & DBH > 0 & Sapwood_area_thermal > 0)
# group by Species and TreeID
options("scipen" = 100, "digits" = 2)













library(dplyr)
library(broom) # for tidy regression results

# Fit log-log model per Species
sapwood_constants <- thermal_data_dry %>%
  group_by(Species) %>%
  do({
    fit <- lm(log(Sapwood_area_thermal) ~ log(DBH), data = .)
    tibble(
      Species = unique(.$Species),
      B1 = coef(fit)[2],                  # slope
      B2 = exp(coef(fit)[1]),              # intercept back-transformed
      R2 = summary(fit)$r.squared,
      n = nrow(.)
    )
  })

sapwood_constants




############################################ Tree Ring and Ground Fusion ####################################################
#SET the path for folder where you have the file stoblue, each function will load the data format
groundA =   readxl::read_excel('//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Metadata/MtBigelow_Metadata/metadata_MtBigelow_Spring_2022.xlsx', sheet=1)
groundB =   readxl::read_excel('//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Metadata/MtBigelow_Metadata/metadata_MtBigelow_Spring_2022.xlsx', sheet=2)
groundC =   readxl::read_excel('//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Metadata/MtBigelow_Metadata/metadata_MtBigelow_Spring_2022.xlsx', sheet=3)
MBA<- "//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Treerings/MtBigelow_Treerings/MBA_revised.rwl"
MBB<- "//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Treerings/MtBigelow_Treerings/MBB.rwl"
MBC<- "//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Treerings/MtBigelow_Treerings/MBC.rwl"

mtB_A_xyz <- "//snre-snow/projects/Babst_Lidar_treering_CLNF/CLNorton/outputs/TLS_Bigelow_A_PCGR.xyz"
mtB_C_xyz <- "//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/outputs/TLS_Bigelow_C_GRPC.xyz"

drone <-  "//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/outputs/Mt.Bigelow_lidarFlight/Mt_Bigelow_2022_Final.las"
climate <- "//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/outputs/Mt.bigelow_prism_metereological.xlsx"
#noaa_climate_csv <- read.csv("https://www.ncei.noaa.gov/access/services/data/v1?dataset=daily-summaries&dataTypes=PRCP,TMAX,TMIN&stations=USC00025737,USC00025734,USC00025735,USC00025733,USC00025732,USC00026202&startDate=1950-01-01&endDate=2022-12-19&format=csv&options=includeStationName:1&includeStationName=true&units=metric")
#rwl<-dplR::read.rwl(tree) #read file path of rwl


###FUNCTION 1 - tree ring data parse and formatting###
#raw tree ring data parsing and mean tree cores
treering_parse <- function(tree) {
  rwl<-dplR::read.rwl(tree) #read file path of rwl
  detrend.rwi <- dplR::detrend(rwl = rwl, method = "Spline") %>% #detrend rwl
    t() %>% #transpose years as columns and row as each core sample 
    as.data.frame() %>% #as data frame
    dplyr::rename_with( ~ paste0(.x)) %>% #adding "year_" to each year column name for unique identifiers
    tibble::rownames_to_column()%>% #converts sample row name into a column to be parsed
    tidyr::separate(rowname,
                    into = c("site", "coreID"),
                    sep = "(?<=[A-ZA-Z])(?=[0-9])") %>% ## separates column of site and treeID
    tidyr::separate(coreID,
                    into = c("treeID", "core"),
                    sep = "(?<=[0-9])(?=[A-ZA-Za-z])")
  return(detrend.rwi)}
MBA_parse <- treering_parse(MBA)
MBB_parse <- treering_parse(MBB)
MBC_parse <- treering_parse(MBC)

#write.csv(data.transpose, file = "datatranspose.csv")
















##RAW GROUND DATA
MBA_coor <-c(-110.72624,32.41585) 	
MBB_coor <-c(-110.72578,32.41596) 	
MBC_coor <-c(-110.72517,32.41548) 	


###FUNCTION 2 - ground data azimuth and distance to coordinates###
ground_parse <- function(coor, ground) { #read xlsx ground metadata file
  dist_azim <- ground
  coor_list <- list() #creates empty list for new coordinates
  
  for (i in 1:nrow(dist_azim)) { #for loop to create a lat long for each azimuth and distance from plot center lat and lon
    stemAzimuth <- dist_azim$Azimuth[i] #grab azimuth
    stemDistance <- dist_azim$Distance[i] #grab distance
    
    new_coordinates <- geosphere::destPoint(coor, stemAzimuth,  stemDistance) %>% #estimates new lat lon coordinate
      as.data.frame() %>% #make as data frame
      add_column(treeID = dist_azim$treeID[i]) %>% #add treeID column of that tree
      add_column(stemAzimuth = stemAzimuth) %>% #add Azimuth column of that tree
      add_column(stemDistance = stemDistance) %>% #add Distance column of that tree
      add_column(Species = dist_azim$Species[i]) %>% #add Species column of that tree
      add_column(Height = dist_azim$Height[i])%>% #add Height column of that tree
      add_column(DBH = dist_azim$DBH[i]) #add DBH column of that tree
    
    
    coor_list[[i]] <- new_coordinates #put result of loop in list
    
  }
  
  coor_df<- coor_list %>% bind_rows() %>% na.omit() #remove NA rows
  
  
  return(coor_df)
}

MBA_ground <- ground_parse(MBA_coor, groundA)%>%
  mutate(Species = case_when(
    Species == "ponderosa"    ~ "PIPO",
    Species == "strobiformis" ~ "PISF",
    Species == "menziesii"    ~ "PSME",
    TRUE                      ~ Species   # keep other Species unchanged
  ))

MBB_ground <- ground_parse(MBB_coor, groundB) %>%
  mutate(Species = case_when(
    Species == "ponderosa"    ~ "PIPO",
    Species == "strobiformis" ~ "PISF",
    Species == "menziesii"    ~ "PSME",
    TRUE                      ~ Species   # keep other Species unchanged
  ))

MBC_ground <- ground_parse(MBC_coor, groundC)  %>%
  mutate(Species = case_when(
    Species == "ponderosa"    ~ "PIPO",
    Species == "strobiformis" ~ "PISF",
    Species == "menziesii"    ~ "PSME",
    TRUE                      ~ Species   # keep other Species unchanged
  ))





all_ground<- MBA_ground%>%bind_rows(MBB_ground)%>%bind_rows(MBC_ground) 























































##FUNCTION 3 -  tree ring and ground data merge###
###TREE rings and ground merge#####
treering_ground <- function (treering_parse, ground_parse, plot) {
  
  treering_ground_join <- treering_parse %>% 
    subset(site == plot) %>% #subsets to input plot
    lapply(as.numeric)%>% #makes numeric
    as.data.frame() %>% #makes into df
    group_by(treeID) %>% #groups data by the treeID
    summarise(across(-c(site, core), mean, na.rm = TRUE))%>% #a mean summary of multiple cores per tree
    left_join(ground_parse, by = "treeID") 
  #joins tree ring data with ground metadat using treeID
  colnames(treering_ground_join)<-gsub("X","",colnames(treering_ground_join))
  
  return(treering_ground_join)
}

MBA_treering_ground <- treering_ground(MBA_parse, MBA_ground, "MBA")
MBB_treering_ground <- treering_ground(MBB_parse, MBB_ground, "MBB")
MBC_treering_ground <- treering_ground(MBC_parse, MBC_ground, "MBC")



# Define DBH values per Species
dbh_list <- list(
  PIPO = c(22.3, 29.2, 43.5, 27.1),
  PISF = c(31.9, 33.4, 49.8),
  PSME = c(20, 37.1)
)


# Tolerance (e.g., 1 cm)
tol <- 10

# Subset trees with DBH within tolerance
MBA_treering_ground <- MBA_treering_ground %>%
  filter(
    (Species == "PIPO" & sapply(DBH, function(x) any(abs(x - dbh_list$PIPO) <= tol))) |
      (Species == "PISF" & sapply(DBH, function(x) any(abs(x - dbh_list$PISF) <= tol))) |
      (Species == "PSME" & sapply(DBH, function(x) any(abs(x - dbh_list$PSME) <= tol)))
  )

MBA_treering_ground



# Subset trees with DBH within tolerance
MBB_treering_ground <- MBB_treering_ground %>%
  filter(
    (Species == "PIPO" & sapply(DBH, function(x) any(abs(x - dbh_list$PIPO) <= tol))) |
      (Species == "PISF" & sapply(DBH, function(x) any(abs(x - dbh_list$PISF) <= tol))) |
      (Species == "PSME" & sapply(DBH, function(x) any(abs(x - dbh_list$PSME) <= tol)))
  )

MBB_treering_ground





# Subset trees with DBH within tolerance
MBC_treering_ground <- MBC_treering_ground %>%
  filter(
    (Species == "PIPO" & sapply(DBH, function(x) any(abs(x - dbh_list$PIPO) <= tol))) |
      (Species == "PISF" & sapply(DBH, function(x) any(abs(x - dbh_list$PISF) <= tol))) |
      (Species == "PSME" & sapply(DBH, function(x) any(abs(x - dbh_list$PSME) <= tol)))
  )

MBC_treering_ground

















###FUNCTION 4 - DBH and Percentage reconstruction###
###DBH RECONSTRUCTIONS
DBH_recon <- function(tree,ground_parse){
  
  rwl<-dplR::read.rwl(tree)%>%#read file path of rwl
    t() %>% #transpose years as columns and row as each core sample 
    as.data.frame() %>% #as data frame
    dplyr::rename_with( ~ paste0("year_",.x)) %>% #adding "year_" to each year column name for unique identifiers
    tibble::rownames_to_column()%>% #converts sample row name into a column to be parsed
    tidyr::separate(rowname,
                    into = c("site", "coreID"),
                    sep = "(?<=[A-ZA-Z])(?=[0-9])") %>% ## separates column of site and treeID
    tidyr::separate(coreID,
                    into = c("treeID", "core"),
                    sep = "(?<=[0-9])(?=[A-ZA-Za-z])")%>%
    rev()
  
  treering_tree <- rwl %>% 
    lapply(as.numeric)%>% #makes numeric
    as.data.frame() %>% #makes into df
    group_by(treeID) %>% #groups data by the treeID
    summarise(across(-c(site, core), mean, na.rm = TRUE))%>%
    left_join(ground_parse, by = "treeID") %>%
    as.data.frame()
  
  treering_tree[is.na(treering_tree)] <- 0  #zero all NA
  year_cols <- which(substr(x = colnames(treering_tree), start = 1, stop = 4) == "year")  # Identify the columns that have data of interest which are year columns
  treering_tree$sum <- rowSums(x = treering_tree[, year_cols], na.rm = TRUE)  #row sums all the years into column, removing NA 
  csum <- treering_tree  # Make a copy of that data frame that we will use for difference
  csum[, year_cols] <- sapply(1:ncol(csum[, year_cols]), function(col){rowSums(csum[, year_cols][1:col], na.rm = TRUE)}) #convert years to cumulative sum of year columns
  dif <- csum # Make a copy of that data frame that we will use for difference
  dif[, year_cols] <- (dif$sum) - (dif[, year_cols]) #calculate difference
  frac <- csum # Make a copy of that data frame that we will use for difference
  frac[, year_cols] <- dif[, year_cols]/dif$sum #calculate difference
  dbh <- csum # Make a copy of that data frame that we will use for difference
  dbh[, year_cols] <- as.numeric(frac$DBH)*frac[, year_cols] #calculate difference
  
  # Reality check to see that percentages add up to 100
  
  return(dbh)
}







MBA_DBH_recon_results <- DBH_recon(MBA, MBA_ground)%>%
  mutate(Species = case_when(
    Species == "ponderosa"    ~ "PIPO",
    Species == "strobiformis" ~ "PISF",
    Species == "menziesii"    ~ "PSME",
    TRUE                      ~ Species   # keep other Species unchanged
  ))
MBB_DBH_recon_results <- DBH_recon(MBB, MBB_ground)%>%
  mutate(Species = case_when(
    Species == "ponderosa"    ~ "PIPO",
    Species == "strobiformis" ~ "PISF",
    Species == "menziesii"    ~ "PSME",
    TRUE                      ~ Species   # keep other Species unchanged
  ))
MBC_DBH_recon_results <- DBH_recon(MBC, MBC_ground)%>%
  mutate(Species = case_when(
    Species == "ponderosa"    ~ "PIPO",
    Species == "strobiformis" ~ "PISF",
    Species == "menziesii"    ~ "PSME",
    TRUE                      ~ Species   # keep other Species unchanged
  ))

library(dplyr)
library(tidyr)
library(stringr)
library(zoo)
library(lubridate)
library(xgboost)
library(ggplot2)

# ----------------------------
# Load sapflow
# ----------------------------
sapflow <- read.csv('//snow/Projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/sapflow/MB_SapFlow_2010-2024.csv')
sapflow[sapflow == "#DIV/0!"] <- NA

# ----------------------------
# Load Ameriflow tower data
# ----------------------------
fluxtower <- read.csv("//snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/sapflow/AMF_US-MtB_FLUXNET_FULLSET_2009-2024_3-7/AMF_US-MtB_FLUXNET_FULLSET_HH_2009-2024_3-7.csv")

fluxtower_edit <- fluxtower %>%
  mutate(
    TIMESTAMP_START = ymd_hm(TIMESTAMP_START),
    YEAR  = year(TIMESTAMP_START),
    DOY   = yday(TIMESTAMP_START),
    DAY   = day(TIMESTAMP_START),
    TOD   = hour(TIMESTAMP_START) + minute(TIMESTAMP_START)/60
  ) %>%
  select(YEAR, DOY, DAY, TOD, SW_IN_F, TA_F, P_F, VPD_F, PA_F, LW_IN_F)

# ----------------------------
# Helper functions
# ----------------------------
get_species <- function(x) str_extract(x, "^[A-Za-z]+")
get_tree_id <- function(x) str_extract(x, "^[A-Za-z]+[0-9]+")
get_orient  <- function(x) str_extract(x, "[A-Za-z]$")

# ----------------------------
# STEP 0: Convert DOY + YEAR to date
# ----------------------------
sapflow <- sapflow %>%
  mutate(date = as.Date(DOY - 1, origin = paste0(YEAR, "-01-01")))

# ----------------------------
# STEP 1: Identify SENSOR columns (numeric, exclude YEAR/DOY/TOD)
# ----------------------------
sensor_cols <- sapflow %>%
  select(where(is.numeric)) %>%
  select(-YEAR, -DOY, -TOD) %>%
  colnames()

# ----------------------------
# STEP 2: Rule-based despiking
# ----------------------------
despike_sapflow <- function(x, k = 3, mad_thresh = 3, max_rate = 5) {
  med <- zoo::rollapply(x, k, median, na.rm = TRUE,
                        fill = NA, align = "center")
  
  madv <- zoo::rollapply(
    x, k,
    function(z) mad(z, constant = 1, na.rm = TRUE),
    fill = NA, align = "center"
  )
  
  spike_median <- abs(x - med) > (mad_thresh * madv)
  spike_rate   <- c(NA, abs(diff(x))) > max_rate
  
  x[spike_median | spike_rate] <- NA
  return(x)
}

sapflow_despiked <- sapflow %>%
  mutate(across(all_of(sensor_cols), ~ despike_sapflow(.x)))

# ----------------------------
# STEP 3: Sensor metadata
# ----------------------------
sensor_meta <- tibble(sensor = sensor_cols) %>%
  mutate(
    species = get_species(sensor),
    tree_id = get_tree_id(sensor),
    orient  = get_orient(sensor)
  )


library(dplyr)
library(xgboost)
library(tidyr)

# ----------------------------
# XGBoost gap-filling function
# ----------------------------
ml_gapfill_xgb <- function(sf_col, predictors_df,
                           nrounds = 300,
                           seed = 42) {
  
  set.seed(seed)
  
  train_idx <- !is.na(sf_col) & complete.cases(predictors_df)
  if (sum(train_idx) < 50) return(sf_col)
  
  dtrain <- xgb.DMatrix(
    data  = as.matrix(predictors_df[train_idx, ]),
    label = sf_col[train_idx]
  )
  
  params <- list(
    objective        = "reg:squarederror",
    eval_metric      = "rmse",
    max_depth        = 4,
    eta              = 0.05,
    subsample        = 0.8,
    colsample_bytree = 0.8
  )
  
  xgb_model <- xgb.train(
    params  = params,
    data    = dtrain,
    nrounds = nrounds,
    verbose = 0
  )
  
  # Predict where predictors are complete
  pred_idx <- complete.cases(predictors_df)
  sf_pred  <- rep(NA, length(sf_col))
  
  if (sum(pred_idx) > 0) {
    dpred <- xgb.DMatrix(
      data = as.matrix(predictors_df[pred_idx, ])
    )
    sf_pred[pred_idx] <- predict(xgb_model, dpred)
  }
  
  # Fill ONLY gaps
  fill_idx <- is.na(sf_col) & pred_idx
  sf_col[fill_idx] <- sf_pred[fill_idx]
  
  # Physical constraint: sap flow cannot be negative
  sf_col[sf_col < 0] <- 0
  
  return(sf_col)
}

# ----------------------------
# Extract species from sensor name
# ----------------------------
get_species <- function(sensor_name) {
  gsub("[0-9]+[ns]?$", "", sensor_name)  # e.g., "PIPO1n" -> "PIPO"
}

# ----------------------------
# Prepare predictor matrix
# ----------------------------
# Align sapflow with tower predictors
sapflow_despiked <- sapflow_despiked %>%
  left_join(fluxtower_edit, by = c("YEAR", "DOY", "TOD"))

# Now predictors exist
predictor_cols <- c("SW_IN_F", "TA_F", "P_F", "VPD_F", "PA_F", "LW_IN_F")
predictors <- sapflow_despiked %>%
  select(all_of(predictor_cols))

# ----------------------------
# Species-level XGBoost gap-filling
# ----------------------------
sapflow_filled <- sapflow_despiked

for (col in sensor_cols) {
  
  species <- get_species(col)
  
  # Other sensors of same species, exclude target
  other_sensors <- sensor_cols[ sensor_cols != col & grepl(species, sensor_cols) ]
  
  # Compute species-level mean sapflow
  if (length(other_sensors) >= 1) {
    species_mean_sf <- rowMeans(
      sapflow_filled[, other_sensors, drop = FALSE],
      na.rm = TRUE
    )
  } else {
    species_mean_sf <- sapflow_filled[[col]]
  }
  
  # Combine predictors with species-level signal
  predictors_species <- predictors %>%
    mutate(species_mean_sf = species_mean_sf)
  
  # Apply XGBoost gap-filling
  sapflow_filled[[col]] <- ml_gapfill_xgb(
    sf_col = sapflow_filled[[col]],
    predictors_df = predictors_species
  )
}

# ----------------------------
# Apply nighttime zero-flow constraint
# ----------------------------
sapflow_filled <- sapflow_filled %>%
  mutate(across(
    all_of(sensor_cols),
    ~ ifelse(SW_IN_F <= 5, 0, .x)
  ))

# ----------------------------
# Save final species-level gap-filled dataset
# ----------------------------
write.csv(
  sapflow_filled,
  "sapflux_final_species_xgb.csv",
  row.names = FALSE
)


















library(dplyr)
library(tidyr)

# Make a long-format table for raw vs XGBoost-filled
raw_vs_xgb <- sapflow %>%
  select(YEAR, DOY, all_of(sensor_cols)) %>%
  pivot_longer(
    cols = all_of(sensor_cols),
    names_to = "sensor",
    values_to = "raw"
  ) %>%
  left_join(
    sapflow_filled %>%
      select(YEAR, DOY, all_of(sensor_cols)) %>%
      pivot_longer(
        cols = all_of(sensor_cols),
        names_to = "sensor",
        values_to = "xgb_filled"
      ),
    by = c("YEAR", "DOY", "sensor")
  )

# Compute MAE per sensor per year
mae_by_sensor_year <- raw_vs_xgb %>%
  group_by(YEAR, sensor) %>%
  summarise(
    MAE = mean(abs(raw - xgb_filled), na.rm = TRUE),
    .groups = "drop"
  )

# View results
head(mae_by_sensor_year)



library(ggplot2)
windows()
ggplot(mae_by_sensor_year, aes(x = factor(YEAR), y = MAE, fill = sensor)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    title = "Mean Absolute Error (MAE) per Year per Sensor — XGBoost Gap-filling",
    x = "Year",
    y = "MAE (Sap Flow Units)",
    fill = "Sensor"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )















library(dplyr)
library(tidyr)

# Long format: raw vs XGBoost-filled
raw_vs_xgb <- sapflow %>%
  select(YEAR, DOY, all_of(sensor_cols)) %>%
  pivot_longer(
    cols = all_of(sensor_cols),
    names_to = "sensor",
    values_to = "raw"
  ) %>%
  left_join(
    sapflow_filled %>%
      select(YEAR, DOY, all_of(sensor_cols)) %>%
      pivot_longer(
        cols = all_of(sensor_cols),
        names_to = "sensor",
        values_to = "xgb_filled"
      ),
    by = c("YEAR", "DOY", "sensor")
  )

# Compute RMSE per sensor per year
rmse_by_sensor_year <- raw_vs_xgb %>%
  group_by(YEAR, sensor) %>%
  summarise(
    RMSE = sqrt(mean((raw - xgb_filled)^2, na.rm = TRUE)),
    .groups = "drop"
  )

# View results
head(rmse_by_sensor_year)




library(ggplot2)
windows()
ggplot(rmse_by_sensor_year, aes(x = factor(YEAR), y = RMSE, fill = sensor)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    title = "RMSE per Year per Sensor — XGBoost Gap-filling",
    x = "Year",
    y = "RMSE (Sap Flow Units)",
    fill = "Sensor"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )



















sapflow_fluxtower <- sapflow %>%
  left_join(fluxtower_edit, by = c("YEAR", "DOY", "TOD"))



# ----------------------------
# Define sapflow sensor columns
# ----------------------------
cols_to_plot <- sapflow_fluxtower %>%
  select(where(is.numeric)) %>%          # all numeric columns
  select(-YEAR, -DOY, -TOD,             # remove time columns
         -SW_IN_F, -TA_F, -P_F, -VPD_F, -PA_F, -LW_IN_F,-DAY) %>%  # remove predictors
  colnames()


# ----------------------------
# Optional: correlation per sensor vs environmental predictors
# ----------------------------
df_corr <- sapflow_filled %>%
  select(YEAR, DOY, TOD, all_of(sensor_cols), all_of(predictor_cols)) %>%
  pivot_longer(cols = all_of(sensor_cols),
               names_to = "sensor",
               values_to = "sapflow")

cor_test_df <- function(x, y) {
  test <- cor.test(x, y, method = "pearson")
  tibble(r = test$estimate, p = test$p.value,
         sig = ifelse(test$p.value < 0.05, "*", ""))
}

cor_results <- df_corr %>%
  group_by(sensor) %>%
  summarize(across(all_of(predictor_cols), ~ cor_test_df(.x, sapflow), .names = "{.col}"),
            .groups = "drop") %>%
  pivot_longer(-sensor, names_to = "predictor", values_to = "res") %>%
  unnest(res) %>%
  mutate(label = paste0(round(r, 2), sig))

ggplot(cor_results, aes(x = predictor, y = sensor, fill = r)) +
  geom_tile() +
  geom_text(aes(label = label), color = "black", size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Pearson Correlation (XGBoost-filled Sap Flow vs Predictors)",
       x = "Predictor", y = "Sensor", fill = "r") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
























library(dplyr)
library(tidyr)
library(ggplot2)

# ----------------------------
# Define columns
# ----------------------------
base_cols <- c("date", "YEAR")
sensor_cols <- cols_to_plot  # numeric sapflow sensors in sapflow_filled

# ----------------------------
# Prepare long-format data for plotting
# ----------------------------
raw_long <- sapflow %>%
  select(all_of(base_cols), all_of(sensor_cols)) %>%
  pivot_longer(
    cols = all_of(sensor_cols),
    names_to = "sensor",
    values_to = "sapflow"
  ) %>%
  mutate(method = "Raw")

rule_long <- sapflow_despiked %>%
  select(all_of(base_cols), all_of(sensor_cols)) %>%
  pivot_longer(
    cols = all_of(sensor_cols),
    names_to = "sensor",
    values_to = "sapflow"
  ) %>%
  mutate(method = "Rule-based")

xgb_long <- sapflow_filled %>%
  select(all_of(base_cols), all_of(sensor_cols)) %>%
  pivot_longer(
    cols = all_of(sensor_cols),
    names_to = "sensor",
    values_to = "sapflow"
  ) %>%
  mutate(method = "XGBoost-filled")

# ----------------------------
# Custom color palette
# ----------------------------
method_colors <- c(
  "Raw" = "#1f78b4",          # blue
  "Rule-based" = "#e31a1c",   # green
  "XGBoost-filled" = "#e31a1c" # red
)

# ----------------------------
# 1️⃣ Plot: Raw vs. Rule-based despiked
# ----------------------------
raw_vs_rule <- bind_rows(raw_long, rule_long)

plot_sensor <- "PISF1n"  # change to sensor of interest
windows(width = 14, height = 10)

ggplot(
  raw_vs_rule %>% filter(sensor == plot_sensor),
  aes(x = date, y = sapflow, color = method)
) +
  geom_line(alpha = 0.8, linewidth = 0.9) +
  facet_wrap(~ YEAR, scales = "free_x") +
  scale_color_manual(values = method_colors[c("Raw", "Rule-based")]) +
  labs(
    title = paste("Raw vs. Rule-based Despiking —", plot_sensor),
    x = "Date",
    y = "Sap Flow",
    color = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )

# ----------------------------
# 2️⃣ Plot: Raw vs. XGBoost-filled
# ----------------------------
raw_vs_xgb <- bind_rows(raw_long, xgb_long)

plot_sensor <- "PISF1n"  # same or another sensor
windows(width = 14, height = 10)

ggplot(
  raw_vs_xgb %>% filter(sensor == plot_sensor),
  aes(x = date, y = sapflow, color = method)
) +
  geom_line(alpha = 0.8, linewidth = 0.9) +
  facet_wrap(~ YEAR, scales = "free_x") +
  scale_color_manual(values = method_colors[c("Raw", "XGBoost-filled")]) +
  labs(
    title = paste("Raw vs. XGBoost-filled (Gap-filled) —", plot_sensor),
    x = "Date",
    y = "Sap Flow",
    color = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )
# ============================================================
# XGBoost-filled PISF1n and SW_IN_F
# ============================================================

# ============================================================
# Dual-axis plot: XGBoost-filled PISF1n and SW_IN_F
# ============================================================

library(dplyr)
library(ggplot2)

# ----------------------------
# Compute scaling factor
# ----------------------------
scale_factor <- max(sapflow_filled$PISF1n, na.rm = TRUE) /
  max(sapflow_filled$TA_F, na.rm = TRUE)

# ----------------------------
# Plot
# ----------------------------
windows(width = 14, height = 10)

ggplot(sapflow_filled, aes(x = date)) +
  
  # PISF1n (left axis)
  geom_line(
    aes(y = PISF1n, color = "PISF1n"),
    linewidth = 0.9,
    alpha = 0.85
  ) +
  
  # SW_IN_F (scaled, right axis)
  geom_line(
    aes(y = TA_F * scale_factor, color = "TA_F"),
    linewidth = 0.9,
    alpha = 0.85
  ) +
  
  facet_wrap(~ YEAR, scales = "free_x") +
  
  scale_y_continuous(
    name = "Sap Flow (PISF1n)",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = expression(SW~(W~m^{-2}))
    )
  ) +
  
  scale_color_manual(
    values = c(
      "PISF1n"  = "#e31a1c",
      "TA_F" = "#1f78b4"
    )
  ) +
  
  labs(
    title = "XGBoost Gap-filled Time Series: PISF1n and TA_F",
    x = "Date",
    color = "Variable"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title.y.left  = element_text(size = 16, color = "#e31a1c"),
    axis.title.y.right = element_text(size = 16, color = "#1f78b4"),
    axis.text          = element_text(size = 12),
    strip.text         = element_text(size = 14, face = "bold"),
    legend.title       = element_text(size = 14),
    legend.text        = element_text(size = 12),
    legend.position    = "bottom",
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )



# ----------------------------
# Define year to plot
# ----------------------------
plot_year <- 2024  # change to desired year

# ----------------------------
# 1️⃣ Raw vs. Rule-based: all sensors for a specific year
# ----------------------------
raw_vs_rule_year <- raw_vs_rule %>%
  filter(YEAR == plot_year)

windows(width = 16, height = 12)
ggplot(raw_vs_rule_year, aes(x = date, y = sapflow, color = method)) +
  geom_line(alpha = 0.9, linewidth = 1) +  # thicker for clarity
  facet_wrap(~ sensor, scales = "free_y") +
  scale_color_manual(values = method_colors[c("Raw", "Rule-based")]) +  # distinct colors
  labs(
    title = paste("Raw vs. Rule-based Despiking —", plot_year),
    x = "Date",
    y = "Sap Flow",
    color = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )

windows()
ggplot(raw_vs_rule_year, aes(x = date, y = sapflow, color = method)) +
  geom_line(alpha = 0.7, linewidth = 0.8) +               # original lines
  facet_wrap(~ sensor, scales = "free_y") +               # separate panel per sensor
  scale_color_manual(values = method_colors[c("Raw", "Rule-based")]) +
  labs(
    title = paste("Raw vs. Rule-based (Despiked) —", plot_year),
    x = "Date",
    y = "Sap Flow",
    color = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )


# ----------------------------
# 2️⃣ Raw vs. XGBoost-filled: all sensors for a specific year
# ----------------------------
raw_vs_xgb_year <- raw_vs_xgb %>%
  filter(YEAR == plot_year)

windows(width = 16, height = 12)
ggplot(raw_vs_xgb_year, aes(x = date, y = sapflow, color = method)) +
  geom_line(alpha = 0.9, linewidth = 1) +  # thicker for clarity
  facet_wrap(~ sensor, scales = "free_y") +
  scale_color_manual(values = method_colors[c("Raw", "XGBoost-filled")]) +  # distinct colors
  labs(
    title = paste("Raw vs. XGBoost-filled (Gap-filled) —", plot_year),
    x = "Date",
    y = "Sap Flow",
    color = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )


windows()
ggplot(raw_vs_xgb_year, aes(x = date, y = sapflow, color = method)) +
  geom_line(alpha = 0.7, linewidth = 0.8) +               # original lines
  facet_wrap(~ sensor, scales = "free_y") +               # separate panel per sensor
  scale_color_manual(values = method_colors[c("Raw", "XGBoost-filled")]) +
  labs(
    title = paste("Raw vs. XGBoost-filled (Gap-filled) —", plot_year),
    x = "Date",
    y = "Sap Flow",
    color = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )






















# ----------------------------
# Define sapflow sensor columns
# ----------------------------
cols_to_plot <- sapflow_filled %>%
  select(where(is.numeric)) %>%          # all numeric columns
  select(-YEAR, -DOY, -TOD,             # remove time columns
         -SW_IN_F, -TA_F, -P_F, -VPD_F, -PA_F, -LW_IN_F,-DAY) %>%  # remove predictors
  colnames()



library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)

# ----------------------------
# Prepare monthly data
# ----------------------------
sapflux_monthly <- sapflow_filled %>%
  select(-TOD) %>%
  mutate(
    YEAR  = factor(YEAR),
    MONTH = month(date),
    MONTH_LABEL = month(date, label = TRUE, abbr = TRUE)
  ) %>%
  pivot_longer(
    cols = all_of(cols_to_plot),
    names_to = "sensor",
    values_to = "sapflow"
  ) %>%
  group_by(YEAR, MONTH, MONTH_LABEL, sensor) %>%
  summarize(
    mean_value = mean(sapflow, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# Plot: monthly mean sap flow per sensor, faceted by year
# ----------------------------
windows(width = 16, height = 12)
library(viridis)  # install.packages("viridis") if needed

ggplot(sapflux_monthly, aes(x = MONTH, y = mean_value, group = sensor, color = sensor)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  facet_wrap(~ YEAR, ncol = 4) +
  scale_x_continuous(
    breaks = 1:12,
    labels = levels(sapflux_monthly$MONTH_LABEL)
  ) +
  scale_color_viridis_d(option = "turbo") +  # <- changed color scheme
  labs(
    title = "Monthly Mean Sap Flow (All Sensors)",
    x = "Month",
    y = "Sap Flow",
    color = "Sensor"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 20, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "right",
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )






library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)

# ----------------------------
# Prepare monthly data
# ----------------------------
sapflux_monthly <- sapflow %>%
  select(-TOD) %>%
  mutate(
    YEAR  = factor(YEAR),
    MONTH = month(date),
    MONTH_LABEL = month(date, label = TRUE, abbr = TRUE)
  ) %>%
  pivot_longer(
    cols = all_of(cols_to_plot),
    names_to = "sensor",
    values_to = "sapflow"
  ) %>%
  group_by(YEAR, MONTH, MONTH_LABEL, sensor) %>%
  summarize(
    mean_value = mean(sapflow, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# Plot: monthly mean sap flow per sensor, faceted by year
# ----------------------------
windows(width = 16, height = 12)
library(viridis)  # install.packages("viridis") if needed

ggplot(sapflux_monthly, aes(x = MONTH, y = mean_value, group = sensor, color = sensor)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  facet_wrap(~ YEAR, ncol = 4) +
  scale_x_continuous(
    breaks = 1:12,
    labels = levels(sapflux_monthly$MONTH_LABEL)
  ) +
  scale_color_viridis_d(option = "turbo") +  # <- changed color scheme
  labs(
    title = "Monthly Mean Sap Flow (All Sensors)",
    x = "Month",
    y = "Sap Flow",
    color = "Sensor"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(size = 20, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = "right",
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )


















library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(viridis)

# ----------------------------
# Compute per-tree yearly mean sap flow
# ----------------------------
sapflow_yearly_tree <- sapflow_filled %>%
  select(YEAR, all_of(sensor_cols)) %>%
  pivot_longer(cols = all_of(sensor_cols),
               names_to = "sensor",
               values_to = "sapflow") %>%
  left_join(sensor_meta %>% select(sensor, tree_id), by = "sensor") %>%
  group_by(YEAR, tree_id) %>%
  summarize(mean_sapflow = mean(sapflow, na.rm = TRUE), .groups = "drop")

# ----------------------------
# Compute yearly mean of SW_IN_F and P_F
# ----------------------------
predictors_yearly <- fluxtower_edit %>%
  group_by(YEAR) %>%
  summarize(
    SW_IN_F = mean(SW_IN_F, na.rm = TRUE),
    P_F     = mean(P_F, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(SW_IN_F, P_F),
               names_to = "predictor",
               values_to = "value")

# ----------------------------
# Colors
# ----------------------------
tree_colors <- viridis(length(unique(sapflow_yearly_tree$tree_id)))
names(tree_colors) <- unique(sapflow_yearly_tree$tree_id)
predictor_colors <- c("SW_IN_F" = "#e31a1c", "P_F" = "#e31a1c")

# ----------------------------
# Separate plots for each predictor with correct secondary axis
# ----------------------------
plots <- lapply(c("SW_IN_F", "P_F"), function(pred) {
  
  pred_data <- predictors_yearly %>% filter(predictor == pred)
  
  sap_min  <- min(sapflow_yearly_tree$mean_sapflow, na.rm = TRUE)
  sap_max  <- max(sapflow_yearly_tree$mean_sapflow, na.rm = TRUE)
  pred_min <- min(pred_data$value, na.rm = TRUE)
  pred_max <- max(pred_data$value, na.rm = TRUE)
  
  pred_data <- pred_data %>%
    mutate(scaled_value = (value - pred_min) / (pred_max - pred_min) * (sap_max - sap_min) + sap_min)
  # ... inside the lapply loop for each predictor
  ggplot() +
    # Sap flow per tree colored
    geom_line(data = sapflow_yearly_tree,
              aes(x = YEAR, y = mean_sapflow, color = tree_id, group = tree_id),
              linewidth = 1.2) +
    geom_point(data = sapflow_yearly_tree,
               aes(x = YEAR, y = mean_sapflow, color = tree_id),
               size = 3) +
    # Predictor line dashed
    geom_line(data = pred_data,
              aes(x = YEAR, y = scaled_value),
              color = predictor_colors[pred], linetype = "dashed", linewidth = 1.2) +
    geom_point(data = pred_data,
               aes(x = YEAR, y = scaled_value),
               color = predictor_colors[pred], shape = 17, size = 3) +
    scale_color_viridis(discrete = TRUE, name = "Tree ID") +
    scale_y_continuous(
      name = "Mean Sap Flow",
      sec.axis = sec_axis(
        trans = ~ (. - sap_min) / (sap_max - sap_min) * (pred_max - pred_min) + pred_min,
        name = pred
      )
    ) +
    labs(title = paste("Per-Tree Yearly Mean Sap Flow with", pred),
         x = "Year") +
    theme_bw(base_size = 28) +
    theme(
      plot.title   = element_text(size = 28, face = "bold"),
      axis.title   = element_text(size = 22),
      axis.text    = element_text(size = 20),
      axis.text.x  = element_text(angle = 45, hjust = 1)
    )
  
})

# ----------------------------
# Display plots in 2x1 layout
# ----------------------------
library(gridExtra)
windows(width = 16, height = 12)
grid.arrange(plots[[1]], plots[[2]], ncol = 1)




















library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(ggpattern)  # for patterns in bars

# ----------------------------
# Compute seasonal mean per species
# ----------------------------
sapflow_season_species <- sapflow_filled %>%
  select(date, all_of(sensor_cols)) %>%
  pivot_longer(
    cols = all_of(sensor_cols),
    names_to = "sensor",
    values_to = "sapflow"
  ) %>%
  left_join(sensor_meta %>% select(sensor, species), by = "sensor") %>%
  mutate(
    month = month(date),
    season = case_when(
      month %in% 7:9  ~ "monsoon",
      month %in% 4:6  ~ "pre-monsoon",
      month == 10     ~ "fall",
      month %in% c(11,12,1,2,3) ~ "winter",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(species, season) %>%
  summarize(
    mean_sapflow = mean(sapflow, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# Define patterns for species
# ----------------------------
species_list <- unique(sapflow_season_species$species)
pattern_types <- c("stripe", "crosshatch", "circle", "wave", "diamond")  # extend if needed
species_patterns <- setNames(pattern_types[1:length(species_list)], species_list)

# ----------------------------
# Plot seasonal sap flow per species with black bars and patterns
# ----------------------------
windows(width = 16, height = 10)

ggplot(sapflow_season_species, aes(x = season, y = mean_sapflow, pattern = species)) +
  geom_col_pattern(
    fill = "black",              # black bars
    color = "black",             # bar border
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 0.9,
    pattern_density = 0.5,
    pattern_spacing = 0.05
  ) +
  geom_text(aes(label = round(mean_sapflow, 1)),
            position = position_dodge(width = 0.8),
            vjust = -0.5,
            size = 4,
            color = "black") +
  scale_pattern_manual(values = species_patterns) +
  labs(
    title = "Seasonal Mean Sap Flow by Species (All Years)",
    x = "Season",
    y = "Mean Sap Flow",
    pattern = "Species"
  ) +
  theme_bw(base_size = 28) +
  theme(
    plot.title   = element_text(size = 28, face = "bold"),
    axis.title   = element_text(size = 22),
    axis.text    = element_text(size = 20),
    axis.text.x  = element_text(angle = 45, hjust = 1))





























library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(lubridate)

# ----------------------------
# Compute yearly mean sap flow per species
# ----------------------------
sapflow_yearly_species <- sapflow_filled %>%
  select(date, all_of(sensor_cols)) %>%
  pivot_longer(
    cols = all_of(sensor_cols),
    names_to = "sensor",
    values_to = "sapflow"
  ) %>%
  left_join(sensor_meta %>% select(sensor, species), by = "sensor") %>%
  mutate(year = year(date)) %>%
  group_by(year, species) %>%
  summarize(mean_sapflow = mean(sapflow, na.rm = TRUE), .groups = "drop")

# ----------------------------
# Define line types per species
# ----------------------------
species_linetypes <- rep_len(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), 
                             length.out = length(unique(sapflow_yearly_species$species)))
names(species_linetypes) <- unique(sapflow_yearly_species$species)

# ----------------------------
# Plot: yearly mean sap flow per species (black lines, different linetypes)
# ----------------------------
windows(width = 14, height = 10)

ggplot(sapflow_yearly_species, aes(x = year, y = mean_sapflow, linetype = species)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(color = "black", size = 3) +
  geom_text(aes(label = round(mean_sapflow, 1)), 
            vjust = -0.7, size = 3.5, check_overlap = TRUE) +
  scale_linetype_manual(values = species_linetypes) +
  scale_x_continuous(breaks = unique(sapflow_yearly_species$year)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(
    title = "Yearly Mean Sap Flow by Species",
    x = "Year",
    y = "Mean Sap Flow",
    linetype = "Species"
  ) +
  theme_bw(base_size = 28) +
  theme(
    plot.title   = element_text(size = 28, face = "bold"),
    axis.title   = element_text(size = 22),
    axis.text    = element_text(size = 20),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )





















































############################################################
# Species-LEVEL YEARLY SAPFLOW MEAN & SD
# - Month-based completeness filtering
# - Remove dead tree sensor PIPO4m
# - Only selected Species: PIPO, PISF, PSME
# - Sensors named like PIPO1s / PIPO1n
# - Species = first 4 letters
# - Aspect  = last letter (n/s)
# - Plot yearly mean ± SD
############################################################

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(ggplot2)

# ----------------------------------------------------------
# 1. LOAD DATA
# ----------------------------------------------------------
sapflow_filled <- read.csv(
  "\\\\snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/sapflow/sapflux_final_Species_xgb.csv",
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------
# 2. IDENTIFY SENSOR COLUMNS
# ----------------------------------------------------------
sensor_cols <- setdiff(names(sapflow_filled), c("date", "YEAR"))

# Remove dead tree sensor
sensor_cols <- setdiff(sensor_cols, "PIPO4m")

# ----------------------------------------------------------
# 3. CONVERT TO LONG FORMAT + PARSE METADATA
# ----------------------------------------------------------
sap_long <- sapflow_filled %>%
  mutate(
    date  = as.Date(date),
    YEAR  = year(date),
    MONTH = month(date)
  ) %>%
  pivot_longer(
    cols = all_of(sensor_cols),
    names_to = "sensor",
    values_to = "sapflow"
  ) %>%
  mutate(
    Species = str_sub(sensor, 1, 4),
    aspect  = str_sub(sensor, -1)
  )

# ----------------------------------------------------------
# 4. KEEP ONLY SELECTED Species
# ----------------------------------------------------------
selected_Species <- c("PIPO", "PISF", "PSME")
sap_long <- sap_long %>%
  filter(Species %in% selected_Species)

# ----------------------------------------------------------
# 5. IDENTIFY COMPLETE SENSOR–YEARS (MONTH-BASED)
#    Require at least 10 months with valid data
# ----------------------------------------------------------
sensor_year_months <- sap_long %>%
  filter(!is.na(sapflow)) %>%
  group_by(sensor, YEAR) %>%
  summarise(
    n_months = n_distinct(MONTH),
    .groups = "drop"
  )

complete_sensor_years <- sensor_year_months %>%
  filter(n_months == 12)

# ----------------------------------------------------------
# 6. PRINT OMITTED SENSORS/YEARS
# ----------------------------------------------------------
omitted <- anti_join(sensor_year_months, complete_sensor_years, by = c("sensor", "YEAR"))

if(nrow(omitted) > 0){
  print("Omitted sensors / years due to insufficient months (<10):")
  print(omitted)
} else {
  print("All sensor-years are complete.")
}

# ----------------------------------------------------------
# 7. FILTER TO COMPLETE DATA ONLY
# ----------------------------------------------------------
sap_complete <- sap_long %>%
  semi_join(complete_sensor_years, by = c("sensor", "YEAR"))

# ----------------------------------------------------------
# 8. Species-LEVEL YEARLY STATISTICS
#    (North + South pooled)
# ----------------------------------------------------------
Species_year_stats <- sap_complete %>%
  group_by(Species, YEAR) %>%
  summarise(
    sum_sapflow = sum(sapflow, na.rm = TRUE),
    sd_sapflow   = sd(sapflow, na.rm = TRUE),
    n_obs        = n(),
    n_sensors    = n_distinct(sensor),
    .groups = "drop"
  ) %>%
  arrange(Species, YEAR)


Species_year_stats_mean <- sap_complete %>%
  group_by(Species, YEAR) %>%
  summarise(
    mean_sapflow = mean(sapflow, na.rm = TRUE),
    sd_sapflow   = sd(sapflow, na.rm = TRUE),
    n_obs        = n(),
    n_sensors    = n_distinct(sensor),
    .groups = "drop"
  ) %>%
  arrange(Species, YEAR)

# ----------------------------------------------------------
# 9. SAVE OUTPUT
# ----------------------------------------------------------
write.csv(
  Species_year_stats,
  "Species_yearly_mean_sd_month_complete_selected.csv",
  row.names = FALSE
)















############################################################
# YEARLY MEAN SAPFLOW PLOT (CLEAN DATA, MONTH-COMPLETE)
# - Only selected Species: PIPO, PISF, PSME
# - Dead tree sensors removed (PIPO4m)
# - Only complete sensor-years (≥10 months)
# - North + South pooled
############################################################

library(dplyr)
library(ggplot2)

# ----------------------------------------------------------
# 1. Filter for selected Species
# ----------------------------------------------------------
selected_Species <- c("PIPO", "PISF", "PSME")

# Species_year_stats already contains clean, month-complete data
sapflow_yearly_clean <- Species_year_stats %>%
  filter(Species %in% selected_Species)
sapflow_yearly_clean_mean <- Species_year_stats_mean %>%
  filter(Species %in% selected_Species)

# ----------------------------------------------------------
# 2. Define line types per Species
# ----------------------------------------------------------
Species_linetypes <- c("solid", "dashed", "dotted")
names(Species_linetypes) <- selected_Species

# ----------------------------------------------------------
# 3. Plot yearly mean ± SD
# ----------------------------------------------------------
windows()
ggplot(sapflow_yearly_clean, aes(x = YEAR, y = sum_sapflow, linetype = Species)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(color = "black", size = 3) +  # SD as shaded band
  geom_text(aes(label = round(sum_sapflow, 1)),
            vjust = -0.7,
            size = 4,
            check_overlap = TRUE) +
  scale_linetype_manual(values = Species_linetypes) +
  scale_fill_manual(values = c("PIPO" = "grey70", "PISF" = "grey50", "PSME" = "grey30")) +
  scale_x_continuous(breaks = unique(sapflow_yearly_clean$YEAR)) +
  labs(
    title = "Yearly Sum Sapflow by Species (Clean Data)",
    x = "Year",
    y = "Mean Sapflow",
    linetype = "Species",
    fill = "Species"
  ) +
  theme_bw(base_size = 20) +
  theme(
    plot.title   = element_text(size = 24, face = "bold"),
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 16),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

windows()
ggplot(sapflow_yearly_clean_mean, aes(x = YEAR, y = mean_sapflow, linetype = Species)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(color = "black", size = 3) +  # SD as shaded band
  geom_text(aes(label = round(mean_sapflow, 1)),
            vjust = -0.7,
            size = 4,
            check_overlap = TRUE) +
  scale_linetype_manual(values = Species_linetypes) +
  scale_fill_manual(values = c("PIPO" = "grey70", "PISF" = "grey50", "PSME" = "grey30")) +
  scale_x_continuous(breaks = unique(sapflow_yearly_clean$YEAR)) +
  labs(
    title = "Yearly Mean Sapflow by Species (Clean Data)",
    x = "Year",
    y = "Mean Sapflow",
    linetype = "Species",
    fill = "Species"
  ) +
  theme_bw(base_size = 20) +
  theme(
    plot.title   = element_text(size = 24, face = "bold"),
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 16),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )






overall_Species_stats <- sapflow_yearly_clean_mean %>%
  group_by(Species) %>%
  summarise(
    overall_mean = mean(mean_sapflow, na.rm = TRUE),
    overall_sd   = sd(mean_sapflow, na.rm = TRUE),
    n_years      = n(),
    .groups = "drop"
  )

print(overall_Species_stats)





############################################################
# SEASONAL MEAN SAPFLOW PLOT (CLEAN DATA, MONTH-COMPLETE)
# - Only selected Species: PIPO, PISF, PSME
# - Dead tree sensors removed (PIPO4m)
# - Only complete sensor-years (≥10 months)
# - North + South pooled
############################################################

library(dplyr)
library(ggplot2)
library(lubridate)
library(ggpattern)

# ----------------------------------------------------------
# 1. Define seasons
# ----------------------------------------------------------
sap_complete_season <- sap_complete %>%
  filter(Species %in% c("PIPO", "PISF", "PSME")) %>%
  mutate(
    season = case_when(
      MONTH %in% 7:9  ~ "monsoon",
      MONTH %in% 4:6  ~ "pre-monsoon",
      MONTH == 10     ~ "fall",
      MONTH %in% c(11,12,1,2,3) ~ "winter",
      TRUE ~ NA_character_
    )
  )

# ----------------------------------------------------------
# 2. Compute seasonal mean ± SD
# ----------------------------------------------------------
seasonal_stats <- sap_complete_season %>%
  group_by(Species, season) %>%
  summarise(
    sum_sapflow = sum(sapflow, na.rm = TRUE),
    sd_sapflow   = sd(sapflow, na.rm = TRUE),
    n_obs        = n(),
    .groups = "drop"
  )

seasonal_stats_mean <- sap_complete_season %>%
  group_by(Species, season) %>%
  summarise(
    mean_sapflow = mean(sapflow, na.rm = TRUE),
    sd_sapflow   = sd(sapflow, na.rm = TRUE),
    n_obs        = n(),
    .groups = "drop"
  )

# ----------------------------------------------------------
# 3. Define patterns per Species
# ----------------------------------------------------------
Species_list <- unique(seasonal_stats$Species)
pattern_types <- c("stripe", "crosshatch", "circle")  # one per Species
Species_patterns <- setNames(pattern_types[1:length(Species_list)], Species_list)

# ----------------------------------------------------------
# 4. Plot seasonal mean ± SD
# ----------------------------------------------------------
windows()
ggplot(seasonal_stats, aes(x = season, y = sum_sapflow, pattern = Species)) +
  geom_col_pattern(
    fill = "black",
    color = "black",
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 0.9,
    pattern_density = 0.5,
    pattern_spacing = 0.05
  ) +
  geom_text(aes(label = round(sum_sapflow, 2)),
            position = position_dodge(width = 0.8),
            vjust = -0.5,
            size = 6,
            color = "black") +
  scale_pattern_manual(values = Species_patterns) +
  labs(
    title = "Seasonal Sum Sapflow by Species (Clean Data)",
    x = "Season",
    y = "Mean Sap Flow",
    pattern = "Species"
  ) +
  theme_bw(base_size = 20) +
  theme(
    plot.title   = element_text(size = 24, face = "bold"),
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 16),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )









windows()
ggplot(seasonal_stats_mean, aes(x = season, y = mean_sapflow, pattern = Species)) +
  geom_col_pattern(
    fill = "black",
    color = "black",
    position = position_dodge(width = 0.8),
    width = 0.6,
    alpha = 0.9,
    pattern_density = 0.5,
    pattern_spacing = 0.05
  ) +
  geom_text(aes(label = round(mean_sapflow, 2)),
            position = position_dodge(width = 0.8),
            vjust = -0.5,
            size = 6,
            color = "black") +
  scale_pattern_manual(values = Species_patterns) +
  labs(
    title = "Seasonal Mean Sapflow by Species (Clean Data)",
    x = "Season",
    y = "Mean Sap Flow",
    pattern = "Species"
  ) +
  theme_bw(base_size = 20) +
  theme(
    plot.title   = element_text(size = 24, face = "bold"),
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 16),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )















































################SAPFLUX-relationship with RWI################


###DBH
#PIPO: 22.3, 29.2,43.5,27.1
#PISF: 31.9, 33.4, 49.8
#PSME: 20, 37.1


# Define DBH values per Species
dbh_list <- list(
  PIPO = c(22.3, 29.2, 43.5, 27.1),
  PISF = c(31.9, 33.4, 49.8),
  PSME = c(20, 37.1)
)

# Tolerance (e.g., ±5 cm)
tol <- 1

# Subset trees with DBH within tolerance
filter_by_Species_dbh <- function(df, dbh_list, tol) {
  df %>%
    filter(Species %in% names(dbh_list)) %>%
    rowwise() %>%
    filter(any(abs(DBH - dbh_list[[Species]]) <= tol)) %>%
    ungroup()
}


MBA_sub <- filter_by_Species_dbh(MBA_treering_ground, dbh_list, tol)
MBB_sub <- filter_by_Species_dbh(MBB_treering_ground, dbh_list, tol)
MBC_sub <- filter_by_Species_dbh(MBC_treering_ground, dbh_list, tol)


sapflow_yearly_clean
sapflow_yearly_clean_mean


library(dplyr)
library(ggplot2)

# ==========================================================
# 1. SAPFLOW: Species-level year availability (background)
# ==========================================================

sapflow_Species_years <- sapflow_yearly_clean_mean %>%
  distinct(Species, YEAR) %>%
  mutate(
    y = 0,                 # fixed background position
    source = "Sapflow"
  )

library(dplyr)
library(tidyr)

transpose_treering <- function(df, site_name) {
  
  # Identify year columns automatically (numeric column names)
  year_cols <- colnames(df)[grepl("^[0-9]{4}$", colnames(df))]
  
  # Pivot years to long format
  df_long <- df %>%
    pivot_longer(
      cols = all_of(year_cols),
      names_to = "YEAR",
      values_to = "value"
    ) %>%
    mutate(
      YEAR = as.integer(YEAR),
      site = site_name
    ) %>%
    select(site, treeID, YEAR, value, Species, DBH) %>%
    filter(!is.na(value))
  
  return(df_long)
}

# =============================
# Apply to your three datasets
# =============================
MBA_long <- transpose_treering(MBA_sub, "MBA")
MBB_long <- transpose_treering(MBB_sub, "MBB")
MBC_long <- transpose_treering(MBC_sub, "MBC")

# Combine all three
treering_tree_years <- bind_rows(MBA_long, MBB_long, MBC_long)



treering_long_std <- treering_tree_years %>%
  group_by(Species, treeID) %>%
  mutate(
    ring_z = scale(value)[,1]
  ) %>%
  ungroup()

sapflow_Species_mean <- sapflow_yearly_clean_mean %>%
  select(Species = Species, YEAR, mean_sapflow) %>%
  filter(!is.na(mean_sapflow))



common_years <- intersect(
  unique(sapflow_Species_mean$YEAR),
  unique(treering_long_std$YEAR)
)

sapflow_Species_mean <- sapflow_Species_mean %>%
  filter(YEAR %in% common_years)

treering_long_std <- treering_long_std %>%
  filter(YEAR %in% common_years)





windows()

ggplot() +
  
  # ---- Individual tree-ring series (thin) ----
geom_line(
  data = treering_long_std,
  aes(x = YEAR, y = ring_z, group = treeID),
  color = "#1b9e77",
  alpha = 0.5,
  linewidth = 0.5
) +
  
  # ---- Species mean sapflow (thick) ----
geom_line(
  data = sapflow_Species_mean,
  aes(x = YEAR, y = mean_sapflow),
  color = "black",
  linewidth = 1.3
) +
  
  facet_wrap(~ Species, scales = "free_y") +
  
  labs(
    x = "Year",
    y = "Sapflow mean (black) / Tree-ring z-score (green)",
    title = "Species-Level Sapflow and Individual Tree-Ring Series"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )







library(dplyr)

cor_stats <- treering_long_std %>%
  left_join(
    sapflow_Species_mean,
    by = c("YEAR", "Species")
  ) %>%
  group_by(Species) %>%
  summarize(
    r = cor(ring_z, mean_sapflow, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("r = ", round(r, 2))
  )
windows()

ggplot() +
  # ---- Individual tree-ring series ----
geom_line(
  data = treering_long_std,
  aes(x = YEAR, y = ring_z, group = treeID),
  color = "#1b9e77",
  alpha = 0.5,
  linewidth = 0.5
) +
  
  # ---- Species mean sapflow ----
geom_line(
  data = sapflow_Species_mean,
  aes(x = YEAR, y = mean_sapflow),
  color = "black",
  linewidth = 1.3
) +
  
  # ---- Correlation annotation ----
geom_text(
  data = cor_stats,
  aes(x = -Inf, y = Inf, label = label),
  hjust = -0.1,
  vjust = 1.2,
  inherit.aes = FALSE,
  size = 3.5,
  fontface = "bold"
) +
  
  facet_wrap(~ Species, scales = "free_y") +
  
  labs(
    x = "Year",
    y = "Sapflow mean (black) / Tree-ring z-score (green)",
    title = "Species-Level Sapflow and Individual Tree-Ring Series"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )











seasonal_sapflow <- sap_complete_season %>%
  filter(Species %in% c("PIPO", "PISF", "PSME")) %>%
  group_by(Species, YEAR, season) %>%
  rename() %>%
  summarise(
    seasonal_mean_sapflow = mean(sapflow, na.rm = TRUE),
    .groups = "drop"
  )

yearly_seasonal_sap <- sap_complete_season %>%
  filter(Species %in% c("PIPO", "PISF", "PSME")) %>%
  group_by(Species, YEAR, season) %>%
  summarise(
    seasonal_mean_sapflow = mean(sapflow, na.rm = TRUE),
    .groups = "drop"
  )
sap_ring_combined <- yearly_seasonal_sap %>%
  left_join(treering_tree_years, by = c("Species", "YEAR"))
library(dplyr)

seasonal_correlations <- sap_ring_combined %>%
  group_by(Species, season) %>%
  summarise(
    correlation = cor(seasonal_mean_sapflow, value, use = "complete.obs"),
    .groups = "drop"
  )

seasonal_correlations

library(ggplot2)

ggplot(seasonal_correlations, aes(x = season, y = correlation, fill = season)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(correlation, 2)), 
            vjust = ifelse(seasonal_correlations$correlation >= 0, -0.3, 1.2),
            size = 4) +
  facet_wrap(~Species) +
  scale_fill_brewer(palette = "Set2") +
  ylim(-1, 1) +
  labs(
    title = "Seasonal Sap Flux vs Tree-Ring Correlation by Species",
    x = "Season",
    y = "Pearson Correlation"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none"
  )







################SAPFLUX-relationship with RWI################




