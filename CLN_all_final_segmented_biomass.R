# =============================================================================
# Point Density Calculator by Scan Type
# =============================================================================
# Calculates mean point density (pts/m²) and standard deviation per scan type
# (nadir, nadir_oblique, TLS, TLS_UAV) across plots (MBA, MBB, MBC).
#
# Requirements:
#   install.packages(c("lidR", "dplyr", "tibble"))
#
# Usage: Open in RStudio and click Source, or run line by line.
# =============================================================================
if (FALSE) {
library(lidR)
library(dplyr)
library(tibble)

# ── CONFIG ────────────────────────────────────────────────────────────────────
base_dir <- "Z:/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/CH2_LiDAR/segmentedClouds"

# Each scan type maps to its three plot folders
scan_types <- list(
  nadir         = c("nadir_MBA",         "nadir_MBB",         "nadir_MBC"),
  nadir_oblique = c("nadir_oblique_MBA",  "nadir_oblique_MBB", "nadir_oblique_MBC"),
  TLS           = c("TLS_MBA",            "TLS_MBB",           "TLS_MBC"),
  TLS_UAV       = c("TLS_UAV_MBA",        "TLS_UAV_MBB",       "TLS_UAV_MBC")
)
# ─────────────────────────────────────────────────────────────────────────────


# ── HELPER: get density for one folder ───────────────────────────────────────
get_folder_density <- function(folder_path) {
  
  if (!dir.exists(folder_path)) {
    warning(paste("Folder not found:", folder_path))
    return(NA_real_)
  }
  
  # Find all .las / .laz files (recursive)
  files <- list.files(folder_path, pattern = "\\.la[sz]$",
                      full.names = TRUE, recursive = TRUE)
  
  if (length(files) == 0) {
    warning(paste("No .las/.laz files found in:", folder_path))
    return(NA_real_)
  }
  
  # Read all files into one LAScatalog or LAS object
  if (length(files) == 1) {
    las <- readLAS(files, select = "xyz")
  } else {
    ctg <- readLAScatalog(folder_path)
    las <- readLAS(ctg, select = "xyz")
  }
  
  if (is.empty(las)) {
    warning(paste("Empty point cloud in:", folder_path))
    return(NA_real_)
  }
  
  # Compute density: total points / convex hull area (m²)
  n_pts  <- nrow(las@data)
  coords <- cbind(las@data$X, las@data$Y)
  
  # Convex hull area via cross-product method
  hull_idx  <- chull(coords)
  hull_pts  <- coords[c(hull_idx, hull_idx[1]), ]  # close the polygon
  # Shoelace formula
  x <- hull_pts[, 1]
  y <- hull_pts[, 2]
  area_m2 <- 0.5 * abs(sum(x[-length(x)] * y[-1]) - sum(x[-1] * y[-length(y)]))
  
  if (area_m2 == 0) return(NA_real_)
  
  density <- n_pts / area_m2
  return(density)
}


# ── MAIN LOOP ─────────────────────────────────────────────────────────────────
cat("=== Point Density Analysis ===\n\n")

results_list <- list()

for (scan_name in names(scan_types)) {
  
  folders <- scan_types[[scan_name]]
  densities <- numeric(length(folders))
  plot_names <- character(length(folders))
  
  cat(sprintf("Processing scan type: %s\n", scan_name))
  
  for (i in seq_along(folders)) {
    folder_path <- file.path(base_dir, folders[i])
    cat(sprintf("  [%d/%d] %s ... ", i, length(folders), folders[i]))
    
    d <- tryCatch(
      get_folder_density(folder_path),
      error = function(e) {
        cat(sprintf("ERROR: %s\n", e$message))
        NA_real_
      }
    )
    
    if (!is.na(d)) cat(sprintf("%.2f pts/m²\n", d))
    densities[i] <- d
    plot_names[i] <- folders[i]
  }
  
  valid <- densities[!is.na(densities)]
  
  results_list[[scan_name]] <- tibble(
    scan_type   = scan_name,
    plot        = plot_names,
    density_m2  = densities
  )
  
  cat(sprintf(
    "  --> Mean: %.2f pts/m²  |  SD: %.2f pts/m²  |  n = %d\n\n",
    ifelse(length(valid) > 0, mean(valid), NA),
    ifelse(length(valid) > 1, sd(valid),   NA),
    length(valid)
  ))
}


# ── SUMMARY TABLE ─────────────────────────────────────────────────────────────
all_results <- bind_rows(results_list)

summary_table <- all_results %>%
  group_by(scan_type) %>%
  summarise(
    n_plots      = sum(!is.na(density_m2)),
    mean_density = mean(density_m2, na.rm = TRUE),
    sd_density   = sd(density_m2,   na.rm = TRUE),
    min_density  = min(density_m2,  na.rm = TRUE),
    max_density  = max(density_m2,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_density))

cat("=== SUMMARY: Mean Point Density by Scan Type ===\n")
print(summary_table, n = Inf)

cat("\n=== PER-PLOT DENSITIES ===\n")
print(all_results, n = Inf)


# ── OPTIONAL: Save results to CSV ─────────────────────────────────────────────
out_dir <- dirname(base_dir)  # saves one level up from segmentedClouds

#write.csv(summary_table, file.path(out_dir, "density_summary.csv"),     row.names = FALSE)
#write.csv(all_results,   file.path(out_dir, "density_per_plot.csv"),     row.names = FALSE)

cat(sprintf("\nResults saved to:\n  %s\n  %s\n",
            file.path(out_dir, "density_summary.csv"),
            file.path(out_dir, "density_per_plot.csv")))




}







library(lidR)
library(VoxR)
library(data.table)
library(geometry)
library(stringr)
library(tidyverse)
library(sf)
library(randomForest)
#remotes::install_github('bi0m3trics/spanner')

Packages <- c('lidR','data.table','geometry','VoxR','stringr','tidyverse','readxl',"dplR","dplyr","sf", "lidR", "sf","tibble","ggplot2","tidyverse",'spanner')
lapply(Packages, library, character.only = TRUE)

setwd("//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/")
gc()
options(scipen = 100, digits = 4)
############################################ Tree Ring and Ground Fusion ####################################################
#SET the path for folder where you have the file stored, each function will load the data format
groundA =   readxl::read_excel('//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Metadata/MtBigelow_Metadata/metadata_MtBigelow_Spring_2022.xlsx', sheet=1)
groundB =   readxl::read_excel('//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Metadata/MtBigelow_Metadata/metadata_MtBigelow_Spring_2022.xlsx', sheet=2)
groundC =   readxl::read_excel('//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/NorthAmerica Data/MtBigelow_Metadata/MtBigelow_Metadata/metadata_MtBigelow_Spring_2022.xlsx', sheet=3) 


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
      add_column(species = dist_azim$Species[i]) %>% #add Species column of that tree
      add_column(Height_m = dist_azim$Height_m[i])%>% #add Height_m column of that tree
      add_column(DBH = dist_azim$DBH[i]) #add DBH column of that tree
    
    
    coor_list[[i]] <- new_coordinates #put result of loop in list
    
  }
  
  coor_df<- coor_list %>% bind_rows() %>% na.omit() #remove NA rows
  
  
  return(coor_df)
}

MBA_ground <- ground_parse(MBA_coor, groundA) 
MBB_ground <- ground_parse(MBB_coor, groundB) 
MBC_ground <- ground_parse(MBC_coor, groundC) 



##Creating ground SPATIALPOINTS
###tree_ring_ground_spatialpoints
#Conversion of data frame to sf object
MBA_xy_ground <- st_as_sf(x = MBA_ground,                         
                          coords = c("lon", "lat"),
                          crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")
#Projection transformation
MBA_tree.points.ground = st_transform(MBA_xy_ground, crs = "+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs")
#Convert it to data frame
MBA_tree.points.ground_df <- MBA_tree.points.ground %>% as.data.frame() %>%
  mutate(X = unlist(map(geometry,1)),
         Y = unlist(map(geometry,2)))
###write_sf(MBA_tree.points.ground, "MBA_tree_points_ground.shp")
###write.csv(MBA_tree.points.ground_df, "MBA_tree_points_ground.csv")

##SPATIALPOINTS
###tree_ring_ground_spatialpoints
#Conversion of data frame to sf object
MBB_xy_ground <- st_as_sf(x = MBB_ground,                         
                          coords = c("lon", "lat"),
                          crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")
#Projection transformation
MBB_tree.points.ground = st_transform(MBB_xy_ground, crs = "+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs")
#Convert it to data frame
MBB_tree.points.ground_df <- MBB_tree.points.ground %>% as.data.frame() %>%
  mutate(X = unlist(map(geometry,1)),
         Y = unlist(map(geometry,2)))
###write_sf(MBB_tree.points.ground, "MBB_tree_points_ground.shp")
###write.csv(MBB_tree.points.ground_df, "MBB_tree_points_ground.csv")


##SPATIALPOINTS
###tree_ring_ground_spatialpoints
#Conversion of data frame to sf object
MBC_xy_ground <- st_as_sf(x = MBC_ground,                         
                          coords = c("lon", "lat"),
                          crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")
#Projection transformation
MBC_tree.points.ground = st_transform(MBC_xy_ground, crs = "+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs")
#Convert it to data frame
MBC_tree.points.ground_df <- MBC_tree.points.ground %>% as.data.frame() %>%
  mutate(X = unlist(map(geometry,1)),
         Y = unlist(map(geometry,2)))

###write_sf(MBC_tree.points.ground, "MBC_tree_points_ground.shp")
###write.csv(MBB_tree.points.ground_df, "MBB_tree_points_ground.csv")


#theisland <- list.files('/gaea/Downloads/',pattern="*.laz$", full.names=TRUE)


all_trees_bigelow <- MBA_tree.points.ground_df%>%
  mutate(site = "MBA") %>%
  rbind(MBB_tree.points.ground_df %>% mutate(site = "MBB")) %>% 
  rbind(MBC_tree.points.ground_df %>% mutate(site = "MBC")) %>% 
  rownames_to_column()%>%rename(ID=rowname) 
























density_measurements <- read_xlsx("//snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/Density Measurments Mt Big.xlsx")%>% select(c(1,4))
density_measurements_summary <- density_measurements %>%
  mutate(
    site = substr(`Sample #`, 1, 1),                                   # first letter
    groundID = as.numeric(gsub("[^0-9]", "", `Sample #`)),            # numbers
    core = sub(".*([A-Za-z])$", "\\1", `Sample #`)                     # last letter
  ) %>%
  group_by(site, groundID) %>%
  summarise(
    mean_density = mean(Density, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    site = dplyr::recode(site,
                         "A" = "MBA",
                         "B" = "MBB",
                         "C" = "MBC")
  )
density_measurements_summary_species <-density_measurements_summary%>% left_join(all_trees_bigelow)%>%
  group_by(species) %>%
  summarise(
    mean_density = mean(mean_density, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

















































if (FALSE) {










library(lidR)
library(data.table)
library(FNN)
library(randomForest)

# ── SETTINGS ──────────────────────────────────────────────────────────────────
parent_dir           <- "//snre-snow.snrenet.arizona.edu/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/CH2_LiDAR/segmentedClouds/"
training_leaf        <- "//snre-snow.snrenet.arizona.edu/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/CH2_LiDAR/training/merged_training/leaf_training_V2.las"
training_woody       <- "//snre-snow.snrenet.arizona.edu/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/CH2_LiDAR/training/merged_training/wood_training_V2.las"

k_nn                 <- 10
k_folds              <- 5
woody_prob_threshold <- 0.35

# ── Crown Base Height_m detection ───────────────────────────────────────────────
# n_slices  : number of horizontal layers to cut the cloud into
# gap_pct   : if a layer has fewer points than this fraction of the max-layer
#             count, it is considered a "gap" layer.  The highest contiguous
#             gap layer from the bottom defines the CBH.
# low_band_m: metres above CBH that are still considered "low canopy"
#             (used for the sparse/non-green override)
cbh_n_slices         <- 50
cbh_gap_pct          <- 0.10
low_canopy_band_m    <- 1.0

# ── Low-canopy no-leaf thresholds ─────────────────────────────────────────────
# A point in the low-canopy band is forced woody when ALL three hold:
#   local density  < low_density_pct  * median density of the whole cloud
#   return ratio   < low_return_thr   (single-return dominant)
#   greenness ExG  < low_exg_thr      (normalised 0-1 scale)
low_density_pct      <- 0.40
low_return_thr       <- 0.55
low_exg_thr          <- 0.30

set.seed(42)

# ── HELPERS ───────────────────────────────────────────────────────────────────
compute_density <- function(pts, k = k_nn) {
  coords <- cbind(pts$X, pts$Y, pts$Z)
  nn     <- FNN::get.knn(coords, k = k)
  1 / (rowMeans(nn$nn.dist) + 1e-6)
}

safe_norm <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

compute_colour_indices <- function(R, G, B) {
  ExG       <-  2*G - R - B
  VARI      <- (G - R)       / (G + R - B + 1e-6)
  GLI       <- (2*G - R - B) / (2*G + R + B + 1e-6)
  GR        <-  G / (R + 1e-6)
  BI        <- (R + G - B)   / (R + G + B + 1e-6)
  NGRDI_inv <- (R - G)       / (R + G + 1e-6)
  RG_excess <-  2*R - G - B
  GBI       <-  ExG - BI
  list(ExG=ExG, VARI=VARI, GLI=GLI, GR=GR,
       BI=BI, NGRDI_inv=NGRDI_inv, RG_excess=RG_excess, GBI=GBI)
}

extract_features <- function(pts) {
  has_rgb <- all(c("R","G","B") %in% names(pts))
  I_raw   <- if ("Intensity" %in% names(pts)) pts$Intensity else rep(0, nrow(pts))
  I_n     <- safe_norm(I_raw)
  Rn      <- if (all(c("ReturnNumber","NumberOfReturns") %in% names(pts)))
    pts$ReturnNumber / pmax(pts$NumberOfReturns, 1) else rep(0, nrow(pts))
  D_n     <- safe_norm(compute_density(pts))
  H_rel   <- safe_norm(pts$Z)
  
  if (has_rgb) {
    R <- safe_norm(pts$R); G <- safe_norm(pts$G); B <- safe_norm(pts$B)
    ci <- compute_colour_indices(R, G, B)
  } else {
    R <- G <- B <- rep(0, nrow(pts))
    ci <- list(ExG=R, VARI=R, GLI=R, GR=R,
               BI=R, NGRDI_inv=R, RG_excess=R, GBI=R)
  }
  
  data.table(
    GBI=ci$GBI, ExG=ci$ExG, BI=ci$BI,
    VARI=ci$VARI, GLI=ci$GLI, GR=ci$GR,
    NGRDI_inv=ci$NGRDI_inv, RG_excess=ci$RG_excess,
    D_n=D_n, H_rel=H_rel, I_n=I_n, Rn=Rn,
    R=R, G=G, B=B
  )
}

feature_cols <- c(
  "GBI","GBI","GBI","GBI","GBI",
  "ExG","ExG","ExG",
  "BI","BI","BI",
  "VARI","GLI","GR",
  "NGRDI_inv","RG_excess",
  "D_n","D_n",
  "H_rel","I_n","Rn",
  "R","G","B"
)

# ── CROWN BASE Height_m DETECTION ───────────────────────────────────────────────
# Returns the Z value (absolute) of the detected crown base.
# Strategy: slice the cloud into n_slices equal Height_m bands, count points
# per band, then walk upward from the bottom and find the first band where
# density exceeds gap_pct * max_band_count — that transition is the CBH.
detect_cbh <- function(z, n_slices = cbh_n_slices, gap_pct = cbh_gap_pct) {
  if (length(z) < 50) return(min(z))         # too few points — treat all as below
  
  breaks  <- seq(min(z), max(z), length.out = n_slices + 1)
  mids    <- (breaks[-1] + breaks[-(n_slices+1)]) / 2
  counts  <- tabulate(findInterval(z, breaks, rightmost.closed = TRUE),
                      nbins = n_slices)
  
  max_cnt <- max(counts)
  is_gap  <- counts < gap_pct * max_cnt
  
  # Walk upward: find the highest contiguous gap layer starting from the bottom
  cbh_idx <- 1L
  for (i in seq_along(is_gap)) {
    if (is_gap[i]) cbh_idx <- i else break
  }
  
  mids[cbh_idx]
}

# ── ACCURACY ASSESSMENT ───────────────────────────────────────────────────────
compute_accuracy <- function(observed, predicted) {
  obs  <- factor(observed,  levels = c("leaf","woody"))
  pred <- factor(predicted, levels = c("leaf","woody"))
  cm   <- table(Observed = obs, Predicted = pred)
  
  oa    <- sum(diag(cm)) / sum(cm)
  n     <- sum(cm)
  p_exp <- sum(rowSums(cm) * colSums(cm)) / n^2
  kappa <- (oa - p_exp) / (1 - p_exp + 1e-9)
  
  metrics <- rbindlist(lapply(rownames(cm), function(cls) {
    tp        <- cm[cls, cls]
    fp        <- sum(cm[, cls]) - tp
    fn        <- sum(cm[cls, ]) - tp
    precision <- tp / (tp + fp + 1e-9)
    recall    <- tp / (tp + fn + 1e-9)
    f1        <- 2 * precision * recall / (precision + recall + 1e-9)
    data.table(Class        = cls,
               TP           = tp,  FP = fp,  FN = fn,
               Producer_Acc = round(recall,    4),
               User_Acc     = round(precision, 4),
               F1           = round(f1,        4))
  }))
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════╗\n")
  cat("║           ACCURACY ASSESSMENT RESULTS           ║\n")
  cat("╠══════════════════════════════════════════════════╣\n")
  cat("║  Confusion Matrix:                               ║\n")
  cat("║                 Predicted                        ║\n")
  cat("║  Observed    leaf     woody                      ║\n")
  cat(sprintf("║    leaf    %6d   %6d                          ║\n",
              cm["leaf","leaf"], cm["leaf","woody"]))
  cat(sprintf("║    woody   %6d   %6d                          ║\n",
              cm["woody","leaf"], cm["woody","woody"]))
  cat("╠══════════════════════════════════════════════════╣\n")
  cat(sprintf("║  Overall Accuracy  : %6.2f %%                  ║\n", oa * 100))
  cat(sprintf("║  Cohen's Kappa     : %7.4f                    ║\n", kappa))
  cat("╠══════════════════════════════════════════════════╣\n")
  cat(sprintf("║  LEAF  — Producer : %6.2f %%  User : %6.2f %%  ║\n",
              metrics[Class=="leaf",  Producer_Acc]*100,
              metrics[Class=="leaf",  User_Acc]*100))
  cat(sprintf("║  LEAF  — F1       : %7.4f                    ║\n",
              metrics[Class=="leaf",  F1]))
  cat("║                                                  ║\n")
  cat(sprintf("║  WOODY — Producer : %6.2f %%  User : %6.2f %%  ║\n",
              metrics[Class=="woody", Producer_Acc]*100,
              metrics[Class=="woody", User_Acc]*100))
  cat(sprintf("║  WOODY — F1       : %7.4f                    ║\n",
              metrics[Class=="woody", F1]))
  cat("╚══════════════════════════════════════════════════╝\n\n")
  
  invisible(list(confusion_matrix=cm, metrics=metrics,
                 overall_accuracy=oa, kappa=kappa))
}

# ── LOAD TRAINING DATA ────────────────────────────────────────────────────────
cat("Loading training files...\n")

las_leaf_train  <- readLAS(training_leaf)
las_woody_train <- readLAS(training_woody)
if (is.empty(las_leaf_train))  stop("Training leaf file is empty.")
if (is.empty(las_woody_train)) stop("Training woody file is empty.")

feat_leaf  <- extract_features(las_leaf_train@data);  feat_leaf[,  label := "leaf"]
feat_woody <- extract_features(las_woody_train@data); feat_woody[, label := "woody"]

n_leaf_raw  <- nrow(feat_leaf)
n_woody_raw <- nrow(feat_woody)

cat("  Raw — leaf:", n_leaf_raw, "| woody:", n_woody_raw, "\n")
cat("  Imbalance ratio:", round(n_leaf_raw / max(n_woody_raw, 1), 2), ":1\n")

# ── BALANCE TRAINING DATA ─────────────────────────────────────────────────────
target_n <- n_leaf_raw

if (n_woody_raw < target_n) {
  idx_over   <- sample(n_woody_raw, target_n, replace = TRUE)
  feat_woody <- feat_woody[idx_over]
  cat("  Oversampled woody:", n_woody_raw, "->", nrow(feat_woody), "\n")
} else {
  feat_woody <- feat_woody[sample(n_woody_raw, target_n, replace = FALSE)]
  cat("  Downsampled woody to:", nrow(feat_woody), "\n")
}

all_labelled <- rbindlist(list(feat_leaf, feat_woody))
all_labelled[is.na(all_labelled)] <- 0
all_labelled[, label := as.factor(label)]

# ── K-FOLD CROSS-VALIDATION ───────────────────────────────────────────────────
cat("\nRunning", k_folds, "-fold cross-validation (stratified)...\n")

idx_leaf  <- which(all_labelled$label == "leaf")
idx_woody <- which(all_labelled$label == "woody")

folds_leaf  <- split(sample(idx_leaf),  rep(1:k_folds, length.out = length(idx_leaf)))
folds_woody <- split(sample(idx_woody), rep(1:k_folds, length.out = length(idx_woody)))

fold_metrics <- vector("list", k_folds)

for (k in seq_len(k_folds)) {
  cat(sprintf("  Fold %d/%d...\n", k, k_folds))
  
  val_idx    <- c(folds_leaf[[k]], folds_woody[[k]])
  train_idx  <- setdiff(seq_len(nrow(all_labelled)), val_idx)
  
  fold_train <- all_labelled[train_idx]
  fold_val   <- all_labelled[val_idx]
  
  fold_model <- randomForest(
    x      = fold_train[, ..feature_cols],
    y      = fold_train$label,
    ntree  = 500,
    importance = FALSE
  )
  
  val_probs <- predict(fold_model, fold_val[, ..feature_cols], type = "prob")
  val_preds <- ifelse(val_probs[, "woody"] >= woody_prob_threshold, "woody", "leaf")
  val_obs   <- as.character(fold_val$label)
  
  obs_f  <- factor(val_obs,   levels = c("leaf","woody"))
  pred_f <- factor(val_preds, levels = c("leaf","woody"))
  cm_f   <- table(obs_f, pred_f)
  
  oa_f    <- sum(diag(cm_f)) / sum(cm_f)
  n_f     <- sum(cm_f)
  p_exp_f <- sum(rowSums(cm_f) * colSums(cm_f)) / n_f^2
  kappa_f <- (oa_f - p_exp_f) / (1 - p_exp_f + 1e-9)
  
  get_f1 <- function(cls) {
    tp <- cm_f[cls, cls]; fp <- sum(cm_f[, cls]) - tp; fn <- sum(cm_f[cls, ]) - tp
    pr <- tp / (tp + fp + 1e-9); re <- tp / (tp + fn + 1e-9)
    2 * pr * re / (pr + re + 1e-9)
  }
  
  fold_metrics[[k]] <- data.table(
    fold        = k,
    overall_acc = oa_f,
    kappa       = kappa_f,
    leaf_f1     = get_f1("leaf"),
    woody_f1    = get_f1("woody")
  )
  
  cat(sprintf("    OA: %.2f%%  Kappa: %.4f  Leaf-F1: %.4f  Woody-F1: %.4f\n",
              oa_f * 100, kappa_f,
              fold_metrics[[k]]$leaf_f1, fold_metrics[[k]]$woody_f1))
}

cv_dt <- rbindlist(fold_metrics)

cat("\n── K-Fold CV Summary ─────────────────────────────────\n")
cat(sprintf("  Mean OA    : %.2f %%  (SD: %.2f %%)\n",
            mean(cv_dt$overall_acc)*100, sd(cv_dt$overall_acc)*100))
cat(sprintf("  Mean Kappa : %.4f   (SD: %.4f)\n",
            mean(cv_dt$kappa), sd(cv_dt$kappa)))
cat(sprintf("  Mean Leaf  F1: %.4f (SD: %.4f)\n",
            mean(cv_dt$leaf_f1),  sd(cv_dt$leaf_f1)))
cat(sprintf("  Mean Woody F1: %.4f (SD: %.4f)\n",
            mean(cv_dt$woody_f1), sd(cv_dt$woody_f1)))

# Aggregated CV confusion matrix (re-predict each fold)
all_obs_cv <- character(0); all_pred_cv <- character(0)

for (k in seq_len(k_folds)) {
  val_idx   <- c(folds_leaf[[k]], folds_woody[[k]])
  fold_val  <- all_labelled[val_idx]
  train_idx <- setdiff(seq_len(nrow(all_labelled)), val_idx)
  fold_train <- all_labelled[train_idx]
  
  fold_model <- randomForest(
    x = fold_train[, ..feature_cols], y = fold_train$label,
    ntree = 500, importance = FALSE
  )
  
  val_probs <- predict(fold_model, fold_val[, ..feature_cols], type = "prob")
  val_preds <- ifelse(val_probs[, "woody"] >= woody_prob_threshold, "woody", "leaf")
  
  all_obs_cv  <- c(all_obs_cv,  as.character(fold_val$label))
  all_pred_cv <- c(all_pred_cv, val_preds)
}

acc_results <- compute_accuracy(all_obs_cv, all_pred_cv)

# Save CV metrics
acc_summary <- data.table(
  cv_mean_oa_pct        = round(mean(cv_dt$overall_acc)*100, 2),
  cv_sd_oa_pct          = round(sd(cv_dt$overall_acc)*100, 2),
  cv_mean_kappa         = round(mean(cv_dt$kappa), 4),
  cv_sd_kappa           = round(sd(cv_dt$kappa), 4),
  cv_mean_leaf_f1       = round(mean(cv_dt$leaf_f1), 4),
  cv_mean_woody_f1      = round(mean(cv_dt$woody_f1), 4),
  agg_overall_acc_pct   = round(acc_results$overall_accuracy*100, 2),
  agg_kappa             = round(acc_results$kappa, 4),
  woody_prob_threshold  = woody_prob_threshold,
  n_leaf_raw            = n_leaf_raw,
  n_woody_raw           = n_woody_raw,
  k_folds               = k_folds
)
fwrite(acc_summary, file.path(dirname(parent_dir), "accuracy_assessment.csv"))
fwrite(cv_dt,       file.path(dirname(parent_dir), "cv_fold_metrics.csv"))
cat("Accuracy saved.\n\n")

# ── TRAIN FINAL MODEL ─────────────────────────────────────────────────────────
cat("Training final Random Forest on 100% of balanced data...\n")

rf_model_full <- randomForest(
  x = all_labelled[, ..feature_cols], y = all_labelled$label,
  ntree = 500, importance = FALSE
)

cat("OOB error:", round(rf_model_full$err.rate[500,"OOB"]*100, 2), "%\n\n")







# ── APPLY TO ALL FILES ────────────────────────────────────────────────────────
platform_folders <- list.dirs(parent_dir, recursive = FALSE)
platform_folders <- platform_folders[grepl("TLS_UAV", platform_folders)]

all_results <- list()

for (folder in platform_folders) {
  cat("\nFolder:", folder, "\n")
  las_files <- list.files(folder, pattern = "\\.las$|\\.laz$", full.names = TRUE)
  
  for (f in las_files) {
    cat("  File:", basename(f), "\n")
    
    las <- readLAS(f)
    if (is.empty(las)) { cat("  [skip – empty]\n"); next }
    
    pts <- las@data
    
    if (nrow(pts) < 20) { cat("  [skip – too few points]\n"); next }
    
    # ── 1. RF CLASSIFICATION ──────────────────────────────────────────────────
    feat_all <- extract_features(pts)
    feat_all[is.na(feat_all)] <- 0
    
    probs_all <- predict(rf_model_full, feat_all[, ..feature_cols], type = "prob")
    preds_all <- ifelse(probs_all[, "woody"] >= woody_prob_threshold, "woody", "leaf")
    
    cat(sprintf("    RF raw — leaf: %d  woody: %d  (%.1f %% leaf)\n",
                sum(preds_all == "leaf"), sum(preds_all == "woody"),
                100 * mean(preds_all == "leaf")))
    
    # ── 2. DETECT CROWN BASE Height_m ───────────────────────────────────────────
    cbh_z <- detect_cbh(pts$Z)
    cat(sprintf("    Detected CBH: %.3f m (abs Z)\n", cbh_z))
    
    # ── 3. HARD OVERRIDE — below CBH → woody ─────────────────────────────────
    below_cbh <- pts$Z <= cbh_z
    n_below   <- sum(below_cbh)
    preds_all[below_cbh] <- "woody"
    cat(sprintf("    Below-CBH override: %d points forced to woody\n", n_below))
    
    # ── 4. HARD OVERRIDE — low canopy, no-leaf zone ───────────────────────────
    # Band: CBH < Z <= CBH + low_canopy_band_m
    in_low_band <- pts$Z > cbh_z & pts$Z <= (cbh_z + low_canopy_band_m)
    
    if (sum(in_low_band) > 0) {
      # Compute raw (un-normalised) signals for thresholding in this band
      raw_density  <- compute_density(pts)          # whole cloud for context
      med_density  <- median(raw_density)
      
      raw_Rn       <- if (all(c("ReturnNumber","NumberOfReturns") %in% names(pts)))
        pts$ReturnNumber / pmax(pts$NumberOfReturns, 1)
      else rep(1, nrow(pts))
      
      has_rgb <- all(c("R","G","B") %in% names(pts))
      if (has_rgb) {
        Rn_s <- safe_norm(pts$R); Gn_s <- safe_norm(pts$G); Bn_s <- safe_norm(pts$B)
        ExG_raw <- 2*Gn_s - Rn_s - Bn_s       # normalised ExG, range ~ -1 to 1
      } else {
        ExG_raw <- rep(0, nrow(pts))
      }
      
      sparse    <- raw_density  < low_density_pct * med_density
      low_ret   <- raw_Rn       < low_return_thr
      not_green <- ExG_raw      < low_exg_thr
      
      # Rescale ExG to 0-1 for consistent thresholding
      not_green_norm <- safe_norm(ExG_raw) < low_exg_thr
      
      force_woody_low <- in_low_band & sparse & low_ret & not_green_norm
      n_low_forced    <- sum(force_woody_low)
      preds_all[force_woody_low] <- "woody"
      cat(sprintf("    Low-canopy no-leaf override: %d points forced to woody\n",
                  n_low_forced))
    }
    
    # ── 5. FINAL COUNTS ───────────────────────────────────────────────────────
    leaf_idx  <- which(preds_all == "leaf")
    woody_idx <- which(preds_all == "woody")
    
    cat(sprintf("    Final — leaf: %d  woody: %d  (%.1f %% leaf)\n",
                length(leaf_idx), length(woody_idx),
                100 * mean(preds_all == "leaf")))
    
    # ── 6. EXPORT SPLIT LAS FILES ─────────────────────────────────────────────
    out_dir   <- file.path(dirname(folder), "classified_split")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    base_name <- paste0(basename(folder), "_",
                        tools::file_path_sans_ext(basename(f)))
    
    all_idx <- seq_len(nrow(las@data))
    
    if (length(leaf_idx) > 0)
      writeLAS(filter_poi(las, all_idx %in% leaf_idx),
               file.path(out_dir, paste0(base_name, "_leaf_032526.las")))
    
    if (length(woody_idx) > 0)
      writeLAS(filter_poi(las, all_idx %in% woody_idx),
               file.path(out_dir, paste0(base_name, "_woody_032526.las")))
    
    # ── 7. STORE RESULTS ──────────────────────────────────────────────────────
    all_results[[length(all_results)+1]] <- data.table(
      folder         = basename(folder),
      file           = basename(f),
      n_total_pts    = nrow(pts),
      n_leaf_pts     = length(leaf_idx),
      n_woody_pts    = length(woody_idx),
      pct_leaf       = round(100 * length(leaf_idx) / max(nrow(pts), 1), 1),
      cbh_z          = round(cbh_z, 3),
      n_below_cbh    = n_below
    )
  }
}

# ── SUMMARY ───────────────────────────────────────────────────────────────────
results_dt <- rbindlist(all_results, fill = TRUE)
print(results_dt)
fwrite(results_dt, file.path(dirname(parent_dir), "leaf_woody_summary_032526.csv"))

saveRDS(rf_model_full, file.path(dirname(parent_dir), "rf_model_full_032526.rds"))

cat("\nDone.\n")
cat("Outputs written to:", dirname(parent_dir), "\n")
cat("  leaf_woody_summary.csv\n")
cat("  accuracy_assessment.csv\n")
cat("  cv_fold_metrics.csv\n")
cat("  rf_model_full.rds\n")
cat("  classified_split/  (leaf + woody .las files)\n")



}





















































if (FALSE) {
# ==============================================================================
# TREE METRICS — all platforms, all metrics
# ==============================================================================

library(lidR)
library(data.table)
library(VoxR)
library(dplyr)
library(stringr)
library(sf)
library(geometry)
library(sp)

# -------- SETTINGS ------------------------------------------------------------
parent_dir <- "//snre-snow.snrenet.arizona.edu/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/CH2_LiDAR/segmentedClouds/"
wood_dir   <- "//snre-snow.snrenet.arizona.edu/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/CH2_LiDAR/classified_split/wood/"
voxel_size <- 0.05
dbh_band   <- 1

# ==============================================================================
# HELPER — 2D slice-based hole filling for wood volume
# ==============================================================================
fill_voxel_holes <- function(dt_vox, voxel_size) {
  if (nrow(dt_vox) < 4) return(dt_vox)
  
  filled_list <- vector("list", length(unique(dt_vox$vz)))
  z_levels    <- sort(unique(dt_vox$vz))
  
  for (i in seq_along(z_levels)) {
    zl    <- z_levels[i]
    slice <- dt_vox[vz == zl, .(vx, vy)]
    
    if (nrow(slice) < 3) {
      filled_list[[i]] <- data.table(vx = slice$vx, vy = slice$vy, vz = zl)
      next
    }
    
    hull_idx <- tryCatch(chull(slice$vx, slice$vy), error = function(e) NULL)
    
    if (is.null(hull_idx)) {
      filled_list[[i]] <- data.table(vx = slice$vx, vy = slice$vy, vz = zl)
      next
    }
    
    hull_pts <- slice[hull_idx, ]
    x_seq    <- seq(min(slice$vx), max(slice$vx))
    y_seq    <- seq(min(slice$vy), max(slice$vy))
    grid     <- data.table(CJ(vx = x_seq, vy = y_seq))
    
    inside <- sp::point.in.polygon(
      point.x = grid$vx,
      point.y = grid$vy,
      pol.x   = hull_pts$vx,
      pol.y   = hull_pts$vy
    )
    
    filled_list[[i]] <- data.table(
      vx = grid$vx[inside >= 1],
      vy = grid$vy[inside >= 1],
      vz = zl
    )
  }
  
  rbindlist(filled_list)
}

# ==============================================================================
# STEP 1 — WOOD VOLUME (classified_split/wood)
# ==============================================================================
# ==============================================================================
# WOOD VOLUME — voxel size stabilization test then full run
# ==============================================================================

library(lidR)
library(data.table)
library(VoxR)
library(dplyr)
library(stringr)
library(ggplot2)

# ==============================================================================
# SETTINGS
# ==============================================================================
wood_dir    <- "Z:/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/CH2_LiDAR/classified_split/wood/"
test_sizes  <- c(0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.15, 0.20)
n_test_trees <- 5   # number of trees to run stabilization test on

wood_files <- list.files(wood_dir, pattern = "\\.las$|\\.laz$",
                         full.names = TRUE, recursive = FALSE)
message("Found ", length(wood_files), " wood files")

# ==============================================================================
# STEP 1 — STABILIZATION TEST on first n_test_trees files
# ==============================================================================
message("\n=== Voxel stabilization test (", n_test_trees, " trees) ===")

test_files   <- wood_files[1:min(n_test_trees, length(wood_files))]
stab_results <- list()

for (f in test_files) {
  
  filename <- tools::file_path_sans_ext(basename(f))
  message("  Testing: ", filename)
  
  las <- tryCatch(readLAS(f, select = "xyz"), error = function(e) NULL)
  if (is.null(las) || is.empty(las)) next
  
  pts   <- as.data.table(las@data)
  dt_xyz <- data.table(x = pts$X, y = pts$Y, z = pts$Z)
  
  for (vs in test_sizes) {
    vol <- tryCatch({
      vox <- VoxR::filled_voxel_cloud(dt_xyz, vs)
      nrow(vox) * vs^3
    }, error = function(e) NA_real_)
    
    stab_results[[length(stab_results) + 1]] <- data.table(
      file       = filename,
      voxel_size = vs,
      volume_m3  = vol
    )
  }
}

stab_df <- rbindlist(stab_results)

# Compute % change between successive voxel sizes per tree
stab_df <- stab_df[order(file, voxel_size)]
stab_df[, pct_change := abs(volume_m3 - data.table::shift(volume_m3)) / data.table::shift(volume_m3) * 100,
        by = file]
# Mean % change across trees at each voxel size step
change_summary <- stab_df[, .(mean_pct_change = mean(pct_change, na.rm = TRUE),
                              mean_vol        = mean(volume_m3,  na.rm = TRUE)),
                          by = voxel_size]

message("\n--- Volume change between successive voxel sizes ---")
print(change_summary)

# Find stable voxel size = first size where mean % change drops below 5%
stable_size <- change_summary[mean_pct_change < 5, voxel_size][1]

if (is.na(stable_size)) {
  message("No size stabilized below 5% — using 0.05 m as default")
  stable_size <- 0.05
} else {
  message("\n>>> Stable voxel size: ", stable_size, " m (volume change < 5%)")
}

# ==============================================================================
# PLOT — stabilization curve
# ==============================================================================
windows(width = 10, Height_m = 6)
ggplot(stab_df, aes(x = voxel_size, y = volume_m3,
                    group = file, color = file)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 2) +
  geom_vline(xintercept = stable_size, linetype = "dashed",
             color = "red", linewidth = 1) +
  annotate("text", x = stable_size + 0.005, y = max(stab_df$volume_m3, na.rm = TRUE),
           label = paste0("stable: ", stable_size, " m"),
           color = "red", hjust = 0, size = 4) +
  scale_x_continuous(breaks = test_sizes) +
  labs(
    title    = "Wood volume vs voxel size (filled voxel cloud)",
    subtitle = paste0("Red dashed = first stable size (< 5% change) = ", stable_size, " m"),
    x        = "Voxel size (m)",
    y        = "Filled voxel volume (m³)",
    color    = "Tree"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 8),
    plot.title      = element_text(face = "bold")
  )

# ==============================================================================
# STEP 2 — FULL RUN at stable voxel size (UPDATED for plotting)
# ==============================================================================
library(data.table)
library(stringr)
library(dplyr)
library(VoxR)
library(alphashape3d)
library(lidR)
library(rgl)

n_trees      <- length(wood_files)
wood_results <- vector("list", n_trees)

cat(sprintf("\n╔══════════════════════════════════════╗\n"))
cat(sprintf("║   Wood Volume Pipeline Starting...   ║\n"))
cat(sprintf("╚══════════════════════════════════════╝\n"))
cat(sprintf("  Trees to process : %d\n", n_trees))
cat(sprintf("  Voxel size       : %.4f m\n\n", stable_size))

# -------------------------------
# RANDOMLY SELECT 2 EXAMPLE TREES
# -------------------------------
set.seed(42)  # reproducibility
example_trees <- sample(seq_along(wood_files), 2)

# ================================
# MAIN LOOP
# ================================
for (i in seq_along(wood_files)) {
  
  progress_bar(i - 1, n_trees)
  
  # --- reset per-tree outputs ---
  wood_vol       <- NA_real_
  occupied_vol   <- NA_real_
  supported_vol  <- NA_real_
  fill_fraction  <- NA_real_
  hull_vol       <- NA_real_
  site_id        <- NA_character_
  platform_id    <- NA_character_
  tree_id        <- NA_character_
  
  tryCatch({
    
    # ================================
    # PARSE IDS
    # ================================
    fname       <- basename(wood_files[i])
    site_id     <- str_match(fname, "^[^_]+_[^_]+_([^_]+)_")[, 2]
    platform_id <- str_match(fname, "^([^_]+_[^_]+)_")[, 2]
    tree_id     <- str_match(fname, "tree(\\d+)")[, 2]
    
    cat(sprintf("\n  ▶ Tree %d/%d  [site: %s | platform: %s | treeID: %s]\n",
                i, n_trees, site_id, platform_id, tree_id))
    
    # ================================
    # 0. LOAD POINT CLOUD
    # ================================
    cat("    [0/4] Loading point cloud... ")
    las <- readLAS(wood_files[i])
    pts <- as.data.table(las@data)[, .(X, Y, Z)]
    cat(sprintf("done  (%d points)\n", nrow(pts)))
    
    # ================================
    # 1. VOXEL (OCCUPIED)
    # ================================
    cat("    [1/4] Building occupied voxels... ")
    
    dt_xyz <- data.table(x = pts$X, y = pts$Y, z = pts$Z)
    voxels <- VoxR::vox(dt_xyz, res = stable_size)
    
    vox_dt <- data.table(
      ix = round(voxels$x / stable_size),
      iy = round(voxels$y / stable_size),
      iz = round(voxels$z / stable_size)
    )
    
    vox_dt <- unique(vox_dt)
    vox_dt[, occupied := TRUE]
    setkey(vox_dt, ix, iy, iz)
    
    cat(sprintf("done  (%d occupied voxels)\n", nrow(vox_dt)))
    
    # ================================
    # 2. SUPPORTED VOXELS
    # ================================
    cat("    [2/4] Finding neighbor-supported voxels... ")
    
    neighbors <- CJ(dx = -1:1, dy = -1:1, dz = -1:1)
    neighbors <- neighbors[!(dx == 0 & dy == 0 & dz == 0)]
    
    candidates <- unique(rbindlist(lapply(1:nrow(neighbors), function(j) {
      vox_dt[, .(ix = ix + neighbors$dx[j],
                 iy = iy + neighbors$dy[j],
                 iz = iz + neighbors$dz[j])]
    })))
    
    candidates <- candidates[!vox_dt, on = .(ix, iy, iz)]
    
    occ_expanded <- rbindlist(lapply(1:nrow(neighbors), function(j) {
      vox_dt[, .(ix = ix + neighbors$dx[j],
                 iy = iy + neighbors$dy[j],
                 iz = iz + neighbors$dz[j])]
    }))
    
    hit_counts <- occ_expanded[
      candidates, on = .(ix, iy, iz), nomatch = 0L
    ][, .N, by = .(ix, iy, iz)]
    
    candidates <- hit_counts[candidates, on = .(ix, iy, iz)]
    candidates[is.na(N), N := 0L]
    setnames(candidates, "N", "n_neighbors")
    
    threshold <- 4
    supported <- candidates[n_neighbors >= threshold]
    
    cat(sprintf("done  (%d supported voxels, threshold = %d)\n",
                nrow(supported), threshold))
    
 
    # ================================
    # 2.5 VISUALIZE (ONLY SELECT TREES)
    # ================================
    if (i %in% example_trees) {
      
      cat("    [viz] Plotting fused voxels...\n")
      
      # ---- COMBINE (FUSE) VOXELS ----
      all_vox <- rbindlist(list(
        vox_dt[, .(ix, iy, iz)],
        supported[, .(ix, iy, iz)]
      ))
      
      # Remove duplicates
      all_vox <- unique(all_vox)
      
      # Convert to coordinates
      fused_coords <- all_vox[, .(
        x = ix * stable_size,
        y = iy * stable_size,
        z = iz * stable_size
      )]
      
      # Create VoxR object
      vox_fused <- VoxR::vox(fused_coords, res = stable_size)
      
      # Plot
      open3d()
      plot_voxels(vox_fused)
    }
    
    # ================================
    # 3. VOXEL VOLUME
    # ================================
    cat("    [3/4] Computing voxel volume... ")
    
    voxel_volume <- stable_size^3
    
    occupied_vol  <- nrow(vox_dt) * voxel_volume
    supported_vol <- nrow(supported) * voxel_volume
    wood_vol      <- occupied_vol + supported_vol
    
    fill_fraction <- ifelse(wood_vol > 0,
                            supported_vol / wood_vol,
                            NA_real_)
    
    cat(sprintf("done  (%.6f m³)\n", wood_vol))
    cat(sprintf("         ↳ occupied:  %.6f m³\n", occupied_vol))
    cat(sprintf("         ↳ supported: %.6f m³ (%.2f%% fill)\n",
                supported_vol, 100 * fill_fraction))
    
    cat(sprintf("    ✔ Tree complete — wood: %.6f m³ | hull: %.6f m³ | ratio: %.3f\n",
                wood_vol, hull_vol, ifelse(hull_vol>0, wood_vol/hull_vol, NA_real_)))
    
  }, error = function(e) {
    cat(sprintf("\n    ✘ Processing failed: %s\n", e$message))
    cat(sprintf("    ✘ Call: %s\n", deparse(e$call)))
  })
  
  # ================================
  # STORE RESULTS
  # ================================
  wood_results[[i]] <- data.table(
    site              = site_id,
    Platform          = platform_id,
    treeID            = tree_id,
    wood_vol_m3       = wood_vol,
    occupied_vol_m3   = occupied_vol,
    supported_vol_m3  = supported_vol,
    fill_fraction     = fill_fraction,
    wood_voxel_size   = stable_size
  )
  
  progress_bar(i, n_trees)
}

# ================================
# COMBINE & SUMMARISE
# ================================
cat("\n══════════════════════════════════════\n")
cat("  Combining results...\n")

wood_df <- rbindlist(wood_results, fill = TRUE) %>%
  mutate(
    treeID = as.integer(treeID),
    Platform = gsub("_MBC$", "", Platform),
    
    vol_ratio = ifelse(hull_vol > 0,
                       wood_vol_m3 / hull_vol,
                       NA_real_),
    
    vol_gap = hull_vol - wood_vol_m3,
    
    occupancy_ratio = ifelse(wood_vol_m3 > 0,
                             occupied_vol_m3 / wood_vol_m3,
                             NA_real_)
  ) %>%
  arrange(treeID)

cat(sprintf("\n✔ Done.  %d trees processed at %.4f m voxel size\n\n",
            nrow(wood_df), stable_size))

print(wood_df)
# ==============================================================================
# STEP 2 — FULL TREE METRICS (segmentedClouds)
# ==============================================================================
message("\n=== Computing tree metrics ===")

platform_folders <- list.dirs(parent_dir, recursive = FALSE)
all_results      <- list()

for (folder in platform_folders) {
  
  las_files <- list.files(folder, pattern = "\\.las$", full.names = TRUE)
  
  for (f in las_files) {
    
    # ── Parse filename ─────────────────────────────────────────────────────────
    filename    <- tools::file_path_sans_ext(basename(f))
    parts       <- str_split(filename, "_", simplify = TRUE)
    plot_id     <- parts[1]
    tree_id     <- gsub("tree", "", parts[length(parts)], ignore.case = TRUE)
    platform_id <- paste(parts[2:(length(parts) - 1)], collapse = "_")
    
    # ── Load LAS ───────────────────────────────────────────────────────────────
    las <- readLAS(f)
    if (is.empty(las)) next
    
    pts <- las@data
    if (nrow(pts) < 10) next
    
    # ── Height_m ─────────────────────────────────────────────────────────────────
    z_min  <- min(pts$Z, na.rm = TRUE)
    z_max  <- quantile(pts$Z, 0.995, na.rm = TRUE)
    Height_m <- z_max - z_min
    mid_z  <- z_min + Height_m / 2
    
    # ── Crown Area (2D convex hull on XY) ──────────────────────────────────────
    crown_area <- NA_real_
    if (nrow(pts) >= 3) {
      tryCatch({
        coords2d   <- as.data.frame(pts[, c("X", "Y")])
        sf_pts     <- st_as_sf(coords2d, coords = c("X", "Y"))
        poly       <- st_convex_hull(st_union(sf_pts))
        crown_area <- as.numeric(st_area(poly))
      }, error = function(e) message("Crown area error in ", f, ": ", e$message))
    }
    
    # ── Voxel Volume (filled voxels × voxel_size³) ────────────────────────────
    total_volume <- NA_real_
    tryCatch({
      dt_xyz <- data.table::data.table(pts[, c("X", "Y", "Z")])
      data.table::setnames(dt_xyz, old = c("X", "Y", "Z"), new = c("x", "y", "z"))
      voxels       <- VoxR::filled_voxel_cloud(dt_xyz, voxel_size)
      total_volume <- nrow(voxels) * (voxel_size^3)
    }, error = function(e) message("Voxel volume error in ", f, ": ", e$message))
    
    # ── Top / Bottom Voxel Density ────────────────────────────────────────────
    bottom_density <- NA_real_
    top_density    <- NA_real_
    
    tryCatch({
      pts_bot <- pts[pts$Z <= mid_z, ]
      pts_top <- pts[pts$Z >  mid_z, ]
      
      xrange <- max(pts$X) - min(pts$X)
      yrange <- max(pts$Y) - min(pts$Y)
      
      possible_half <- (xrange / voxel_size) *
        (yrange / voxel_size) *
        ((Height_m / 2) / voxel_size)
      
      if (possible_half > 0 && nrow(pts_bot) >= 3) {
        dt_bot         <- data.table::data.table(x = pts_bot$X, y = pts_bot$Y, z = pts_bot$Z)
        vox_bot        <- vox(dt_bot, res = voxel_size)
        bottom_density <- nrow(vox_bot) / possible_half
      }
      
      if (possible_half > 0 && nrow(pts_top) >= 3) {
        dt_top      <- data.table::data.table(x = pts_top$X, y = pts_top$Y, z = pts_top$Z)
        vox_top     <- vox(dt_top, res = voxel_size)
        top_density <- nrow(vox_top) / possible_half
      }
    }, error = function(e) {})
    
    # ── Leaf / Woody Voxel Volume (density + intensity) ───────────────────────
    leaf_voxel_volume  <- NA_real_
    woody_voxel_volume <- NA_real_
    
    tryCatch({
      pts_canopy <- pts[pts$Z > z_min + (Height_m * 3), ]
      
      if (nrow(pts_canopy) >= 10) {
        
        dt_canopy <- data.table::data.table(
          x         = pts_canopy$X,
          y         = pts_canopy$Y,
          z         = pts_canopy$Z,
          intensity = if ("Intensity" %in% names(pts_canopy)) pts_canopy$Intensity else NA_real_
        )
        
        dt_canopy[, `:=`(
          vox_x = floor((x - min(x)) / voxel_size),
          vox_y = floor((y - min(y)) / voxel_size),
          vox_z = floor((z - min(z)) / voxel_size)
        )]
        
        vox_stats <- dt_canopy[, .(
          N              = .N,
          mean_intensity = mean(intensity, na.rm = TRUE)
        ), by = .(vox_x, vox_y, vox_z)]
        
        density_thresh   <- quantile(vox_stats$N,              0.5, na.rm = TRUE)
        intensity_thresh <- quantile(vox_stats$mean_intensity, 0.5, na.rm = TRUE)
        
        vox_stats[, is_foliage := N              >= density_thresh &
                    mean_intensity <= intensity_thresh]
        
        leaf_voxel_volume  <- nrow(vox_stats[is_foliage == TRUE])  * (voxel_size^3)
        woody_voxel_volume <- nrow(vox_stats[is_foliage == FALSE]) * (voxel_size^3)
      }
      
    }, error = function(e) message("Leaf voxel error in ", f, ": ", e$message))
    
    # ── Convex Hull Volume (3D) ───────────────────────────────────────────────
    hull_volume <- NA_real_
    tryCatch({
      pts_mat <- as.matrix(pts[, c("X", "Y", "Z")])
      pts_mat <- unique(pts_mat)
      if (nrow(pts_mat) >= 4) {
        hull        <- geometry::convhulln(pts_mat, options = "FA")
        hull_volume <- hull$vol
      }
    }, error = function(e) message("Hull3D volume error in ", f, ": ", e$message))
    
    # ── DBH ───────────────────────────────────────────────────────────────────
    dbh   <- NA_real_
    slice <- pts[pts$Z > (1.3 - dbh_band) &
                   pts$Z < (1.3 + dbh_band), ]
    
    if (nrow(slice) > 30) {
      xm  <- mean(slice$X)
      ym  <- mean(slice$Y)
      dbh <- 2 * mean(sqrt((slice$X - xm)^2 + (slice$Y - ym)^2)) * 100
    }
    
    # ── Collect ───────────────────────────────────────────────────────────────
    all_results[[length(all_results) + 1]] <- data.frame(
      site                 = plot_id,
      Platform             = platform_id,
      treeID               = tree_id,
      Height_m_m             = Height_m,
      CrownArea_m2         = crown_area,
      VoxelVolume_m3       = total_volume,
      HullVolume_m3        = hull_volume,
      BottomVoxelDensity   = bottom_density,
      TopVoxelDensity      = top_density,
      DBH_cm               = dbh
    )
  }
}

# ==============================================================================
# STEP 3 — JOIN & WRITE
# ==============================================================================
results_df <- dplyr::bind_rows(all_results) %>%
  mutate(treeID = as.double(treeID)) %>%
  # Join field measurements
  left_join(all_trees_bigelow %>% select(treeID, site, species, Height_m, DBH),
            by = c("treeID" = "treeID", "site" = "site")) %>%
  # Join wood volumes from classified_split/wood
  left_join(wood_df,
            by = c("site" = "site", "treeID" = "treeID", "Platform" = "Platform"))



}



results_df <- read_csv(
  "\\\\snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/all_platform_tree_metrics_05mvox_033126.csv"
)
final <- results_df %>%
  arrange(site, treeID) %>%
  mutate(tree_index = as.integer(factor(paste(site, treeID))))

setorder(final, site, treeID, Platform)

#write.csv(final, "all_platform_tree_metrics_05mvox_033126.csv", row.names = FALSE)
message("Done. Rows written: ", nrow(final))































































































































































# ==============================================================================
# STEP 4 — PLOTS
# ==============================================================================
platform_colors <- c(
  "nadir"         = "#2196F3",
  "nadir_oblique" = "#FF9800",
  "TLS"           = "#4CAF50",
  "TLS_UAV"       = "red"
)

# ── Main metrics ──────────────────────────────────────────────────────────────
final_long <- final %>%
  select(site, treeID, tree_index, Platform,
         Height_m_m, CrownArea_m2, VoxelVolume_m3) %>%
  pivot_longer(cols = c(Height_m_m, CrownArea_m2, VoxelVolume_m3),
               names_to  = "metric",
               values_to = "value")



tls_ref <- final %>%
  filter(Platform == "TLS") %>%
  select(tree_index, tls_volume = VoxelVolume_m3)

diff_df <- final %>%
  left_join(tls_ref, by = "tree_index") %>%
  mutate(volume_diff_tls = VoxelVolume_m3 - tls_volume)
diff_summary <- diff_df %>%
  filter(Platform != "TLS") %>%  # exclude TLS vs itself
  group_by(Platform) %>%
  summarise(
    mean_diff = mean(volume_diff_tls, na.rm = TRUE),
    sd_diff   = sd(volume_diff_tls, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      Platform, ": ",
      round(mean_diff, 2), " ± ", round(sd_diff, 2), " m³"
    )
  )
# get top of voxel volume panel
y_max <- final_long %>%
  filter(metric == "VoxelVolume_m3") %>%
  summarise(max_val = max(value, na.rm = TRUE)) %>%
  pull(max_val)

diff_summary <- diff_summary %>%
  mutate(
    x = 1,  # left side of plot
    y = y_max * seq(0.95, 0.75, length.out = n())
  )


windows(width = 18, Height_m = 20)
ggplot(final_long, aes(x = factor(tree_index), y = value,
                       group = Platform, color = Platform)) +
  geom_line(alpha = 0.8, linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_color_manual(values = platform_colors) +
  labs(title = "Per-tree metrics across platforms",
       x = "Tree Index", y = "Value", color = "Platform") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    strip.text  = element_text(face = "bold", size = 13),
    plot.title  = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 12)
  )+ geom_text(
    data = diff_summary,
    aes(x = x, y = y, label = label, color = Platform),
    inherit.aes = FALSE,
    hjust = 0,
    size = 5,
    fontface = "bold"
  )






















# ── Voxel density ─────────────────────────────────────────────────────────────
density_long <- results_df %>%
  select(treeID, site, Platform, BottomVoxelDensity, TopVoxelDensity) %>%
  mutate(tree_index = as.integer(factor(paste(site, treeID)))) %>%
  pivot_longer(cols = c(BottomVoxelDensity, TopVoxelDensity),
               names_to  = "metric",
               values_to = "value")
density_diff <- results_df %>% select(treeID, site, Platform, BottomVoxelDensity, TopVoxelDensity) %>% 
  left_join( results_df %>% filter(Platform == "TLS") %>% 
               select(treeID, site, TLS_bottom = BottomVoxelDensity, TLS_top = TopVoxelDensity), by = c("treeID", "site") ) %>% 
  filter(Platform != "TLS") %>% 
  pivot_longer( cols = c(BottomVoxelDensity, TopVoxelDensity), names_to = "position", values_to = "lidar_density" ) %>% 
  mutate( tls_density = ifelse(position == "BottomVoxelDensity", TLS_bottom, TLS_top), position = dplyr::recode(position, "BottomVoxelDensity" = "Bottom", "TopVoxelDensity" = "Top") )
# ── Convert to percentage density (×100) ─────────────────────────────────────
density_diff <- density_diff %>%
  mutate(
    lidar_density = lidar_density * 100,
    tls_density   = tls_density   * 100
  )

# ── Compute paired t-tests + mean & SD (in % units) ──────────────────────────
ttest_labels <- density_diff %>%
  group_by(Platform, position) %>%
  summarise(
    p_val     = tryCatch(
      t.test(lidar_density, tls_density, paired = TRUE)$p.value,
      error = function(e) NA_real_
    ),
    mean_diff = mean(lidar_density - tls_density, na.rm = TRUE),
    sd_diff   = sd(lidar_density   - tls_density, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      is.na(p_val)  ~ "p = NA",
      p_val < 0.001 ~ "p < 0.001***",
      p_val < 0.01  ~ paste0("p = ", formatC(p_val, digits = 3, format = "f"), "**"),
      p_val < 0.05  ~ paste0("p = ", formatC(p_val, digits = 3, format = "f"), "*"),
      TRUE          ~ paste0("p = ", formatC(p_val, digits = 3, format = "f"), " ns")
    ),
    stats_label = paste0(
      Platform, "\n",
      "  ", p_label, "\n",
      "  \u0394mean = ", formatC(mean_diff, digits = 3, format = "f"), "%\n",
      "  \u0394SD = ",   formatC(sd_diff,   digits = 3, format = "f"), "%"
    )
  )

# ── Annotation positions: BOTTOM-RIGHT corner to avoid points ─────────────────
# Place labels at bottom-right, stacking upward per platform
# ── Annotation positions: TOP-LEFT corner to avoid points ─────────────────────
# Place labels at top-left, stacking downward per platform
ttest_annotations <- ttest_labels %>%
  group_by(position) %>%
  mutate(row_idx = row_number()) %>%        # 1 = top-most label
  ungroup() %>%
  left_join(
    density_diff %>%
      group_by(position) %>%
      summarise(
        x_min = min(tls_density,   na.rm = TRUE),
        y_max = max(lidar_density, na.rm = TRUE),
        y_rng = diff(range(lidar_density, na.rm = TRUE)),
        .groups = "drop"
      ),
    by = "position"
  ) %>%
  mutate(
    x_ann = x_min,                                      # left-align
    y_ann = y_max - (row_idx - 1) * y_rng * 0.18       # stack downward
  )
# ── Plot ──────────────────────────────────────────────────────────────────────
windows(width = 16, Height_m = 12)

ggplot(density_diff, aes(x = tls_density, y = lidar_density,
                         color = Platform, shape = position)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, aes(group = Platform)) +
  geom_text(
    data        = ttest_annotations,
    aes(x = x_ann, y = y_ann, label = stats_label, color = Platform),
    hjust       = 0,          # left-justify so text grows rightward from x_min
    vjust       = 1,          # anchor top of text block at y_ann
    size        = 3.5,
    lineHeight_m  = 1.15,
    fontface    = "italic",
    inherit.aes = FALSE
  )  +
  facet_wrap(~ position, scales = "free") +
  scale_color_manual(values = platform_colors) +
  scale_shape_manual(values = c("Bottom" = 16, "Top" = 17)) +
  scale_x_continuous(labels = scales::label_number(suffix = "%", accuracy = 0.01)) +
  scale_y_continuous(labels = scales::label_number(suffix = "%", accuracy = 0.01)) +
  labs(
    title    = "Voxel density: TLS vs other platforms",
    subtitle = "Dashed line = 1:1 | \u0394 = Platform \u2212 TLS (%) | Paired t-test vs TLS",
    x        = "TLS Voxel Density (%)", y = "Platform Voxel Density (%)",
    color    = "Platform", shape = "Crown Position"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text  = element_text(face = "bold", size = 14),
    plot.title  = element_text(face = "bold", size = 15),
    axis.text   = element_text(size = 12),
    legend.text = element_text(size = 12)
  )











































































# ── Pivot volume metrics to long then wide for joining ────────────────────────
volume_long <- results_df %>%
  select(treeID, site, Platform, VoxelVolume_m3, HullVolume_m3,wood_vol_m3) %>%
  rename(volume = VoxelVolume_m3, volume_hull = HullVolume_m3, wood_volume = wood_vol_m3)%>%
  
  mutate(tree_index = row_number())
# ── Species metrics with biomass estimates ────────────────────────────────────
species_metrics <- results_df %>%
  left_join(volume_long) %>%
  mutate(
    # Voxel volume biomass (m3 × wood density kg/m3 → metric ton)
    biomass_allometry_volume = case_when(
      species == "ponderosa"    ~ (wood_vol_m3     * 0.5308) ,
      species == "menziesii"   ~ (wood_vol_m3      * 0.5327) ,
      species == "strobiformis"~ (wood_vol_m3    * 0.5299) ,
      TRUE ~ NA_real_
    ),
    # DBH allometric biomass (log-linear, result in kg → metric ton)
    biomass_allometry_ground = case_when(
      species == "ponderosa"    ~ (exp(-2.6177 + 2.4638 * log(DBH))/1000) ,
      species == "menziesii"   ~ (exp(-2.4623 + 2.4852 * log(DBH))/1000) ,
      species == "strobiformis"~ (exp(-2.6177 + 2.4638 * log(DBH))/1000) ,
      TRUE ~ NA_real_
    ),
    # Crown area × Height_m allometric biomass
    biomass_allometry_crownHeight_m = case_when(
      species == "ponderosa"    ~ (0.0776  * (CrownArea_m2^0.934)  * (Height_m_m^1.018))/1000,
      species == "menziesii"   ~ (0.0474  * (CrownArea_m2^0.879)  * (Height_m_m^1.119))/1000,
      species == "strobiformis"~ (0.0618 * (CrownArea_m2^0.950)  * (Height_m_m^0.988))/1000,
      TRUE ~ NA_real_
    )
  )

#write.csv(species_metrics,"species_metrics_05mvox_033126.csv")
getwd()




library(ggplot2)
library(dplyr)
library(tidyr)

# Filter for TLS_UAV platform only
species_metrics_tls_uav <- species_metrics %>%
  filter(Platform == "TLS_UAV") %>%
  mutate(
    diff_volume_vs_ground = biomass_allometry_volume - biomass_allometry_ground,
    diff_crown_vs_ground  = biomass_allometry_crownHeight_m - biomass_allometry_ground
  ) %>%
  select(diff_volume_vs_ground, diff_crown_vs_ground)

# Pivot to long format for ggplot
diff_long <- species_metrics_tls_uav %>%
  pivot_longer(
    everything(),
    names_to = "method",
    values_to = "diff"
  ) 
# Compute means for labeling
means <- diff_long %>%
  group_by(method) %>%
  summarize(mean_diff = mean(diff, na.rm = TRUE))

# Boxplot with mean labels
windows()
ggplot(diff_long, aes(x = method, y = diff, fill = method)) +
  geom_boxplot(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_text(
    data = means,
    aes(x = method, y = mean_diff, label = round(mean_diff, 2)),
    vjust = -0.5,
    fontface = "bold",
    color = "black"
  ) +
  labs(
    x = "Method",
    y = "Difference from Ground Allometry (t)",
    title = "Comparison of Biomass Estimates vs Ground Allometry (TLS_UAV)"
  ) +
  theme_minimal(base_size = 16) +             # base font size up from default 11
  theme(
    legend.position = "none",
    plot.title  = element_text(face = "bold", size = 18),
    axis.title  = element_text(size = 16),
    axis.text   = element_text(size = 14)
  )










library(ggplot2)
library(dplyr)
library(tidyr)

# ── Filter TLS_UAV and compute differences ────────────────────────────────────
species_metrics_tls_uav <- species_metrics %>%
  filter(Platform == "TLS_UAV") %>%
  mutate(
    diff_volume_vs_ground = biomass_allometry_volume - biomass_allometry_ground,
    diff_crown_vs_ground  = biomass_allometry_crownHeight_m - biomass_allometry_ground
  )

# ── FIT CORRECTION MODEL: crown-Height_m ~ voxel (where both exist) ─────────────
correction_df <- species_metrics_tls_uav %>%
  filter(!is.na(biomass_allometry_volume), !is.na(biomass_allometry_crownHeight_m))

correction_model <- lm(biomass_allometry_volume ~ biomass_allometry_crownHeight_m, 
                       data = correction_df)

slope     <- round(coef(correction_model)[2], 3)
intercept <- round(coef(correction_model)[1], 3)
r2        <- round(summary(correction_model)$r.squared, 3)
rmse      <- round(sqrt(mean(residuals(correction_model)^2)), 3)

cat("Correction model: Biomass_voxel =", intercept, "+", slope, 
    "× Biomass_crown\n")
cat("R² =", r2, "| RMSE =", rmse, "t\n")

# ── Apply correction to get calibrated crown estimates ────────────────────────
species_metrics_tls_uav <- species_metrics_tls_uav %>%
  mutate(
    biomass_corrected_crown = intercept + slope * biomass_allometry_crownHeight_m,
    diff_corrected_vs_ground = biomass_corrected_crown - biomass_allometry_ground
  )

# ── Pivot to long format ──────────────────────────────────────────────────────
diff_long <- species_metrics_tls_uav %>%
  select(diff_volume_vs_ground, diff_crown_vs_ground, diff_corrected_vs_ground) %>%
  pivot_longer(everything(), names_to = "method", values_to = "diff")

# ── Means for labels ──────────────────────────────────────────────────────────
means <- diff_long %>%
  group_by(method) %>%
  summarize(mean_diff = mean(diff, na.rm = TRUE))

# ── Equation label for plot ───────────────────────────────────────────────────
eq_label <- paste0(
  "Correction model (voxel plots only):\n",
  "Biomass_voxel = ", intercept, " + ", slope, " × Biomass_crown\n",
  "R² = ", r2, "  |  RMSE = ", rmse, " t"
)

# ── PLOT ──────────────────────────────────────────────────────────────────────
windows()
ggplot(diff_long, aes(x = method, y = diff, fill = method)) +
  geom_boxplot(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_text(
    data = means,
    aes(x = method, y = mean_diff, label = round(mean_diff, 2)),
    vjust = -0.5,
    fontface = "bold",
    color = "black"
  ) +
  # ── Correction equation annotation ─────────────────────────────────────────
  annotate("label",
           x = 0.6, y = Inf,
           label = eq_label,
           hjust = 0, vjust = 1.1,
           size = 4, fill = "white", color = "black",
           label.padding = unit(0.4, "lines"),
           label.border  = unit(0.3, "lines")
  ) +
  scale_x_discrete(labels = c(
    diff_corrected_vs_ground = "diff_corrected\n(crown × voxel model)",
    diff_crown_vs_ground     = "diff_crown\nvs ground",
    diff_volume_vs_ground    = "diff_volume\nvs ground"
  )) +
  scale_fill_manual(values = c(
    diff_corrected_vs_ground = "#9b59b6",
    diff_crown_vs_ground     = "#e07070",
    diff_volume_vs_ground    = "#5bbfbf"
  )) +
  labs(
    x     = "Method",
    y     = "Difference from Ground Allometry (t)",
    title = "Comparison of Biomass Estimates vs Ground Allometry (TLS_UAV)",
    caption = paste0(
      "Corrected crown estimates apply the voxel-derived linear model to plots lacking voxel data.\n",
      "Correction fitted on n = ", nrow(correction_df), " trees with both voxel and crown-Height_m estimates."
    )
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title    = element_text(face = "bold", size = 18),
    plot.caption  = element_text(size = 11, color = "gray40", hjust = 0),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 13)
  )






























library(ggplot2)
library(dplyr)
library(tidyr)

# ── Filter TLS_UAV ────────────────────────────────────────────────────────────
# ── Clean Species: convert literal "N/A" strings to real NA ──────────────────
species_metrics_tls_uav <- species_metrics_tls_uav %>%
  mutate(species = na_if(species, "N/A"))

# ── FIT LOG-LOG CORRECTION MODEL ─────────────────────────────────────────────
correction_df <- species_metrics_tls_uav %>%
  filter(!is.na(biomass_allometry_volume), !is.na(Height_m),
         !is.na(CrownArea_m2), !is.na(species),
         biomass_allometry_volume > 0, Height_m > 0, CrownArea_m2 > 0)

correction_model <- lm(log(biomass_allometry_volume) ~ species + log(Height_m) * log(CrownArea_m2),
                       data = correction_df)

r2   <- round(summary(correction_model)$r.squared, 3)
rmse_log <- round(sqrt(mean(residuals(correction_model)^2)), 3)
cat("R² =", r2, "| RMSE (log scale) =", rmse_log, "| n =", nrow(correction_df), "\n")
print(summary(correction_model))

# ── Apply correction: predict in log space, back-transform ───────────────────
# Smearing correction (Duan 1983) to correct for retransformation bias
smearing <- mean(exp(residuals(correction_model)))
cat("Smearing factor:", round(smearing, 4), "\n")

species_metrics_tls_uav <- species_metrics_tls_uav %>%
  mutate(
    biomass_corrected = if_else(
      !is.na(species) & !is.na(Height_m) & !is.na(CrownArea_m2) &
        Height_m > 0 & CrownArea_m2 > 0,
      smearing * exp(predict(correction_model, newdata = pick(everything()))),
      NA_real_
    ),
    diff_volume_vs_ground    = biomass_allometry_volume - biomass_allometry_ground,
    diff_corrected_vs_ground = biomass_corrected        - biomass_allometry_ground
  )

# ── RMSE on original scale ────────────────────────────────────────────────────
rmse_orig <- correction_df %>%
  mutate(
    pred = smearing * exp(predict(correction_model)),
    sq_err = (biomass_allometry_volume - pred)^2
  ) %>%
  summarize(rmse = round(sqrt(mean(sq_err)), 3)) %>%
  pull(rmse)

cat("RMSE (original scale) =", rmse_orig, "t\n")

# ── Pivot to long format ──────────────────────────────────────────────────────
diff_long <- species_metrics_tls_uav %>%
  select(Species, diff_volume_vs_ground, diff_corrected_vs_ground) %>%
  pivot_longer(-Species, names_to = "method", values_to = "diff")

# ── Means per species + method ────────────────────────────────────────────────
means <- diff_long %>%
  group_by(Species, method) %>%
  summarize(mean_diff = mean(diff, na.rm = TRUE), .groups = "drop")

# ── Equation label ────────────────────────────────────────────────────────────
eq_label <- paste0(
  "Model: log(Biomass_voxel) ~ Species + log(Height_m) * log(CrownArea_m2)\n",
  "R² = ", r2, "  |  RMSE = ", rmse_orig, " t (original scale)",
  "  |  n = ", nrow(correction_df), "\n",
  "Back-transformed with Duan smearing factor = ", round(smearing, 4)
)

# ── PLOT ──────────────────────────────────────────────────────────────────────
windows()
ggplot(diff_long, aes(x = method, y = diff, fill = method)) +
  geom_boxplot(alpha = 0.6, outlier.shape = 21, outlier.size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_text(
    data = means,
    aes(x = method, y = mean_diff, label = round(mean_diff, 2)),
    vjust = -0.5,
    fontface = "bold",
    size = 3.5,
    color = "black"
  ) +
  annotate("label",
           x = 0.6, y = Inf,
           label = eq_label,
           hjust = 0, vjust = 1.1,
           size = 3.2, fill = "white", color = "black",
           label.padding = unit(0.3, "lines"),
           label.border  = unit(0.3, "lines")
  ) +
  scale_x_discrete(labels = c(
    diff_corrected_vs_ground = "Corrected\n(log-log Species + H×CA)",
    diff_volume_vs_ground    = "Voxel Volume\nvs Ground"
  )) +
  scale_fill_manual(values = c(
    diff_corrected_vs_ground = "#9b59b6",
    diff_volume_vs_ground    = "#5bbfbf"
  )) +
  facet_wrap(~ Species, scales = "free_y") +
  labs(
    x       = "Method",
    y       = "Difference from Ground Allometry (t)",
    title   = "Biomass Estimate Differences by Species — Log-Log Model (TLS_UAV)",
    caption = paste0(
      "Log-log model with Duan smearing back-transformation.\n",
      "Fitted on n = ", nrow(correction_df), " trees with voxel biomass > 0, Height_m > 0, CrownArea_m2 > 0."
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", size = 16),
    plot.caption     = element_text(size = 10, color = "gray40", hjust = 0),
    strip.text       = element_text(face = "bold", size = 12),
    axis.title       = element_text(size = 13),
    axis.text        = element_text(size = 11),
    axis.text.x      = element_text(size = 9)
  )





















library(ggplot2)
library(dplyr)
library(tidyr)

# Filter for TLS_UAV platform only
species_metrics_tls_uav <- species_metrics %>%
  filter(Platform == "TLS_UAV") %>%
  select(biomass_allometry_ground,
         biomass_allometry_volume)
biomass_long <- species_metrics_tls_uav %>%
  mutate(Tree = row_number()) %>%  # sequential x-axis
  pivot_longer(
    cols = c(biomass_allometry_ground, biomass_allometry_volume),
    names_to = "method",
    values_to = "biomass"
  )

ggplot(biomass_long, aes(x = Tree, y = biomass, color = method, group = method)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.5) +
  labs(
    x = "Tree",
    y = "Biomass (t)",
    title = "Biomass Estimates by Method (TLS_UAV)"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")








































####BIOMASS RECONSTRUCTION########
##LiDAR TREES BIOMASS current date 
#subset the DBH of the biomassI have available for Sap Flow
library(dplyr)
library(readr)

lidar_biomass <- read_csv(
  "\\\\snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/species_metrics_05mvox_033126.csv"
)

dbh_tol <- 100  # cm tolerance

lidar_biomass_subset <- lidar_biomass %>%
  mutate(
    Species = case_when(
      species == "ponderosa"    ~ "PIPO",
      species == "strobiformis" ~ "PISF",
      species == "menziesii"    ~ "PSME",
      TRUE                      ~ species
    )
  ) %>%
  rowwise() %>%
  filter(
    Species %in% names(dbh_list) &&
      any(abs(DBH - dbh_list[[Species]]) <= dbh_tol)
  ) %>%
  ungroup()




library(ggplot2)
library(broom)

# ── Subset by platform ────────────────────────────────────────────────────────
# Change "TLS" to whatever platform value(s) you want to keep
lidar_biomass_subset <- lidar_biomass_subset %>%
  filter(Platform == "TLS_UAV")   # e.g. "TLS", "ALS", "UAS" — adjust as needed

# ── Fit power-law per species ─────────────────────────────────────────────────
models <- lidar_biomass_subset %>%
  group_by(Species) %>%
  do(fit = lm(log(biomass_allometry_volume) ~ log(DBH), data = .)) %>%
  mutate(
    intercept = exp(coef(fit)[1]),    # a  (back-transformed)
    exponent  = coef(fit)[2],         # b
    r2        = summary(fit)$r.squared
  )

print(models %>% select(Species, intercept, exponent, r2))

# ── Visualize ─────────────────────────────────────────────────────────────────
windows()
ggplot(lidar_biomass_subset, 
       aes(x = DBH, y = biomass_allometry_volume, color = Species)) +
  geom_point(alpha = 0.6) +
  stat_smooth(method = "lm", formula = y ~ x,
              aes(group = Species), se = TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "DBH–Biomass Relationship by Species",
    subtitle = paste("Platform:", unique(lidar_biomass_subset$Platform)),
    x = "DBH (cm)",
    y = "Biomass — allometry volume (ton)"
  ) +
  theme_minimal()







# ── 1. Combine all plot DBH reconstructions ───────────────────────────────────
all_DBH_recon <- bind_rows(
  MBA_DBH_recon_results %>% mutate(Plot = "MBA"),
  MBB_DBH_recon_results %>% mutate(Plot = "MBB"),
  MBC_DBH_recon_results %>% mutate(Plot = "MBC")
) %>% select(-DBH)

# ── 2. Extract model coefficients ────────────────────────────────────────────
coef_table <- models %>%
  select(Species, intercept, exponent)

# ── 3. Reconstruct biomass ────────────────────────────────────────────────────
biomass_recon <- all_DBH_recon %>%
  pivot_longer(
    cols = starts_with("year_"),
    names_to = "Year",
    values_to = "DBH_recon"
  ) %>%
  mutate(Year = as.integer(gsub("year_", "", Year))) %>%
  left_join(coef_table, by = "Species") %>%
  mutate(
    DBH_recon  = as.numeric(DBH_recon),
    Biomass_t = intercept * DBH_recon ^ exponent
  )

# ── 4. Summarise across all plots combined ────────────────────────────────────
biomass_summary <- biomass_recon %>%
  filter(Year >= 2009 & Year <= 2021) %>%
  filter(!is.na(DBH_recon) & DBH_recon != 0 & !is.na(Biomass_t))%>%
  group_by(Species, Year) %>%
  summarise(
    N_trees       = n(),
    Total_biomass = sum(Biomass_t, na.rm = TRUE),
    Mean_biomass  = mean(Biomass_t, na.rm = TRUE),
    .groups = "drop"
  )

# ── 5. Visualise ──────────────────────────────────────────────────────────────
windows()
ggplot(biomass_summary, aes(x = Year, y = Total_biomass, color = Species)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(
    title    = "Reconstructed Stand Biomass 2008–Present",
    subtitle = "Power-law allometry from TLS–UAV DBH–biomass relationship",
    x        = "Year",
    y        = "Total Biomass (ton)"
  ) +
  theme_minimal()


####BIOMASS RECONSTRUCTION########






























####BIOMASS RECONSTRUCTION########
##LiDAR TREES BIOMASS current date 
#subset the DBH of the biomassI have available for Sap Flow
library(dplyr)
library(readr)

lidar_biomass <- read_csv(
  "\\\\snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/species_metrics_05mvox_033126.csv"
)

dbh_tol <- 100  # cm tolerance

lidar_biomass_subset <- lidar_biomass %>%
  mutate(
    Species = case_when(
      species == "ponderosa"    ~ "PIPO",
      species == "strobiformis" ~ "PISF",
      species == "menziesii"    ~ "PSME",
      TRUE                      ~ species
    )
  ) %>%
  rowwise() %>%
  filter(
    Species %in% names(dbh_list) &&
      any(abs(DBH - dbh_list[[Species]]) <= dbh_tol)
  ) %>%
  ungroup()




library(ggplot2)
library(broom)

# ── Subset by platform ────────────────────────────────────────────────────────
# Change "TLS" to whatever platform value(s) you want to keep
lidar_biomass_subset <- lidar_biomass_subset %>%
  filter(Platform == "TLS_UAV")   # e.g. "TLS", "ALS", "UAS" — adjust as needed

# ── Fit power-law per species ─────────────────────────────────────────────────
models <- lidar_biomass_subset %>%
  group_by(Species) %>%
  do(fit = lm(log(biomass_allometry_volume) ~ log(DBH), data = .)) %>%
  mutate(
    intercept = exp(coef(fit)[1]),    # a  (back-transformed)
    exponent  = coef(fit)[2],         # b
    r2        = summary(fit)$r.squared
  )

print(models %>% select(Species, intercept, exponent, r2))

# ── Visualize ─────────────────────────────────────────────────────────────────
windows()
ggplot(lidar_biomass_subset, 
       aes(x = DBH, y = biomass_allometry_volume, color = Species)) +
  geom_point(alpha = 0.6) +
  stat_smooth(method = "lm", formula = y ~ x,
              aes(group = Species), se = TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "DBH–Biomass Relationship by Species",
    subtitle = paste("Platform:", unique(lidar_biomass_subset$Platform)),
    x = "DBH (cm)",
    y = "Biomass — allometry volume (ton)"
  ) +
  theme_minimal()







# ── 1. Combine all plot DBH reconstructions ───────────────────────────────────
all_DBH_recon <- bind_rows(
  MBA_DBH_recon_results %>% mutate(Plot = "MBA"),
  MBB_DBH_recon_results %>% mutate(Plot = "MBB"),
  MBC_DBH_recon_results %>% mutate(Plot = "MBC")
) %>% select(-DBH)

# ── 2. Extract model coefficients ────────────────────────────────────────────
coef_table <- models %>%
  select(Species, intercept, exponent)

# ── 3. Reconstruct biomass ────────────────────────────────────────────────────
biomass_recon <- all_DBH_recon %>%
  pivot_longer(
    cols = starts_with("year_"),
    names_to = "Year",
    values_to = "DBH_recon"
  ) %>%
  mutate(Year = as.integer(gsub("year_", "", Year))) %>%
  left_join(coef_table, by = "Species") %>%
  mutate(
    DBH_recon  = as.numeric(DBH_recon),
    Biomass_t = intercept * DBH_recon ^ exponent
  )

# ── 4. Summarise across all plots combined ────────────────────────────────────
biomass_summary <- biomass_recon %>%
  filter(Year >= 2009 & Year <= 2021) %>%
  filter(!is.na(DBH_recon) & DBH_recon != 0 & !is.na(Biomass_t))%>%
  group_by(Species, Year) %>%
  summarise(
    N_trees       = n(),
    Total_biomass = sum(Biomass_t, na.rm = TRUE),
    Mean_biomass  = mean(Biomass_t, na.rm = TRUE),
    .groups = "drop"
  )

# ── 5. Visualise ──────────────────────────────────────────────────────────────
windows()
ggplot(biomass_summary, aes(x = Year, y = Total_biomass, color = Species)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(
    title    = "Reconstructed Stand Biomass 2008–Present",
    subtitle = "Power-law allometry from TLS–UAV DBH–biomass relationship",
    x        = "Year",
    y        = "Total Biomass (ton)"
  ) +
  theme_minimal()


####BIOMASS RECONSTRUCTION########
















####BIOMASS RECONSTRUCTION########
##LiDAR TREES BIOMASS current date 
#subset the DBH of the biomassI have available for Sap Flow
library(dplyr)
library(readr)

lidar_biomass <- read_csv(
  "\\\\snre-snow/projects/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/species_metrics_05mvox_033126.csv"
)

dbh_tol <- 100  # cm tolerance

lidar_biomass_subset <- lidar_biomass %>%
  mutate(
    Species = case_when(
      species == "ponderosa"    ~ "PIPO",
      species == "strobiformis" ~ "PISF",
      species == "menziesii"    ~ "PSME",
      TRUE                      ~ species
    )
  ) %>%
  rowwise() %>%
  filter(
    Species %in% names(dbh_list) &&
      any(abs(DBH - dbh_list[[Species]]) <= dbh_tol)
  ) %>%
  ungroup()




library(ggplot2)
library(broom)

# ── Subset by platform ────────────────────────────────────────────────────────
# Change "TLS" to whatever platform value(s) you want to keep
lidar_biomass_subset <- lidar_biomass_subset %>%
  filter(Platform == "TLS_UAV")   # e.g. "TLS", "ALS", "UAS" — adjust as needed

# ── Fit power-law per species ─────────────────────────────────────────────────
models <- lidar_biomass_subset %>%
  group_by(Species) %>%
  do(fit = lm(log(biomass_allometry_volume) ~ log(DBH), data = .)) %>%
  mutate(
    intercept = exp(coef(fit)[1]),    # a  (back-transformed)
    exponent  = coef(fit)[2],         # b
    r2        = summary(fit)$r.squared
  )

print(models %>% select(Species, intercept, exponent, r2))

# ── Visualize ─────────────────────────────────────────────────────────────────
windows()
ggplot(lidar_biomass_subset, 
       aes(x = DBH, y = biomass_allometry_volume, color = Species)) +
  geom_point(alpha = 0.6) +
  stat_smooth(method = "lm", formula = y ~ x,
              aes(group = Species), se = TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "DBH–Biomass Relationship by Species",
    subtitle = paste("Platform:", unique(lidar_biomass_subset$Platform)),
    x = "DBH (cm)",
    y = "Biomass — allometry volume (ton)"
  ) +
  theme_minimal()







# ── 1. Combine all plot DBH reconstructions ───────────────────────────────────
all_DBH_recon <- bind_rows(
  MBA_DBH_recon_results %>% mutate(Plot = "MBA"),
  MBB_DBH_recon_results %>% mutate(Plot = "MBB"),
  MBC_DBH_recon_results %>% mutate(Plot = "MBC")
) %>% select(-DBH)

# ── 2. Extract model coefficients ────────────────────────────────────────────
coef_table <- models %>%
  select(Species, intercept, exponent)

# ── 3. Reconstruct biomass ────────────────────────────────────────────────────
biomass_recon <- all_DBH_recon %>%
  pivot_longer(
    cols = starts_with("year_"),
    names_to = "Year",
    values_to = "DBH_recon"
  ) %>%
  mutate(Year = as.integer(gsub("year_", "", Year))) %>%
  left_join(coef_table, by = "Species") %>%
  mutate(
    DBH_recon  = as.numeric(DBH_recon),
    Biomass_t = intercept * DBH_recon ^ exponent
  )

# ── 4. Summarise across all plots combined ────────────────────────────────────
biomass_summary <- biomass_recon %>%
  filter(Year >= 2009 & Year <= 2021) %>%
  filter(!is.na(DBH_recon) & DBH_recon != 0 & !is.na(Biomass_t))%>%
  group_by(Species, Year) %>%
  summarise(
    N_trees       = n(),
    Total_biomass = sum(Biomass_t, na.rm = TRUE),
    Mean_biomass  = mean(Biomass_t, na.rm = TRUE),
    .groups = "drop"
  )

# ── 5. Visualise ──────────────────────────────────────────────────────────────
windows()
ggplot(biomass_summary, aes(x = Year, y = Total_biomass, color = Species)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(
    title    = "Reconstructed Stand Biomass 2008–Present",
    subtitle = "Power-law allometry from TLS–UAV DBH–biomass relationship",
    x        = "Year",
    y        = "Total Biomass (ton)"
  ) +
  theme_minimal()


####BIOMASS RECONSTRUCTION########




# Volume ~ Crown Area: simple log-log fit for upscaling
windows()
ggplot(lidar_biomass_subset,
       aes(x = CrownArea_m2, y = VoxelVolume_m3, color = Species)) +
  geom_point(alpha = 0.6) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE, aes(group = Species)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Volume ~ Crown Area (TLS-UAV)",
    x     = "Crown area (m²)",
    y     = "Biomass allometry volume (ton)"
  ) +
  theme_minimal()

# Coefficients per species
crown_models <- lidar_biomass_subset %>%
  filter(!is.na(CrownArea_m2) & CrownArea_m2 > 0) %>%
  group_by(Species) %>%
  do(fit = lm(log(VoxelVolume_m3) ~ log(CrownArea_m2), data = .)) %>%
  mutate(
    intercept = exp(coef(fit)[1]),
    exponent  = coef(fit)[2],
    r2        = summary(fit)$r.squared,
    p_value   = summary(fit)$coefficients[2, 4]  # p-value for the slope
  )

print(crown_models %>% select(Species, intercept, exponent, r2, p_value))





crown_models <- lidar_biomass_subset %>%
  filter(!is.na(CrownArea_m2) & CrownArea_m2 > 0) %>%
  summarise(fit = list(lm(log(VoxelVolume_m3) ~ log(CrownArea_m2), data = cur_data()))) %>%
  mutate(
    intercept = exp(coef(fit[[1]])[1]),
    exponent  = coef(fit[[1]])[2],
    r2        = summary(fit[[1]])$r.squared,
    p_value   = summary(fit[[1]])$coefficients[2, 4]
  )

print(crown_models %>% select(intercept, exponent, r2, p_value))






# ── All pairwise relationships against VoxelVolume_m3 ─────────────────────────
library(patchwork)

p1 <- ggplot(lidar_biomass_subset, aes(x = DBH, y = VoxelVolume_m3, color = Species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE, aes(group = Species)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "DBH (cm)", y = "Voxel volume (m³)") +
  theme_minimal() + theme(legend.position = "none")

p2 <- ggplot(lidar_biomass_subset, aes(x = Height, y = VoxelVolume_m3, color = Species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE, aes(group = Species)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Height (m)", y = "Voxel volume (m³)") +
  theme_minimal() + theme(legend.position = "none")

p3 <- ggplot(lidar_biomass_subset, aes(x = CrownArea_m2, y = VoxelVolume_m3, color = Species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE, aes(group = Species)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Crown area (m²)", y = "Voxel volume (m³)") +
  theme_minimal() + theme(legend.position = "none")

p4 <- ggplot(lidar_biomass_subset, aes(x = CrownArea_m2, y = Height, color = Species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE, aes(group = Species)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Crown area (m²)", y = "Height (m)") +
  theme_minimal() + theme(legend.position = "none")

p5 <- ggplot(lidar_biomass_subset, aes(x = VoxelVolume_m3, y = biomass_allometry_volume, color = Species)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ x, se = TRUE, aes(group = Species)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Voxel volume (m³)", y = "Biomass allometry volume (ton)") +
  theme_minimal() + theme(legend.position = "none")

# shared legend from one plot
legend <- ggplot(lidar_biomass_subset, aes(x = DBH, y = VoxelVolume_m3, color = Species)) +
  geom_point() + theme_minimal() +
  theme(legend.position = "bottom")
shared_legend <- cowplot::get_legend(legend)

windows()
cowplot::plot_grid(
  cowplot::plot_grid(p1, p2, p3, p4, p5, ncol = 3),
  shared_legend,
  ncol = 1, rel_heights = c(1, 0.05)
)






fit_stats <- function(data, x_var, y_var, label) {
  data %>%
    filter(!is.na(.data[[x_var]]) & .data[[x_var]] > 0 &
             !is.na(.data[[y_var]]) & .data[[y_var]] > 0) %>%
    group_by(Species) %>%
    summarise(
      fit = list(lm(log(.data[[y_var]]) ~ log(.data[[x_var]]), data = cur_data())),
      n   = n(),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      relationship = label,
      intercept    = exp(coef(fit)[1]),
      exponent     = coef(fit)[2],
      r2           = summary(fit)$r.squared,
      p_value      = summary(fit)$coefficients[2, 4]
    ) %>%
    ungroup() %>%
    select(relationship, Species, intercept, exponent, r2, p_value, n)
}

stats_all <- bind_rows(
  fit_stats(lidar_biomass_subset, "DBH",            "VoxelVolume_m3",           "DBH ~ VoxelVol"),
  fit_stats(lidar_biomass_subset, "Height",          "VoxelVolume_m3",           "Height ~ VoxelVol"),
  fit_stats(lidar_biomass_subset, "CrownArea_m2",    "VoxelVolume_m3",           "CrownArea ~ VoxelVol"),
  fit_stats(lidar_biomass_subset, "CrownArea_m2",    "Height",                   "CrownArea ~ Height"),
  fit_stats(lidar_biomass_subset, "VoxelVolume_m3",  "biomass_allometry_volume", "VoxelVol ~ BiomassAllom"),
  
  # Wood volume checks
  fit_stats(lidar_biomass_subset, "DBH",             "wood_vol_m3",            "DBH ~ WoodVol"),
  fit_stats(lidar_biomass_subset, "Height",           "wood_vol_m3",            "Height ~ WoodVol"),
  fit_stats(lidar_biomass_subset, "CrownArea_m2",     "wood_vol_m3",            "CrownArea ~ WoodVol"),
  fit_stats(lidar_biomass_subset, "VoxelVolume_m3",   "wood_vol_m3",            "VoxelVol ~ WoodVol"),
  fit_stats(lidar_biomass_subset, "wood_vol_m3",    "biomass_allometry_volume", "WoodVol ~ BiomassAllom")
)

print(stats_all, n = Inf)
















































load_tree <- function(platform, site, parent_dir, tree) {
  folder <- file.path(parent_dir, paste0(platform, "_", site))
  f <- list.files(folder, pattern = paste0("^", site, "_.*tree", tree, "\\.las$"),
                  full.names = TRUE)
  if (length(f) == 0) { message("Not found: ", platform, "_", site, " tree", tree); return(NULL) }
  las <- readLAS(f[1])
  if (is.empty(las)) return(NULL)
  data.frame(X = las@data$X, Y = las@data$Y, Z = las@data$Z, Platform = platform)
}

# ── Load all 4 platforms ──────────────────────────────────────────────────────
target_site <- "MBA"
target_tree <- "6"
plot_name   <- paste(target_site, "Tree", target_tree)

tls           <- load_tree("TLS",           target_site, parent_dir, target_tree)
nadir         <- load_tree("nadir",         target_site, parent_dir, target_tree)
nadir_oblique <- load_tree("nadir_oblique", target_site, parent_dir, target_tree)
tls_uav       <- load_tree("TLS_UAV",       target_site, parent_dir, target_tree)

# ── Common bounding box across all ───────────────────────────────────────────
all_pts <- rbind(tls, nadir, nadir_oblique, tls_uav)
xlim <- range(all_pts$X, na.rm = TRUE)
ylim <- range(all_pts$Y, na.rm = TRUE)
zlim <- range(all_pts$Z, na.rm = TRUE)

# ── Colors per platform ───────────────────────────────────────────────────────
col_tls           <- "#4CAF50"
col_nadir         <- "#2196F3"
col_nadir_oblique <- "#FF9800"
col_tls_uav       <- "red"

# ── Helper: overlay two clouds in current panel ───────────────────────────────
overlay_plot <- function(df1, col1, lab1, df2, col2, lab2, title) {
  rgl::plot3d(df1$X, df1$Y, df1$Z,
              col = col1, size = 2, alpha = 0.6,
              xlim = xlim, ylim = ylim, zlim = zlim,
              main = title, xlab = "X", ylab = "Y", zlab = "Z")
  rgl::points3d(df2$X, df2$Y, df2$Z,
                col = col2, size = 2, alpha = 0.6)
  rgl::grid3d(c("x", "y", "z"))
  rgl::legend3d("topright", legend = c(lab1, lab2),
                col = c(col1, col2), pch = 16, cex = 1.2)
}

# ── 3 panels: TLS (green) overlaid with each platform ────────────────────────
rgl::open3d()
rgl::par3d(windowRect = c(50, 50, 1800, 700))
rgl::mfrow3d(1, 3)

overlay_plot(tls, col_tls, "TLS",
             nadir, col_nadir, "Nadir",
             paste0(plot_name, " | TLS vs Nadir"))

rgl::next3d()
overlay_plot(tls, col_tls, "TLS",
             nadir_oblique, col_nadir_oblique, "Nadir Oblique",
             paste0(plot_name, " | TLS vs Nadir Oblique"))

rgl::next3d()
overlay_plot(tls, col_tls, "TLS",
             tls_uav, col_tls_uav, "TLS-UAV",
             paste0(plot_name, " | TLS vs TLS-UAV"))







target_site <- "MBA"
plot_name   <- target_site

# ── Load entire plot (all trees combined) ─────────────────────────────────────
load_plot <- function(platform, site, parent_dir) {
  folder <- file.path(parent_dir, paste0(platform, "_", site))
  files  <- list.files(folder, pattern = "\\.las$", full.names = TRUE)
  if (length(files) == 0) { message("Not found: ", platform, "_", site); return(NULL) }
  pts <- lapply(files, function(f) {
    las <- readLAS(f)
    if (is.empty(las)) return(NULL)
    data.frame(X = las@data$X, Y = las@data$Y, Z = las@data$Z)
  })
  df <- do.call(rbind, Filter(Negate(is.null), pts))
  df$Platform <- platform
  df
}

tls           <- load_plot("TLS",           target_site, parent_dir)
nadir         <- load_plot("nadir",         target_site, parent_dir)
nadir_oblique <- load_plot("nadir_oblique", target_site, parent_dir)
tls_uav       <- load_plot("TLS_UAV",       target_site, parent_dir)

# ── Common bounding box ───────────────────────────────────────────────────────
all_pts <- rbind(tls, nadir, nadir_oblique, tls_uav)
xlim <- range(all_pts$X, na.rm = TRUE)
ylim <- range(all_pts$Y, na.rm = TRUE)
zlim <- range(all_pts$Z, na.rm = TRUE)

# ── Colors ────────────────────────────────────────────────────────────────────
col_tls           <- "#4CAF50"
col_nadir         <- "#2196F3"
col_nadir_oblique <- "#FF9800"
col_tls_uav       <- "red"

# ── Helper ────────────────────────────────────────────────────────────────────
overlay_plot <- function(df1, col1, lab1, df2, col2, lab2, title) {
  rgl::plot3d(df1$X, df1$Y, df1$Z,
              col = col1, size = 1, alpha = 0.5,
              xlim = xlim, ylim = ylim, zlim = zlim,
              main = title, xlab = "X", ylab = "Y", zlab = "Z")
  rgl::points3d(df2$X, df2$Y, df2$Z,
                col = col2, size = 1, alpha = 0.5)
  rgl::grid3d(c("x", "y", "z"))
  rgl::legend3d("topright", legend = c(lab1, lab2),
                col = c(col1, col2), pch = 16, cex = 1.2)
}

# ── 3 panels ──────────────────────────────────────────────────────────────────
rgl::open3d()
rgl::par3d(windowRect = c(50, 50, 1800, 700))
rgl::mfrow3d(1, 3)

overlay_plot(tls, col_tls, "TLS",
             nadir, col_nadir, "Nadir",
             paste0(plot_name, " | TLS vs Nadir"))
rgl::next3d()
overlay_plot(tls, col_tls, "TLS",
             nadir_oblique, col_nadir_oblique, "Nadir Oblique",
             paste0(plot_name, " | TLS vs Nadir Oblique"))
rgl::next3d()
overlay_plot(tls, col_tls, "TLS",
             tls_uav, col_tls_uav, "TLS-UAV",
             paste0(plot_name, " | TLS vs TLS-UAV"))















library(lidR)
library(rgl)
library(data.table)

# ── PATHS ─────────────────────────────────────────────────────────────────────
leaf_folder  <- "Z:/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/CH2_LiDAR/classified_split/leaf"
woody_folder <- "Z:/Babst_Lidar_treering_CLN/CLNorton/TLiDAR/data/CH2_LiDAR/classified_split/wood"

# ── READ + PARSE PLOT FROM FILENAME ───────────────────────────────────────────
read_folder <- function(folder, class_label) {
  files <- list.files(folder, pattern = "\\.las$|\\.laz$", full.names = TRUE)
  rbindlist(lapply(files, function(f) {
    las <- readLAS(f, select = "xyz")
    if (is.empty(las)) return(NULL)
    dt       <- as.data.table(las@data)
    dt[, class := class_label]
    dt[, plot  := stringr::str_extract(basename(f), "MB[ABC]")]
    dt
  }))
}

leaf_pts  <- read_folder(leaf_folder,  "leaf")
woody_pts <- read_folder(woody_folder, "woody")
all_pts   <- rbindlist(list(leaf_pts, woody_pts))

# ── THIN PER PLOT ─────────────────────────────────────────────────────────────
max_pts <- 300000
all_pts <- all_pts[, .SD[sample(.N, min(.N, max_pts))], by = plot]

# ── COLORS ────────────────────────────────────────────────────────────────────
all_pts[, color := ifelse(class == "leaf", "#2E8B57", "#8B4513")]

plots <- c("MBA", "MBB", "MBC")

# ── PLOT: 1 row x 3 panels ────────────────────────────────────────────────────
par3d(windowRect = c(0, 0, 2400, 800))
mfrow3d(1, 3, sharedMouse = TRUE)

for (plt in plots) {
  sub <- all_pts[plot == plt]
  next3d()
  bg3d("white")
  points3d(sub$X, sub$Y, sub$Z, col = sub$color, size = 3)
  bbox3d(col = "black", emission = "black", specular = "black", alpha = 0.3,
         xlen = 5, ylen = 5, zlen = 5)
  title3d(main = plt, xlab = "X", ylab = "Y", zlab = "Z", col = "black")
}

legend3d("topright",
         legend = c("Leaf", "Woody"),
         col    = c("#2E8B57", "#8B4513"),
         pch    = 16, cex = 1.2)

