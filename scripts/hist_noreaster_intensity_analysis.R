# ==============================================================================

# Analysis of historical nor'easter intensity characteristics in CFSR, ERA5, 
#   CRCM5, and CRCM5-CMIP6

# ==============================================================================

# Associated publication: Weaver, M.M. et al. (2026). Historical nor'easter
#                         intensity trends in reanalysis and regional climate
#                         simulations: A cross-dataset perspective.

# This file contains all code necessary to recreate the main text figures, 
#   tables, and supplementary figures in the aforementioned publication.

# ==============================================================================

# Load necessary .R script(s): contains functions used throughout the analysis
source("/Path/to/functions/file/historical_noreaster_analysis_functions.R")

# Load necessary R packages:
library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library("RNetCDF")
library("MASS")
library("fields")
library("maps")

# Define base directory
base_dir <- "/Path/to/all/data"

# ==============================================================================

# Load CFSR data, create CFSR lat/lon grid, apply post-processing filters

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Load CFSR data

# Initialize storage list for CFSR nor'easters
cfsr_noreasters <- list()

# File paths
cfsr_path <- file.path("/Path/to/CFSR/data")
csv_files <- list.files(cfsr_path, pattern = "\\.csv$", full.names = TRUE)

# Load cfsr data
for (csv_file in csv_files) {
  year <- extract_year(csv_file)
  
  if (!is.na(year)) {
    tryCatch({
      df <- read_storm_csv(csv_file)
      cfsr_noreasters[[as.character(year)]] <- df
      cat(paste("  ", year, ":", nrow(df), "storms\n"))
    }, error = function(e) {
      cat(paste("  Error reading", csv_file, ":", e$message, "\n"))
    })
  } else {
    cat(paste("  Warning: Could not extract year from", basename(csv_file), "\n"))
  }
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Recreate the cfsr grid
cfsr_grid <- create_cfsr_grid()
cfsr_lon <- cfsr_grid$lon
cfsr_lat <- cfsr_grid$lat

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Apply post-processing filters to CFSR nor'easters: 24 hr lifetime, 1000 km 
#   travel distance, remove tropical cyclones, U.S. east coast proximity, 
#   intensity subsets, common time period(s), seasons

# Define wind_1979_2014 data directory (needed for tropical cyclone filtering)
cfsr_wind_dir <- file.path(base_dir, "path/to/CFSR/data/CFSR_wind_data")

cfsr_noreasters <- filter_cfsr_noreasters(cfsr_noreasters, cfsr_lat, cfsr_lon, 
                                          wind_dirs = cfsr_wind_dir,
                                          apply_region_filter = FALSE, 
                                          apply_TC_filter = TRUE,
                                          TC_min_lat = 30, TC_min_slp = 98000, 
                                          TC_max_wind = 20)

# Filter for just strong nor'easters (lifetime min. sea-level pressure ≤ 98000 Pa)
cfsr_strong_noreasters <- apply_SLP_filter(cfsr_noreasters, limit = 98000, 
                                           lim_min = FALSE, era5 = TRUE, 
                                           dataframe = TRUE)
# Filter for just weak nor'easters (lifetime min. sea-level pressure > 98000 Pa)
cfsr_weak_noreasters <- apply_SLP_filter(cfsr_noreasters, limit = 98000, 
                                         lim_min = TRUE, era5 = TRUE, 
                                         dataframe = TRUE)
                                         # era5 = TRUE indicates non-nested data structure

# Filter CFSR storms to common time period w/ other datasets
cfsr_strong_noreasters_1979_2014 <- filter_years_era5(cfsr_strong_noreasters,   # CFSR uses the same function as ERA5 for this
                                                      start_year = 1980,  # CFSR nor'easter seasons are
                                                      end_year = 2014)    # keyed to the year they end

# Filter CFSR strong storms by season:

# Strong dataset
cfsr_strong_noreasters_SON <- filter_season_era5(cfsr_strong_noreasters, "SON")
cfsr_strong_noreasters_DJF <- filter_season_era5(cfsr_strong_noreasters, "DJF")
cfsr_strong_noreasters_MA  <- filter_season_era5(cfsr_strong_noreasters, "MA")

# 1979-2014
cfsr_strong_noreasters_1979_2014_SON <- filter_season_era5(cfsr_strong_noreasters_1979_2014, "SON")
cfsr_strong_noreasters_1979_2014_DJF <- filter_season_era5(cfsr_strong_noreasters_1979_2014, "DJF")
cfsr_strong_noreasters_1979_2014_MA  <- filter_season_era5(cfsr_strong_noreasters_1979_2014, "MA")

# ==============================================================================

# Load ERA5 data, create ERA5 lat/lon grid, apply post-processing filters

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Load ERA5 data

# Initialize storage list
era5_noreasters <- list()

# File paths
era5_path <- file.path("/Path/to/ERA5/data")
csv_files <- list.files(era5_path, pattern = "\\.csv$", full.names = TRUE)

# Load ERA5 data
for (csv_file in csv_files) {
  year <- extract_year(csv_file)
  
  if (!is.na(year)) {
    tryCatch({
      df <- read_storm_csv(csv_file)
      era5_noreasters[[as.character(year)]] <- df
      cat(paste("  ", year, ":", nrow(df), "storms\n"))
    }, error = function(e) {
      cat(paste("  Error reading", csv_file, ":", e$message, "\n"))
    })
  } else {
    cat(paste("  Warning: Could not extract year from", basename(csv_file), "\n"))
  }
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Recreate the ERA5 grid

era5_grid <- create_era5_grid()
era5_lon <- era5_grid$lon
era5_lat <- era5_grid$lat

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Apply post-processing to ERA5 nor'easters: tropical cyclone removal, intensity
#   subsets, and common time period(s)
#   ***since ERA5 nor'easters are from an existing tracked cyclone dataset
#      (Chen et al., 2025), the 24 hr lifetime, 1000 km travel distance, and
#      U.S. east coast post-processing filters have already been applied***

# Filter for just strong nor'easters (lifetime min. sea-level pressure ≤ 98000 Pa)
era5_strong_noreasters <- apply_SLP_filter(era5_noreasters, limit = 98000, 
                                           lim_min = FALSE, era5 = TRUE, 
                                           dataframe = TRUE)
# Filter for just weak nor'easters (lifetime min. sea-level pressure > 98000 Pa)
era5_weak_noreasters <- apply_SLP_filter(era5_noreasters, limit = 98000, 
                                         lim_min = TRUE, era5 = TRUE, 
                                         dataframe = TRUE)

# Apply tropical cyclone filter to strong ERA5 nor'easters (no era5_wind_dir directory 
#   needed, ERA5 wind data is integrated into storm track data (Chen et al., 2025))
era5_strong_noreasters <- TC_filter_ERA5(era5_strong_noreasters, 
                                         era5_lon, era5_lat, 
                                         min_lat = 30, min_slp = 98000, 
                                         max_wind = 20)
# Apply tropical cyclone filter to weak ERA5 nor'easters
era5_weak_noreasters <- TC_filter_ERA5(era5_weak_noreasters, era5_lon, era5_lat, 
                                       min_lat = 30, min_slp = 98000, 
                                       max_wind = 20)
# Apply tropical cyclone filter to full ERA5 nor'easter dataset
era5_noreasters <- TC_filter_ERA5(era5_noreasters, era5_lon, era5_lat, 
                                  min_lat = 30, min_slp = 98000, 
                                  max_wind = 20)

# Filter ERA5 strong nor'easters to common time period with CRCM5-CMIP6 (1950-2014)
era5_strong_noreasters_1950_2014 <- filter_years_era5(era5_strong_noreasters, 
                                                      start_year = 1950,   # ERA5 nor'easter seasons are
                                                      end_year = 2013)     # keyed to the year they begin
# Filter ERA5 strong nor'easters to common time period with CFSR and CRCM5-ERA5 (1979-2014)
era5_strong_noreasters_1979_2014 <- filter_years_era5(era5_strong_noreasters, 
                                                      start_year = 1979, 
                                                      end_year = 2013)

# Filter ERA5 strong storms by season:

# Strong dataset
era5_strong_noreasters_SON <- filter_season_era5(era5_strong_noreasters, "SON")
era5_strong_noreasters_DJF <- filter_season_era5(era5_strong_noreasters, "DJF")
era5_strong_noreasters_MA  <- filter_season_era5(era5_strong_noreasters, "MA")

# 1950-2014
era5_strong_noreasters_1950_2014_SON <- filter_season_era5(era5_strong_noreasters_1950_2014, "SON")
era5_strong_noreasters_1950_2014_DJF <- filter_season_era5(era5_strong_noreasters_1950_2014, "DJF")
era5_strong_noreasters_1950_2014_MA  <- filter_season_era5(era5_strong_noreasters_1950_2014, "MA")

# 1979-2014
era5_strong_noreasters_1979_2014_SON <- filter_season_era5(era5_strong_noreasters_1979_2014, "SON")
era5_strong_noreasters_1979_2014_DJF <- filter_season_era5(era5_strong_noreasters_1979_2014, "DJF")
era5_strong_noreasters_1979_2014_MA  <- filter_season_era5(era5_strong_noreasters_1979_2014, "MA")

# ==============================================================================

# Load CRCM5 data (CRCM5-ERA5 & CRCM5-CMIP6), create CRCM5 lat/lon grid, apply 
#   post-processing filters

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Load CRCM5 data

# Define pilot GCMs (non-ERA5 CRCM5 simulations)
gcm_models <- c("canesm5", "cnrm", "noresm2", "mpi")

# Define model paths
models <- list(
  canesm5 = file.path(base_dir, "canesm5_historical_unfiltered_tracks"),
  cnrm = file.path(base_dir, "cnrm_historical_unfiltered_tracks"),
  noresm2 = file.path(base_dir, "noresm2_historical_unfiltered_tracks"),
  mpi = file.path(base_dir, "mpi_historical_unfiltered_tracks/sim1"),
  era5 = file.path(base_dir, "crcm5-era5_unfiltered_tracks/v1-r1")
)

# Initialize storage list
crcm5_noreasters <- list()

for (model in names(models)) {
  cat(paste("Loading", model, "...\n"))
  
  model_path <- models[[model]]
  csv_files <- list.files(model_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Create a list to store this model's data
  model_data <- list()
  
  for (csv_file in csv_files) {
    year <- extract_year(csv_file)
    
    if (!is.na(year)) {
      tryCatch({
        df <- read_storm_csv(csv_file)
        model_data[[as.character(year)]] <- df
        cat(paste("  ", year, ":", nrow(df), "storms\n"))
      }, error = function(e) {
        cat(paste("  Error reading", csv_file, ":", e$message, "\n"))
      })
    } else {
      cat(paste("  Warning: Could not extract year from", basename(csv_file), "\n"))
    }
  }
  
  crcm5_noreasters[[model]] <- model_data
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Create lat/lon grids based on CRCM5 spatial extent/resolution (0.11 degrees)

# Load data file for lon/lat arrays: unlike CFSR & ERA5, this requires a data file 
#   (any CRCM5 file will do). Available via Thredds (see manuscript Data Availability)
psl.data = open.nc("/Path/to/raw/CRCM5/data/file/psl_NAM-12_CanESM5_historical_r1i1p2f1_OURANOS_CRCM5_v1-r1_3hr_195001010000-195012312100.nc")

# Pull out regular lon and lat variables from psl.data:
lon = var.get.nc(psl.data, 'lon')
lat = var.get.nc(psl.data, 'lat')

# Transpose to fix R's interpretation
crcm5_lat <- t(lat)
crcm5_lon <- t(lon)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Apply post-processing filters to CRCM5 nor'easters: 24 hr lifetime, 1000 km 
#   travel distance, remove tropical cyclones, U.S. east coast proximity, 
#   intensity subsets, common time period(s), seasons

# Define wind_1979_2014 data directories
crcm5_wind_dirs <- list(
  canesm5 = file.path(base_dir, 
                      "path/to/wind/data/CanESM5_wind_data/canesm5_historical"),
  cnrm = file.path(base_dir, 
                   "path/to/wind/data/CNRM_wind_data/cnrm_historical"),
  noresm2 = file.path(base_dir, 
                      "path/to/wind/data/NorESM2_wind_data/noresm2_historical"),
  mpi = file.path(base_dir, 
                  "path/to/wind/data/MPI_wind_data/mpi_historical1"),
  era5 = file.path(base_dir, 
                   "path/to/wind/data/CRCM5-ERA5_wind_data/era5_v1-r1")
)

# Apply main post-processing function (USEC, time, distance, hurricanes, 
#   genesis, NO target region)
crcm5_noreasters <- filter_crcm5_noreasters(crcm5_noreasters, 
                                            crcm5_lat, crcm5_lon, 
                                            wind_dirs = crcm5_wind_dirs,
                                            apply_region_filter = FALSE, 
                                            apply_TC_filter = TRUE,
                                            TC_min_lat = 30, TC_min_slp = 98000, 
                                            TC_max_wind = 20)

# Filter for just strong nor'easters (lifetime min. sea-level pressure ≤ 98000 Pa)
strong_crcm5_noreasters <- apply_SLP_filter(crcm5_noreasters, limit = 98000, 
                                            lim_min = FALSE, era5 = FALSE, 
                                            dataframe = TRUE)
# Filter for just weak nor'easters (lifetime min. sea-level pressure > 98000 Pa)
weak_crcm5_noreasters <- apply_SLP_filter(crcm5_noreasters, limit = 98000, 
                                          lim_min = TRUE, era5 = FALSE, 
                                          dataframe = TRUE)

# Filter CRCM5 storms to common time period w/ other datasets
strong_crcm5_noreasters_1979_2014 <- filter_years_multi_model(strong_crcm5_noreasters, 
                                                        start_year = 1980,  # CRCM5 nor'easter seasons are
                                                        end_year = 2014)    # keyed to the year they end

# Filter CRCM5 STRONG dataset by season

# 1950-2014
strong_crcm5_noreasters_SON <- filter_season_multi_model(strong_crcm5_noreasters, "SON")
strong_crcm5_noreasters_DJF <- filter_season_multi_model(strong_crcm5_noreasters, "DJF")
strong_crcm5_noreasters_MA  <- filter_season_multi_model(strong_crcm5_noreasters, "MA")

# 1988-89 to 2013-14
strong_crcm5_noreasters_1979_2014_SON <- filter_season_multi_model(strong_crcm5_noreasters_1979_2014, "SON")
strong_crcm5_noreasters_1979_2014_DJF <- filter_season_multi_model(strong_crcm5_noreasters_1979_2014, "DJF")
strong_crcm5_noreasters_1979_2014_MA  <- filter_season_multi_model(strong_crcm5_noreasters_1979_2014, "MA")

# ==============================================================================

# Calculate nor'easter frequencies for each dataset (full dataset, strong & weak
#   subsets, & strong nor'easters for each season)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# CFSR frequencies

# Define years for time series
cfsr_years <- (as.numeric(names(cfsr_noreasters))) + 1

# All
cfsr_noreaster_counts <- sapply(cfsr_noreasters, nrow)
cfsr_noreaster_freq <- data.frame(year = cfsr_years, count = cfsr_noreaster_counts)

# Strong
cfsr_strong_noreaster_counts <- sapply(cfsr_strong_noreasters, nrow)
cfsr_strong_noreaster_freq <- data.frame(year = cfsr_years, 
                                         count = cfsr_strong_noreaster_counts)

# Weak
cfsr_weak_noreaster_counts <- sapply(cfsr_weak_noreasters, nrow)
cfsr_weak_noreaster_freq <- data.frame(year = cfsr_years, 
                                       count = cfsr_weak_noreaster_counts)

# Strong SON
cfsr_strong_noreaster_SON_counts <- sapply(cfsr_strong_noreasters_SON, nrow)
cfsr_strong_noreaster_SON_freq <- data.frame(year = cfsr_years, 
                                             count = cfsr_strong_noreaster_SON_counts)

# Strong DJF
cfsr_strong_noreaster_DJF_counts <- sapply(cfsr_strong_noreasters_DJF, nrow)
cfsr_strong_noreaster_DJF_freq <- data.frame(year = cfsr_years, 
                                             count = cfsr_strong_noreaster_DJF_counts)

# Strong MA
cfsr_strong_noreaster_MA_counts <- sapply(cfsr_strong_noreasters_MA, nrow)
cfsr_strong_noreaster_MA_freq <- data.frame(year = cfsr_years, 
                                            count = cfsr_strong_noreaster_MA_counts)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# ERA5 frequencies

# Define years for time series
era5_years <- (as.numeric(names(era5_noreasters))) + 1

# All
era5_noreaster_counts <- sapply(era5_noreasters, nrow)
era5_noreaster_freq <- data.frame(year = era5_years, 
                                  count = era5_noreaster_counts)

# Strong
era5_strong_noreaster_counts <- sapply(era5_strong_noreasters, nrow)
era5_strong_noreaster_freq <- data.frame(year = era5_years, 
                                         count = era5_strong_noreaster_counts)

# Weak
era5_weak_noreaster_counts <- sapply(era5_weak_noreasters, nrow)
era5_weak_noreaster_freq <- data.frame(year = era5_years, 
                                       count = era5_weak_noreaster_counts)

# Strong SON
era5_strong_noreaster_SON_counts <- sapply(era5_strong_noreasters_SON, nrow)
era5_strong_noreaster_SON_freq <- data.frame(year = era5_years, 
                                             count = era5_strong_noreaster_SON_counts)

# Strong DJF
era5_strong_noreaster_DJF_counts <- sapply(era5_strong_noreasters_DJF, nrow)
era5_strong_noreaster_DJF_freq <- data.frame(year = era5_years, 
                                             count = era5_strong_noreaster_DJF_counts)

# Strong MA
era5_strong_noreaster_MA_counts <- sapply(era5_strong_noreasters_MA, nrow)
era5_strong_noreaster_MA_freq <- data.frame(year = era5_years, 
                                            count = era5_strong_noreaster_MA_counts)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# CRCM5-ERA5 & CRCM5-CMIP6 frequencies

# Full dataset
crcm5_noreaster_frequencies <- list()
for (model in names(crcm5_noreasters)) {
  years <- as.numeric(names(crcm5_noreasters[[model]]))
  counts <- sapply(crcm5_noreasters[[model]], nrow)
  crcm5_noreaster_frequencies[[model]] <- data.frame(year = years, 
                                                     count = counts)
}

# Weak
weak_crcm5_noreaster_frequencies <- list()
for (model in names(weak_crcm5_noreasters)) {
  years <- as.numeric(names(weak_crcm5_noreasters[[model]]))
  counts <- sapply(weak_crcm5_noreasters[[model]], nrow)
  weak_crcm5_noreaster_frequencies[[model]] <- data.frame(year = years, 
                                                          count = counts)
}

# Strong
strong_crcm5_noreaster_frequencies <- list()
for (model in names(strong_crcm5_noreasters)) {
  years <- as.numeric(names(strong_crcm5_noreasters[[model]]))
  counts <- sapply(strong_crcm5_noreasters[[model]], nrow)
  strong_crcm5_noreaster_frequencies[[model]] <- data.frame(year = years, 
                                                            count = counts)
}

# SON
strong_crcm5_noreaster_SON_freq <- list()
for (model in names(strong_crcm5_noreasters_SON)) {
  years <- as.numeric(names(strong_crcm5_noreasters_SON[[model]]))
  counts <- sapply(strong_crcm5_noreasters_SON[[model]], nrow)
  strong_crcm5_noreaster_SON_freq[[model]] <- data.frame(year = years, 
                                                         count = counts)
}
# Print summary
cat("\n=== SUMMARY ===\n")
for (model in names(strong_crcm5_noreaster_SON_freq)) {
  freq_data <- strong_crcm5_noreaster_SON_freq[[model]]
  cat(paste(toupper(model), ":", nrow(freq_data), "seasons,",
            "mean =", round(mean(freq_data$count), 1),
            "std =", round(sd(freq_data$count), 1), "\n"))
}

# DJF
strong_crcm5_noreaster_DJF_freq <- list()
for (model in names(strong_crcm5_noreasters_DJF)) {
  years <- as.numeric(names(strong_crcm5_noreasters_DJF[[model]]))
  counts <- sapply(strong_crcm5_noreasters_DJF[[model]], nrow)
  strong_crcm5_noreaster_DJF_freq[[model]] <- data.frame(year = years, 
                                                         count = counts)
}
# Print summary
cat("\n=== SUMMARY ===\n")
for (model in names(strong_crcm5_noreaster_DJF_freq)) {
  freq_data <- strong_crcm5_noreaster_DJF_freq[[model]]
  cat(paste(toupper(model), ":", nrow(freq_data), "seasons,",
            "mean =", round(mean(freq_data$count), 1),
            "std =", round(sd(freq_data$count), 1), "\n"))
}

# MA
strong_crcm5_noreaster_MA_freq <- list()
for (model in names(strong_crcm5_noreasters_MA)) {
  years <- as.numeric(names(strong_crcm5_noreasters_MA[[model]]))
  counts <- sapply(strong_crcm5_noreasters_MA[[model]], nrow)
  strong_crcm5_noreaster_MA_freq[[model]] <- data.frame(year = years, 
                                                        count = counts)
}
# Print summary
cat("\n=== SUMMARY ===\n")
for (model in names(strong_crcm5_noreaster_MA_freq)) {
  freq_data <- strong_crcm5_noreaster_MA_freq[[model]]
  cat(paste(toupper(model), ":", nrow(freq_data), "seasons,",
            "mean =", round(mean(freq_data$count), 1),
            "std =", round(sd(freq_data$count), 1), "\n"))
}

# ==============================================================================
#
# Create master CRCM5 frequency dataframes and calculate CRCM5 multi-model mean 
#   frequencies (full dataset, strong & weak subsets, & strong nor'easters 
#   for each season)

# Define years for time series
crcm5_years <- sort(unique(c(
  crcm5_noreaster_frequencies$canesm5$year,
  crcm5_noreaster_frequencies$cnrm$year,
  crcm5_noreaster_frequencies$noresm2$year,
  crcm5_noreaster_frequencies$mpi$year,
  crcm5_noreaster_frequencies$era5$year
)))

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Full dataset

# Build a combined data frame
crcm5_noreaster_mme_df <- data.frame(
  year = crcm5_years,
  canesm5 = NA,
  cnrm = NA,
  noresm2 = NA,
  mpi = NA,
  era5 = NA
)

# Fill in storm counts for each model
crcm5_noreaster_mme_df$canesm5 <- crcm5_noreaster_frequencies$canesm5$count[match(crcm5_years, 
                                                                                  crcm5_noreaster_frequencies$canesm5$year)]
crcm5_noreaster_mme_df$cnrm     <- crcm5_noreaster_frequencies$cnrm$count[match(crcm5_years, 
                                                                                crcm5_noreaster_frequencies$cnrm$year)]
crcm5_noreaster_mme_df$noresm2  <- crcm5_noreaster_frequencies$noresm2$count[match(crcm5_years, 
                                                                                   crcm5_noreaster_frequencies$noresm2$year)]
crcm5_noreaster_mme_df$mpi <- crcm5_noreaster_frequencies$mpi$count[match(crcm5_years, 
                                                                          crcm5_noreaster_frequencies$mpi$year)]
crcm5_noreaster_mme_df$era5 <- crcm5_noreaster_frequencies$era5$count[match(crcm5_years, 
                                                                            crcm5_noreaster_frequencies$mpi$year)]

# Calculate multi-model ensemble mean and std dev across the 4 models
crcm5_noreaster_mme_df$mean_count <- rowMeans(crcm5_noreaster_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                              na.rm = TRUE)
crcm5_noreaster_mme_df$std_count  <- apply(crcm5_noreaster_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                           1, sd, na.rm = TRUE)
#
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Weak subset
#
# Build a combined data frame
weak_crcm5_noreaster_mme_df <- data.frame(
  year = crcm5_years,
  canesm5 = NA,
  cnrm = NA,
  noresm2 = NA,
  mpi = NA,
  era5 = NA
)

# Fill in storm counts for each model
weak_crcm5_noreaster_mme_df$canesm5 <- weak_crcm5_noreaster_frequencies$canesm5$count[match(crcm5_years, 
                                                                                            weak_crcm5_noreaster_frequencies$canesm5$year)]
weak_crcm5_noreaster_mme_df$cnrm     <- weak_crcm5_noreaster_frequencies$cnrm$count[match(crcm5_years, 
                                                                                          weak_crcm5_noreaster_frequencies$cnrm$year)]
weak_crcm5_noreaster_mme_df$noresm2  <- weak_crcm5_noreaster_frequencies$noresm2$count[match(crcm5_years, 
                                                                                             weak_crcm5_noreaster_frequencies$noresm2$year)]
weak_crcm5_noreaster_mme_df$mpi <- weak_crcm5_noreaster_frequencies$mpi$count[match(crcm5_years, 
                                                                                    weak_crcm5_noreaster_frequencies$mpi$year)]
weak_crcm5_noreaster_mme_df$era5 <- weak_crcm5_noreaster_frequencies$era5$count[match(crcm5_years, 
                                                                                      weak_crcm5_noreaster_frequencies$mpi$year)]

# Calculate multi-model ensemble mean and std dev across the 4 models
weak_crcm5_noreaster_mme_df$mean_count <- rowMeans(weak_crcm5_noreaster_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                                   na.rm = TRUE)
weak_crcm5_noreaster_mme_df$std_count  <- apply(weak_crcm5_noreaster_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                                1, sd, na.rm = TRUE)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Strong subset

# Build a combined data frame
strong_crcm5_noreaster_mme_df <- data.frame(
  year = crcm5_years,
  canesm5 = NA,
  cnrm = NA,
  noresm2 = NA,
  mpi = NA,
  era5 = NA
)

# Fill in storm counts for each model
strong_crcm5_noreaster_mme_df$canesm5 <- strong_crcm5_noreaster_frequencies$canesm5$count[match(crcm5_years, 
                                                                                                strong_crcm5_noreaster_frequencies$canesm5$year)]
strong_crcm5_noreaster_mme_df$cnrm     <- strong_crcm5_noreaster_frequencies$cnrm$count[match(crcm5_years, 
                                                                                              strong_crcm5_noreaster_frequencies$cnrm$year)]
strong_crcm5_noreaster_mme_df$noresm2  <- strong_crcm5_noreaster_frequencies$noresm2$count[match(crcm5_years, 
                                                                                                 strong_crcm5_noreaster_frequencies$noresm2$year)]
strong_crcm5_noreaster_mme_df$mpi <- strong_crcm5_noreaster_frequencies$mpi$count[match(crcm5_years, 
                                                                                        strong_crcm5_noreaster_frequencies$mpi$year)]
strong_crcm5_noreaster_mme_df$era5 <- strong_crcm5_noreaster_frequencies$era5$count[match(crcm5_years, 
                                                                                          strong_crcm5_noreaster_frequencies$mpi$year)]

# Calculate multi-model ensemble mean and std dev across the 4 models
strong_crcm5_noreaster_mme_df$mean_count <- rowMeans(strong_crcm5_noreaster_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                                     na.rm = TRUE)
strong_crcm5_noreaster_mme_df$std_count  <- apply(strong_crcm5_noreaster_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                                  1, sd, na.rm = TRUE)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Strong SON subset

# Build a combined data frame
strong_crcm5_noreaster_SON_mme_df <- data.frame(
  year = crcm5_years,
  canesm5 = NA,
  cnrm = NA,
  noresm2 = NA,
  mpi = NA,
  era5 = NA
)

# Fill in storm counts for each model
strong_crcm5_noreaster_SON_mme_df$canesm5 <- strong_crcm5_noreaster_SON_freq$canesm5$count[match(crcm5_years, 
                                                                                                 strong_crcm5_noreaster_SON_freq$canesm5$year)]
strong_crcm5_noreaster_SON_mme_df$cnrm     <- strong_crcm5_noreaster_SON_freq$cnrm$count[match(crcm5_years, 
                                                                                               strong_crcm5_noreaster_SON_freq$cnrm$year)]
strong_crcm5_noreaster_SON_mme_df$noresm2  <- strong_crcm5_noreaster_SON_freq$noresm2$count[match(crcm5_years, 
                                                                                                  strong_crcm5_noreaster_SON_freq$noresm2$year)]
strong_crcm5_noreaster_SON_mme_df$mpi <- strong_crcm5_noreaster_SON_freq$mpi$count[match(crcm5_years, 
                                                                                         strong_crcm5_noreaster_SON_freq$mpi$year)]
strong_crcm5_noreaster_SON_mme_df$era5 <- strong_crcm5_noreaster_SON_freq$era5$count[match(crcm5_years, 
                                                                                           strong_crcm5_noreaster_SON_freq$mpi$year)]

# Calculate multi-model ensemble mean and std dev across the 4 models
strong_crcm5_noreaster_SON_mme_df$mean_count <- rowMeans(strong_crcm5_noreaster_SON_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                                         na.rm = TRUE)
strong_crcm5_noreaster_SON_mme_df$std_count  <- apply(strong_crcm5_noreaster_SON_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                                      1, sd, na.rm = TRUE)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Strong DJF subset

# Build a combined data frame
strong_crcm5_noreaster_DJF_mme_df <- data.frame(
  year = crcm5_years,
  canesm5 = NA,
  cnrm = NA,
  noresm2 = NA,
  mpi = NA,
  era5 = NA
)

# Fill in storm counts for each model
strong_crcm5_noreaster_DJF_mme_df$canesm5 <- strong_crcm5_noreaster_DJF_freq$canesm5$count[match(crcm5_years, 
                                                                                                 strong_crcm5_noreaster_DJF_freq$canesm5$year)]
strong_crcm5_noreaster_DJF_mme_df$cnrm     <- strong_crcm5_noreaster_DJF_freq$cnrm$count[match(crcm5_years, 
                                                                                               strong_crcm5_noreaster_DJF_freq$cnrm$year)]
strong_crcm5_noreaster_DJF_mme_df$noresm2  <- strong_crcm5_noreaster_DJF_freq$noresm2$count[match(crcm5_years, 
                                                                                                  strong_crcm5_noreaster_DJF_freq$noresm2$year)]
strong_crcm5_noreaster_DJF_mme_df$mpi <- strong_crcm5_noreaster_DJF_freq$mpi$count[match(crcm5_years, 
                                                                                         strong_crcm5_noreaster_DJF_freq$mpi$year)]
strong_crcm5_noreaster_DJF_mme_df$era5 <- strong_crcm5_noreaster_DJF_freq$era5$count[match(crcm5_years, 
                                                                                           strong_crcm5_noreaster_DJF_freq$mpi$year)]

# Calculate multi-model ensemble mean and std dev across the 4 models
strong_crcm5_noreaster_DJF_mme_df$mean_count <- rowMeans(strong_crcm5_noreaster_DJF_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                                         na.rm = TRUE)
strong_crcm5_noreaster_DJF_mme_df$std_count  <- apply(strong_crcm5_noreaster_DJF_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                                      1, sd, na.rm = TRUE)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Strong MA subset

# Build a combined data frame
strong_crcm5_noreaster_MA_mme_df <- data.frame(
  year = crcm5_years,
  canesm5 = NA,
  cnrm = NA,
  noresm2 = NA,
  mpi = NA,
  era5 = NA
)

# Fill in storm counts for each model
strong_crcm5_noreaster_MA_mme_df$canesm5 <- strong_crcm5_noreaster_MA_freq$canesm5$count[match(crcm5_years, 
                                                                                               strong_crcm5_noreaster_MA_freq$canesm5$year)]
strong_crcm5_noreaster_MA_mme_df$cnrm     <- strong_crcm5_noreaster_MA_freq$cnrm$count[match(crcm5_years, 
                                                                                             strong_crcm5_noreaster_MA_freq$cnrm$year)]
strong_crcm5_noreaster_MA_mme_df$noresm2  <- strong_crcm5_noreaster_MA_freq$noresm2$count[match(crcm5_years, 
                                                                                                strong_crcm5_noreaster_MA_freq$noresm2$year)]
strong_crcm5_noreaster_MA_mme_df$mpi <- strong_crcm5_noreaster_MA_freq$mpi$count[match(crcm5_years, 
                                                                                       strong_crcm5_noreaster_MA_freq$mpi$year)]
strong_crcm5_noreaster_MA_mme_df$era5 <- strong_crcm5_noreaster_MA_freq$era5$count[match(crcm5_years, 
                                                                                         strong_crcm5_noreaster_MA_freq$mpi$year)]

# Calculate multi-model ensemble mean and std dev across the 4 models
strong_crcm5_noreaster_MA_mme_df$mean_count <- rowMeans(strong_crcm5_noreaster_MA_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                                        na.rm = TRUE)
strong_crcm5_noreaster_MA_mme_df$std_count  <- apply(strong_crcm5_noreaster_MA_mme_df[, c("canesm5", "cnrm", "noresm2", "mpi")], 
                                                     1, sd, na.rm = TRUE)

# ==============================================================================

# Extract intensity values (lifetime min. SLP & max. surface wind) from each dataset

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CFSR, 1979-2014

# Full dataset
cfsr_min_slp_1979_2014 <- process_era5_for_intensity(cfsr_strong_noreasters_1979_2014, 
                                                     metric = "min_slp")
cfsr_max_wind_1979_2014 <- process_era5_for_intensity(cfsr_strong_noreasters_1979_2014, 
                                                      metric = "max_wind")
# SON
cfsr_min_slp_1979_2014_SON <- process_era5_for_intensity(cfsr_strong_noreasters_1979_2014_SON, 
                                                         metric = "min_slp")
cfsr_max_wind_1979_2014_SON <- process_era5_for_intensity(cfsr_strong_noreasters_1979_2014_SON, 
                                                          metric = "max_wind")
# DJF
cfsr_min_slp_1979_2014_DJF <- process_era5_for_intensity(cfsr_strong_noreasters_1979_2014_DJF, 
                                                         metric = "min_slp")
cfsr_max_wind_1979_2014_DJF <- process_era5_for_intensity(cfsr_strong_noreasters_1979_2014_DJF, 
                                                          metric = "max_wind")
# MA
cfsr_min_slp_1979_2014_MA <- process_era5_for_intensity(cfsr_strong_noreasters_1979_2014_MA, 
                                                        metric = "min_slp")
cfsr_max_wind_1979_2014_MA <- process_era5_for_intensity(cfsr_strong_noreasters_1979_2014_MA, 
                                                         metric = "max_wind")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# ERA5

# 1950-2014
era5_min_slp_1950_2014 <- process_era5_for_intensity(era5_strong_noreasters_1950_2014, 
                                                            metric = "min_slp")
era5_max_wind_1950_2014 <- process_era5_for_intensity(era5_strong_noreasters_1950_2014, 
                                                             metric = "max_wind")
# SON
era5_min_slp_1950_2014_SON <- process_era5_for_intensity(era5_strong_noreasters_1950_2014_SON, 
                                                                metric = "min_slp")
era5_max_wind_1950_2014_SON <- process_era5_for_intensity(era5_strong_noreasters_1950_2014_SON, 
                                                                 metric = "max_wind")
# DJF
era5_min_slp_1950_2014_DJF <- process_era5_for_intensity(era5_strong_noreasters_1950_2014_DJF, 
                                                                metric = "min_slp")
era5_max_wind_1950_2014_DJF <- process_era5_for_intensity(era5_strong_noreasters_1950_2014_DJF, 
                                                                 metric = "max_wind")
# MA
era5_min_slp_1950_2014_MA <- process_era5_for_intensity(era5_strong_noreasters_1950_2014_MA, 
                                                               metric = "min_slp")
era5_max_wind_1950_2014_MA <- process_era5_for_intensity(era5_strong_noreasters_1950_2014_MA, 
                                                                metric = "max_wind")

# 1979-2014
era5_min_slp_1979_2014 <- process_era5_for_intensity(era5_strong_noreasters_1979_2014, 
                                                            metric = "min_slp")
era5_max_wind_1979_2014 <- process_era5_for_intensity(era5_strong_noreasters_1979_2014, 
                                                             metric = "max_wind")
# SON
era5_min_slp_1979_2014_SON <- process_era5_for_intensity(era5_strong_noreasters_1979_2014_SON, 
                                                                metric = "min_slp")
era5_max_wind_1979_2014_SON <- process_era5_for_intensity(era5_strong_noreasters_1979_2014_SON, 
                                                                 metric = "max_wind")
# DJF
era5_min_slp_1979_2014_DJF <- process_era5_for_intensity(era5_strong_noreasters_1979_2014_DJF, 
                                                                metric = "min_slp")
era5_max_wind_1979_2014_DJF <- process_era5_for_intensity(era5_strong_noreasters_1979_2014_DJF, 
                                                                 metric = "max_wind")
# MA
era5_min_slp_1979_2014_MA <- process_era5_for_intensity(era5_strong_noreasters_1979_2014_MA, 
                                                               metric = "min_slp")
era5_max_wind_1979_2014_MA <- process_era5_for_intensity(era5_strong_noreasters_1979_2014_MA, 
                                                                metric = "max_wind")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-ERA5, 1979-2014

# Full dataset
crcm5_era5_min_slp_1979_2014 <- extract_intensity_values(strong_crcm5_noreasters_1979_2014["era5"], 
                                                         metric = "min_slp")
crcm5_era5_min_slp_1979_2014 <- crcm5_era5_min_slp_1979_2014$era5     # remove unnecessary data frame level
crcm5_era5_max_wind_1979_2014 <- extract_intensity_values(strong_crcm5_noreasters_1979_2014["era5"], 
                                                          metric = "max_wind")
crcm5_era5_max_wind_1979_2014 <- crcm5_era5_max_wind_1979_2014$era5
# SON
crcm5_era5_min_slp_1979_2014_SON <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_SON["era5"], 
                                                             metric = "min_slp")
crcm5_era5_min_slp_1979_2014_SON <- crcm5_era5_min_slp_1979_2014_SON$era5
crcm5_era5_max_wind_1979_2014_SON <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_SON["era5"], 
                                                              metric = "max_wind")
crcm5_era5_max_wind_1979_2014_SON <- crcm5_era5_max_wind_1979_2014_SON$era5
# DJF
crcm5_era5_min_slp_1979_2014_DJF <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_DJF["era5"], 
                                                             metric = "min_slp")
crcm5_era5_min_slp_1979_2014_DJF <- crcm5_era5_min_slp_1979_2014_DJF$era5
crcm5_era5_max_wind_1979_2014_DJF <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_DJF["era5"], 
                                                              metric = "max_wind")
crcm5_era5_max_wind_1979_2014_DJF <- crcm5_era5_max_wind_1979_2014_DJF$era5
# MA
crcm5_era5_min_slp_1979_2014_MA <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_MA["era5"], 
                                                            metric = "min_slp")
crcm5_era5_min_slp_1979_2014_MA <- crcm5_era5_min_slp_1979_2014_MA$era5
crcm5_era5_max_wind_1979_2014_MA <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_MA["era5"], 
                                                             metric = "max_wind")
crcm5_era5_max_wind_1979_2014_MA <- crcm5_era5_max_wind_1979_2014_MA$era5

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-CMIP6

# 1950-2014
crcm5_cmip6_min_slp_1950_2014 <- extract_intensity_values(strong_crcm5_noreasters[gcm_models], 
                                                                metric = "min_slp")
crcm5_cmip6_max_wind_1950_2014 <- extract_intensity_values(strong_crcm5_noreasters[gcm_models], 
                                                                 metric = "max_wind")
# SON
crcm5_cmip6_min_slp_1950_2014_SON <- extract_intensity_values(strong_crcm5_noreasters_SON[gcm_models], 
                                                                    metric = "min_slp")
crcm5_cmip6_max_wind_1950_2014_SON <- extract_intensity_values(strong_crcm5_noreasters_SON[gcm_models], 
                                                                     metric = "max_wind")
# DJF
crcm5_cmip6_min_slp_1950_2014_DJF <- extract_intensity_values(strong_crcm5_noreasters_DJF[gcm_models], 
                                                                    metric = "min_slp")
crcm5_cmip6_max_wind_1950_2014_DJF <- extract_intensity_values(strong_crcm5_noreasters_DJF[gcm_models], 
                                                                     metric = "max_wind")
# MA
crcm5_cmip6_min_slp_1950_2014_MA <- extract_intensity_values(strong_crcm5_noreasters_MA[gcm_models], 
                                                                   metric = "min_slp")
crcm5_cmip6_max_wind_1950_2014_MA <- extract_intensity_values(strong_crcm5_noreasters_MA[gcm_models], 
                                                                    metric = "max_wind")

# 1979-2014
crcm5_cmip6_min_slp_1979_2014 <- extract_intensity_values(strong_crcm5_noreasters_1979_2014[gcm_models], 
                                                                metric = "min_slp")
crcm5_cmip6_max_wind_1979_2014 <- extract_intensity_values(strong_crcm5_noreasters_1979_2014[gcm_models], 
                                                                 metric = "max_wind")
# SON
crcm5_cmip6_min_slp_1979_2014_SON <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_SON[gcm_models], 
                                                                    metric = "min_slp")
crcm5_cmip6_max_wind_1979_2014_SON <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_SON[gcm_models], 
                                                                     metric = "max_wind")
# DJF
crcm5_cmip6_min_slp_1979_2014_DJF <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_DJF[gcm_models], 
                                                                    metric = "min_slp")
crcm5_cmip6_max_wind_1979_2014_DJF <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_DJF[gcm_models], 
                                                                     metric = "max_wind")
# MA
crcm5_cmip6_min_slp_1979_2014_MA <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_MA[gcm_models], 
                                                                   metric = "min_slp")
crcm5_cmip6_max_wind_1979_2014_MA <- extract_intensity_values(strong_crcm5_noreasters_1979_2014_MA[gcm_models], 
                                                                    metric = "max_wind")

# ==============================================================================

# Calculate pooled intensity quantiles for each dataset, 1979-2014

my.probs = seq(from = 0, to = 1, by = 0.01)

# CFSR
cfsr_min_slp_quants_1979_2014 <- calculate_cfsr_intensity_quantiles(cfsr_min_slp_1979_2014, metric = "min_slp", my.probs)
cfsr_max_wind_quants_1979_2014 <- calculate_cfsr_intensity_quantiles(cfsr_max_wind_1979_2014, metric = "max_wind", my.probs)

# ERA5
era5_min_slp_quants_1979_2014 <- calculate_era5_intensity_quantiles(era5_min_slp_1979_2014, metric = "min_slp", my.probs)
era5_max_wind_quants_1979_2014 <- calculate_era5_intensity_quantiles(era5_max_wind_1979_2014, metric = "max_wind", my.probs)

# CRCM5-ERA5 (same function as CFSR)
crcm5_era5_min_slp_quants_1979_2014 <- calculate_cfsr_intensity_quantiles(crcm5_era5_min_slp_1979_2014, metric = "min_slp", my.probs)
crcm5_era5_max_wind_quants_1979_2014 <- calculate_cfsr_intensity_quantiles(crcm5_era5_max_wind_1979_2014, metric = "max_wind", my.probs)

# CRCM5-CMIP6
crcm5_cmip6_min_slp_quants_1979_2014 <- calculate_ensemble_intensity_quantiles(crcm5_cmip6_min_slp_1979_2014, metric = "min_slp", my.probs)
crcm5_cmip6_max_wind_quants_1979_2014 <- calculate_ensemble_intensity_quantiles(crcm5_cmip6_max_wind_1979_2014, metric = "max_wind", my.probs)

# ==============================================================================

# Calculate intensity quantiles over time & trends for each dataset

my.probs = seq(from = 0, to = 1, by = 0.01)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CFSR, 1979-2014

# Full dataset
cfsr_slp_quant_trends_1979_2014 <- plot_era5_intensity_quantile_trends(cfsr_min_slp_1979_2014, # plots time series of each quantile
                                                                       "min_slp", 
                                                                       quantiles = my.probs)
cfsr_wind_quant_trends_1979_2014 <- plot_era5_intensity_quantile_trends(cfsr_max_wind_1979_2014, 
                                                                        "max_wind", 
                                                                        quantiles = my.probs)
# SON
cfsr_slp_quant_trends_1979_2014_SON <- plot_era5_intensity_quantile_trends(cfsr_min_slp_1979_2014_SON, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
cfsr_wind_quant_trends_1979_2014_SON <- plot_era5_intensity_quantile_trends(cfsr_max_wind_1979_2014_SON, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# DJF
cfsr_slp_quant_trends_1979_2014_DJF <- plot_era5_intensity_quantile_trends(cfsr_min_slp_1979_2014_DJF, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
cfsr_wind_quant_trends_1979_2014_DJF <- plot_era5_intensity_quantile_trends(cfsr_max_wind_1979_2014_DJF, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# MA
cfsr_slp_quant_trends_1979_2014_MA <- plot_era5_intensity_quantile_trends(cfsr_min_slp_1979_2014_MA, 
                                                                          "min_slp", 
                                                                          quantiles = my.probs)
cfsr_wind_quant_trends_1979_2014_MA <- plot_era5_intensity_quantile_trends(cfsr_max_wind_1979_2014_MA, 
                                                                           "max_wind", 
                                                                           quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# ERA5

# 1950-2014
era5_slp_quant_trends_1950_2014 <- plot_era5_intensity_quantile_trends(era5_min_slp_1950_2014,
                                                                       "min_slp", 
                                                                       quantiles = my.probs)
era5_wind_quant_trends_1950_2014 <- plot_era5_intensity_quantile_trends(era5_max_wind_1950_2014, 
                                                                        "max_wind", 
                                                                        quantiles = my.probs)
# SON
era5_slp_quant_trends_1950_2014_SON <- plot_era5_intensity_quantile_trends(era5_min_slp_1950_2014_SON, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
era5_wind_quant_trends_1950_2014_SON <- plot_era5_intensity_quantile_trends(era5_max_wind_1950_2014_SON, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# DJF
era5_slp_quant_trends_1950_2014_DJF <- plot_era5_intensity_quantile_trends(era5_min_slp_1950_2014_DJF, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
era5_wind_quant_trends_1950_2014_DJF <- plot_era5_intensity_quantile_trends(era5_max_wind_1950_2014_DJF, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# MA
era5_slp_quant_trends_1950_2014_MA <- plot_era5_intensity_quantile_trends(era5_min_slp_1950_2014_MA, 
                                                                          "min_slp", 
                                                                          quantiles = my.probs)
era5_wind_quant_trends_1950_2014_MA <- plot_era5_intensity_quantile_trends(era5_max_wind_1950_2014_MA, 
                                                                           "max_wind", 
                                                                           quantiles = my.probs)

# 1979-2014
era5_slp_quant_trends_1979_2014 <- plot_era5_intensity_quantile_trends(era5_min_slp_1979_2014,
                                                                       "min_slp", 
                                                                       quantiles = my.probs)
era5_wind_quant_trends_1979_2014 <- plot_era5_intensity_quantile_trends(era5_max_wind_1979_2014, 
                                                                        "max_wind", 
                                                                        quantiles = my.probs)
# SON
era5_slp_quant_trends_1979_2014_SON <- plot_era5_intensity_quantile_trends(era5_min_slp_1979_2014_SON, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
era5_wind_quant_trends_1979_2014_SON <- plot_era5_intensity_quantile_trends(era5_max_wind_1979_2014_SON, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# DJF
era5_slp_quant_trends_1979_2014_DJF <- plot_era5_intensity_quantile_trends(era5_min_slp_1979_2014_DJF, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
era5_wind_quant_trends_1979_2014_DJF <- plot_era5_intensity_quantile_trends(era5_max_wind_1979_2014_DJF, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# MA
era5_slp_quant_trends_1979_2014_MA <- plot_era5_intensity_quantile_trends(era5_min_slp_1979_2014_MA, 
                                                                          "min_slp", 
                                                                          quantiles = my.probs)
era5_wind_quant_trends_1979_2014_MA <- plot_era5_intensity_quantile_trends(era5_max_wind_1979_2014_MA, 
                                                                           "max_wind", 
                                                                           quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-ERA5, 1979-2014

# Full dataset
crcm5_era5_slp_quant_trends_1979_2014 <- plot_era5_intensity_quantile_trends(crcm5_era5_min_slp_1979_2014, # plots time series of each quantile
                                                                       "min_slp", 
                                                                       quantiles = my.probs)
crcm5_era5_wind_quant_trends_1979_2014 <- plot_era5_intensity_quantile_trends(crcm5_era5_max_wind_1979_2014, 
                                                                        "max_wind", 
                                                                        quantiles = my.probs)
# SON
crcm5_era5_slp_quant_trends_1979_2014_SON <- plot_era5_intensity_quantile_trends(crcm5_era5_min_slp_1979_2014_SON, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
crcm5_era5_wind_quant_trends_1979_2014_SON <- plot_era5_intensity_quantile_trends(crcm5_era5_max_wind_1979_2014_SON, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# DJF
crcm5_era5_slp_quant_trends_1979_2014_DJF <- plot_era5_intensity_quantile_trends(crcm5_era5_min_slp_1979_2014_DJF, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
crcm5_era5_wind_quant_trends_1979_2014_DJF <- plot_era5_intensity_quantile_trends(crcm5_era5_max_wind_1979_2014_DJF, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# MA
crcm5_era5_slp_quant_trends_1979_2014_MA <- plot_era5_intensity_quantile_trends(crcm5_era5_min_slp_1979_2014_MA, 
                                                                          "min_slp", 
                                                                          quantiles = my.probs)
crcm5_era5_wind_quant_trends_1979_2014_MA <- plot_era5_intensity_quantile_trends(crcm5_era5_max_wind_1979_2014_MA, 
                                                                           "max_wind", 
                                                                           quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-CMIP6 ensemble (pooled pilot models)

# 1950-2014
crcm5_cmip6_ensemble_slp_quant_trends_1950_2014 <- plot_ensemble_quantile_trends(crcm5_cmip6_min_slp_1950_2014, 
                                                                                 metric = "min_slp", 
                                                                                 quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_1950_2014 <- plot_ensemble_quantile_trends(crcm5_cmip6_max_wind_1950_2014, 
                                                                                  metric = "max_wind", 
                                                                                  quantiles = my.probs)
# SON
crcm5_cmip6_ensemble_slp_quant_trends_1950_2014_SON <- plot_ensemble_quantile_trends(crcm5_cmip6_min_slp_1950_2014_SON, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_1950_2014_SON <- plot_ensemble_quantile_trends(crcm5_cmip6_max_wind_1950_2014_SON, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# DJF
crcm5_cmip6_ensemble_slp_quant_trends_1950_2014_DJF <- plot_ensemble_quantile_trends(crcm5_cmip6_min_slp_1950_2014_DJF, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_1950_2014_DJF <- plot_ensemble_quantile_trends(crcm5_cmip6_max_wind_1950_2014_DJF, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# MA
crcm5_cmip6_ensemble_slp_quant_trends_1950_2014_MA <- plot_ensemble_quantile_trends(crcm5_cmip6_min_slp_1950_2014_MA, 
                                                                                    metric = "min_slp", 
                                                                                    quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_1950_2014_MA <- plot_ensemble_quantile_trends(crcm5_cmip6_max_wind_1950_2014_MA, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)

# 1979-2014
crcm5_cmip6_ensemble_slp_quant_trends_1979_2014 <- plot_ensemble_quantile_trends(crcm5_cmip6_min_slp_1979_2014, 
                                                                                 metric = "min_slp", 
                                                                                 quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_1979_2014 <- plot_ensemble_quantile_trends(crcm5_cmip6_max_wind_1979_2014, 
                                                                                  metric = "max_wind", 
                                                                                  quantiles = my.probs)
# SON
crcm5_cmip6_ensemble_slp_quant_trends_1979_2014_SON <- plot_ensemble_quantile_trends(crcm5_cmip6_min_slp_1979_2014_SON, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_1979_2014_SON <- plot_ensemble_quantile_trends(crcm5_cmip6_max_wind_1979_2014_SON, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# DJF
crcm5_cmip6_ensemble_slp_quant_trends_1979_2014_DJF <- plot_ensemble_quantile_trends(crcm5_cmip6_min_slp_1979_2014_DJF, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_1979_2014_DJF <- plot_ensemble_quantile_trends(crcm5_cmip6_max_wind_1979_2014_DJF, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# MA
crcm5_cmip6_ensemble_slp_quant_trends_1979_2014_MA <- plot_ensemble_quantile_trends(crcm5_cmip6_min_slp_1979_2014_MA, 
                                                                                    metric = "min_slp", 
                                                                                    quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_1979_2014_MA <- plot_ensemble_quantile_trends(crcm5_cmip6_max_wind_1979_2014_MA, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-CMIP6 individual models (saved in data frames)

# 1950-2014
crcm5_cmip6_models_slp_quant_trends_1950_2014 <- plot_intensity_quantile_trends(crcm5_cmip6_min_slp_1950_2014, 
                                                                                 metric = "min_slp", 
                                                                                 quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_1950_2014 <- plot_intensity_quantile_trends(crcm5_cmip6_max_wind_1950_2014, 
                                                                                  metric = "max_wind", 
                                                                                  quantiles = my.probs)
# SON
crcm5_cmip6_models_slp_quant_trends_1950_2014_SON <- plot_intensity_quantile_trends(crcm5_cmip6_min_slp_1950_2014_SON, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_1950_2014_SON <- plot_intensity_quantile_trends(crcm5_cmip6_max_wind_1950_2014_SON, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# DJF
crcm5_cmip6_models_slp_quant_trends_1950_2014_DJF <- plot_intensity_quantile_trends(crcm5_cmip6_min_slp_1950_2014_DJF, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_1950_2014_DJF <- plot_intensity_quantile_trends(crcm5_cmip6_max_wind_1950_2014_DJF, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# MA
crcm5_cmip6_models_slp_quant_trends_1950_2014_MA <- plot_intensity_quantile_trends(crcm5_cmip6_min_slp_1950_2014_MA, 
                                                                                    metric = "min_slp", 
                                                                                    quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_1950_2014_MA <- plot_intensity_quantile_trends(crcm5_cmip6_max_wind_1950_2014_MA, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)

# 1979-2014
crcm5_cmip6_models_slp_quant_trends_1979_2014 <- plot_intensity_quantile_trends(crcm5_cmip6_min_slp_1979_2014, 
                                                                                 metric = "min_slp", 
                                                                                 quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_1979_2014 <- plot_intensity_quantile_trends(crcm5_cmip6_max_wind_1979_2014, 
                                                                                  metric = "max_wind", 
                                                                                  quantiles = my.probs)
# SON
crcm5_cmip6_models_slp_quant_trends_1979_2014_SON <- plot_intensity_quantile_trends(crcm5_cmip6_min_slp_1979_2014_SON, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_1979_2014_SON <- plot_intensity_quantile_trends(crcm5_cmip6_max_wind_1979_2014_SON, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# DJF
crcm5_cmip6_models_slp_quant_trends_1979_2014_DJF <- plot_intensity_quantile_trends(crcm5_cmip6_min_slp_1979_2014_DJF, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_1979_2014_DJF <- plot_intensity_quantile_trends(crcm5_cmip6_max_wind_1979_2014_DJF, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# MA
crcm5_cmip6_models_slp_quant_trends_1979_2014_MA <- plot_intensity_quantile_trends(crcm5_cmip6_min_slp_1979_2014_MA, 
                                                                                    metric = "min_slp", 
                                                                                    quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_1979_2014_MA <- plot_intensity_quantile_trends(crcm5_cmip6_max_wind_1979_2014_MA, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)

# ==============================================================================

# Calculate intensification rates for each dataset

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CFSR, 1979-2014

# Full dataset
cfsr_intens_rates_slp_1979_2014 <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_1979_2014, 
                                                                           metric ="min_slp", 
                                                                           time_window_hours = 24, 
                                                                           timestep_hours = 6)
cfsr_intens_rates_wind_1979_2014 <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_1979_2014, 
                                                                            metric ="max_wind", 
                                                                            time_window_hours = 24, 
                                                                            timestep_hours = 6)
# SON
cfsr_intens_rates_slp_1979_2014_SON <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_1979_2014_SON, 
                                                                               metric ="min_slp", 
                                                                               time_window_hours = 24, 
                                                                               timestep_hours = 6)
cfsr_intens_rates_wind_1979_2014_SON <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_1979_2014_SON, 
                                                                                metric ="max_wind", 
                                                                                time_window_hours = 24, 
                                                                                timestep_hours = 6)
# DJF
cfsr_intens_rates_slp_1979_2014_DJF <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_1979_2014_DJF, 
                                                                               metric ="min_slp", 
                                                                               time_window_hours = 24, 
                                                                               timestep_hours = 6)
cfsr_intens_rates_wind_1979_2014_DJF <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_1979_2014_DJF, 
                                                                                metric ="max_wind", 
                                                                                time_window_hours = 24, 
                                                                                timestep_hours = 6)
# MA
cfsr_intens_rates_slp_1979_2014_MA <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_1979_2014_MA, 
                                                                              metric ="min_slp", 
                                                                              time_window_hours = 24, 
                                                                              timestep_hours = 6)
cfsr_intens_rates_wind_1979_2014_MA <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_1979_2014_MA, 
                                                                               metric ="max_wind", 
                                                                               time_window_hours = 24, 
                                                                               timestep_hours = 6)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# ERA5

# 1950-2014
era5_intens_rates_slp_1950_2014 <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1950_2014, 
                                                                           metric ="min_slp", 
                                                                           time_window_hours = 24, 
                                                                           timestep_hours = 1)
era5_intens_rates_wind_1950_2014 <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1950_2014, 
                                                                            metric ="max_wind", 
                                                                            time_window_hours = 24, 
                                                                            timestep_hours = 1)
# SON
era5_intens_rates_slp_1950_2014_SON <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1950_2014_SON, 
                                                                               metric ="min_slp", 
                                                                               time_window_hours = 24, 
                                                                               timestep_hours = 1)
era5_intens_rates_wind_1950_2014_SON <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1950_2014_SON, 
                                                                                metric ="max_wind", 
                                                                                time_window_hours = 24, 
                                                                                timestep_hours = 1)
# DJF
era5_intens_rates_slp_1950_2014_DJF <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1950_2014_DJF, 
                                                                               metric ="min_slp", 
                                                                               time_window_hours = 24, 
                                                                               timestep_hours = 1)
era5_intens_rates_wind_1950_2014_DJF <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1950_2014_DJF, 
                                                                                metric ="max_wind", 
                                                                                time_window_hours = 24, 
                                                                                timestep_hours = 1)
# MA
era5_intens_rates_slp_1950_2014_MA <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1950_2014_MA, 
                                                                              metric ="min_slp", 
                                                                              time_window_hours = 24, 
                                                                              timestep_hours = 1)
era5_intens_rates_wind_1950_2014_MA <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1950_2014_MA, 
                                                                               metric ="max_wind", 
                                                                               time_window_hours = 24, 
                                                                               timestep_hours = 1)

# 1979-2014
era5_intens_rates_slp_1979_2014 <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1979_2014, 
                                                                           metric ="min_slp", 
                                                                           time_window_hours = 24, 
                                                                           timestep_hours = 1)
era5_intens_rates_wind_1979_2014 <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1979_2014, 
                                                                            metric ="max_wind", 
                                                                            time_window_hours = 24, 
                                                                            timestep_hours = 1)
# SON
era5_intens_rates_slp_1979_2014_SON <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1979_2014_SON, 
                                                                               metric ="min_slp", 
                                                                               time_window_hours = 24, 
                                                                               timestep_hours = 1)
era5_intens_rates_wind_1979_2014_SON <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1979_2014_SON, 
                                                                                metric ="max_wind", 
                                                                                time_window_hours = 24, 
                                                                                timestep_hours = 1)
# DJF
era5_intens_rates_slp_1979_2014_DJF <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1979_2014_DJF, 
                                                                               metric ="min_slp", 
                                                                               time_window_hours = 24, 
                                                                               timestep_hours = 1)
era5_intens_rates_wind_1979_2014_DJF <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1979_2014_DJF, 
                                                                                metric ="max_wind", 
                                                                                time_window_hours = 24, 
                                                                                timestep_hours = 1)
# MA
era5_intens_rates_slp_1979_2014_MA <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1979_2014_MA, 
                                                                              metric ="min_slp", 
                                                                              time_window_hours = 24, 
                                                                              timestep_hours = 1)
era5_intens_rates_wind_1979_2014_MA <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1979_2014_MA, 
                                                                               metric ="max_wind", 
                                                                               time_window_hours = 24, 
                                                                               timestep_hours = 1)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-ERA5, 1979-2014

crcm5_era5_intens_rates_slp_1979_2014 <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014["era5"], 
                                                                            metric ="min_slp", time_window_hours = 24)
crcm5_era5_intens_rates_slp_1979_2014 <- crcm5_era5_intens_rates_slp_1979_2014$era5
crcm5_era5_intens_rates_wind_1979_2014 <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014["era5"], 
                                                                            metric ="max_wind", time_window_hours = 24)
crcm5_era5_intens_rates_wind_1979_2014 <- crcm5_era5_intens_rates_wind_1979_2014$era5
# SON
crcm5_era5_intens_rates_slp_1979_2014_SON <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_SON["era5"], 
                                                                            metric ="min_slp", time_window_hours = 24)
crcm5_era5_intens_rates_slp_1979_2014_SON <- crcm5_era5_intens_rates_slp_1979_2014_SON$era5
crcm5_era5_intens_rates_wind_1979_2014_SON <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_SON["era5"], 
                                                                             metric ="max_wind", time_window_hours = 24)
crcm5_era5_intens_rates_wind_1979_2014_SON <- crcm5_era5_intens_rates_wind_1979_2014_SON$era5
# DJF
crcm5_era5_intens_rates_slp_1979_2014_DJF <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_DJF["era5"], 
                                                                                metric ="min_slp", time_window_hours = 24)
crcm5_era5_intens_rates_slp_1979_2014_DJF <- crcm5_era5_intens_rates_slp_1979_2014_DJF$era5
crcm5_era5_intens_rates_wind_1979_2014_DJF <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_DJF["era5"], 
                                                                                 metric ="max_wind", time_window_hours = 24)
crcm5_era5_intens_rates_wind_1979_2014_DJF <- crcm5_era5_intens_rates_wind_1979_2014_DJF$era5
# MA
crcm5_era5_intens_rates_slp_1979_2014_MA <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_MA["era5"], 
                                                                                metric ="min_slp", time_window_hours = 24)
crcm5_era5_intens_rates_slp_1979_2014_MA <- crcm5_era5_intens_rates_slp_1979_2014_MA$era5
crcm5_era5_intens_rates_wind_1979_2014_MA <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_MA["era5"], 
                                                                                 metric ="max_wind", time_window_hours = 24)
crcm5_era5_intens_rates_wind_1979_2014_MA <- crcm5_era5_intens_rates_wind_1979_2014_MA$era5

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-CMIP6 ensemble (pooled pilot models)

# 1950-2014
crcm5_cmip6_intens_rates_slp_1950_2014 <- calculate_max_intensification_rate(strong_crcm5_noreasters[gcm_models], 
                                                                              metric ="min_slp", 
                                                                              time_window_hours = 24)
crcm5_cmip6_intens_rates_wind_1950_2014 <- calculate_max_intensification_rate(strong_crcm5_noreasters[gcm_models], 
                                                                              metric ="max_wind", 
                                                                              time_window_hours = 24)
# SON
crcm5_cmip6_intens_rates_slp_1950_2014_SON <- calculate_max_intensification_rate(strong_crcm5_noreasters_SON[gcm_models], 
                                                                             metric ="min_slp", 
                                                                             time_window_hours = 24)
crcm5_cmip6_intens_rates_wind_1950_2014_SON <- calculate_max_intensification_rate(strong_crcm5_noreasters_SON[gcm_models], 
                                                                              metric ="max_wind", 
                                                                              time_window_hours = 24)
# DJF
crcm5_cmip6_intens_rates_slp_1950_2014_DJF <- calculate_max_intensification_rate(strong_crcm5_noreasters_DJF[gcm_models], 
                                                                                 metric ="min_slp", 
                                                                                 time_window_hours = 24)
crcm5_cmip6_intens_rates_wind_1950_2014_DJF <- calculate_max_intensification_rate(strong_crcm5_noreasters_DJF[gcm_models], 
                                                                                  metric ="max_wind", 
                                                                                  time_window_hours = 24)
# MA
crcm5_cmip6_intens_rates_slp_1950_2014_MA <- calculate_max_intensification_rate(strong_crcm5_noreasters_MA[gcm_models], 
                                                                                 metric ="min_slp", 
                                                                                 time_window_hours = 24)
crcm5_cmip6_intens_rates_wind_1950_2014_MA <- calculate_max_intensification_rate(strong_crcm5_noreasters_MA[gcm_models], 
                                                                                  metric ="max_wind", 
                                                                                  time_window_hours = 24)

# 1979-2014
crcm5_cmip6_intens_rates_slp_1979_2014 <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014[gcm_models], 
                                                                             metric ="min_slp", 
                                                                             time_window_hours = 24)
crcm5_cmip6_intens_rates_wind_1979_2014 <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014[gcm_models], 
                                                                              metric ="max_wind", 
                                                                              time_window_hours = 24)
# SON
crcm5_cmip6_intens_rates_slp_1979_2014_SON <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_SON[gcm_models], 
                                                                                 metric ="min_slp", 
                                                                                 time_window_hours = 24)
crcm5_cmip6_intens_rates_wind_1979_2014_SON <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_SON[gcm_models], 
                                                                                  metric ="max_wind", 
                                                                                  time_window_hours = 24)
# DJF
crcm5_cmip6_intens_rates_slp_1979_2014_DJF <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_DJF[gcm_models], 
                                                                                 metric ="min_slp", 
                                                                                 time_window_hours = 24)
crcm5_cmip6_intens_rates_wind_1979_2014_DJF <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_DJF[gcm_models], 
                                                                                  metric ="max_wind", 
                                                                                  time_window_hours = 24)
# MA
crcm5_cmip6_intens_rates_slp_1979_2014_MA <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_MA[gcm_models], 
                                                                                metric ="min_slp", 
                                                                                time_window_hours = 24)
crcm5_cmip6_intens_rates_wind_1979_2014_MA <- calculate_max_intensification_rate(strong_crcm5_noreasters_1979_2014_MA[gcm_models], 
                                                                                 metric ="max_wind", 
                                                                                 time_window_hours = 24)

# ==============================================================================

# Calculate intensification rate quantiles

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CFSR, 1979-2014

# Full dataset
cfsr_slp_intens_quant_trends_1979_2014 <- plot_era5_intensification_quantile_trends(cfsr_intens_rates_slp_1979_2014, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
cfsr_wind_intens_quant_trends_1979_2014 <- plot_era5_intensification_quantile_trends(cfsr_intens_rates_wind_1979_2014, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# SON
cfsr_slp_intens_quant_trends_1979_2014_SON <- plot_era5_intensification_quantile_trends(cfsr_intens_rates_slp_1979_2014_SON, 
                                                                                            "min_slp", 
                                                                                            quantiles = my.probs)
cfsr_wind_intens_quant_trends_1979_2014_SON <- plot_era5_intensification_quantile_trends(cfsr_intens_rates_wind_1979_2014_SON, 
                                                                                             "max_wind", 
                                                                                             quantiles = my.probs)
# DJF
cfsr_slp_intens_quant_trends_1979_2014_DJF <- plot_era5_intensification_quantile_trends(cfsr_intens_rates_slp_1979_2014_DJF, 
                                                                                            "min_slp", 
                                                                                            quantiles = my.probs)
cfsr_wind_intens_quant_trends_1979_2014_DJF <- plot_era5_intensification_quantile_trends(cfsr_intens_rates_wind_1979_2014_DJF, 
                                                                                             "max_wind", 
                                                                                             quantiles = my.probs)
# MA
cfsr_slp_intens_quant_trends_1979_2014_MA <- plot_era5_intensification_quantile_trends(cfsr_intens_rates_slp_1979_2014_MA, 
                                                                                           "min_slp", 
                                                                                           quantiles = my.probs)
cfsr_wind_intens_quant_trends_1979_2014_MA <- plot_era5_intensification_quantile_trends(cfsr_intens_rates_wind_1979_2014_MA, 
                                                                                            "max_wind", 
                                                                                            quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# ERA5

# 1950-2014
era5_slp_intens_quant_trends_1950_2014 <- plot_era5_intensification_quantile_trends(era5_intens_rates_slp_1950_2014, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
era5_wind_intens_quant_trends_1950_2014 <- plot_era5_intensification_quantile_trends(era5_intens_rates_wind_1950_2014, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# SON
era5_slp_intens_quant_trends_1950_2014_SON <- plot_era5_intensification_quantile_trends(era5_intens_rates_slp_1950_2014_SON, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
era5_wind_intens_quant_trends_1950_2014_SON <- plot_era5_intensification_quantile_trends(era5_intens_rates_wind_1950_2014_SON, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# DJF
era5_slp_intens_quant_trends_1950_2014_DJF <- plot_era5_intensification_quantile_trends(era5_intens_rates_slp_1950_2014_DJF, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
era5_wind_intens_quant_trends_1950_2014_DJF <- plot_era5_intensification_quantile_trends(era5_intens_rates_wind_1950_2014_DJF, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# MA
era5_slp_intens_quant_trends_1950_2014_MA <- plot_era5_intensification_quantile_trends(era5_intens_rates_slp_1950_2014_MA, 
                                                                                       "min_slp", 
                                                                                       quantiles = my.probs)
era5_wind_intens_quant_trends_1950_2014_MA <- plot_era5_intensification_quantile_trends(era5_intens_rates_wind_1950_2014_MA, 
                                                                                        "max_wind", 
                                                                                        quantiles = my.probs)

# 1979-2014
era5_slp_intens_quant_trends_1979_2014 <- plot_era5_intensification_quantile_trends(era5_intens_rates_slp_1979_2014, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
era5_wind_intens_quant_trends_1979_2014 <- plot_era5_intensification_quantile_trends(era5_intens_rates_wind_1979_2014, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# SON
era5_slp_intens_quant_trends_1979_2014_SON <- plot_era5_intensification_quantile_trends(era5_intens_rates_slp_1979_2014_SON, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
era5_wind_intens_quant_trends_1979_2014_SON <- plot_era5_intensification_quantile_trends(era5_intens_rates_wind_1979_2014_SON, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# DJF
era5_slp_intens_quant_trends_1979_2014_DJF <- plot_era5_intensification_quantile_trends(era5_intens_rates_slp_1979_2014_DJF, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
era5_wind_intens_quant_trends_1979_2014_DJF <- plot_era5_intensification_quantile_trends(era5_intens_rates_wind_1979_2014_DJF, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# MA
era5_slp_intens_quant_trends_1979_2014_MA <- plot_era5_intensification_quantile_trends(era5_intens_rates_slp_1979_2014_MA, 
                                                                                       "min_slp", 
                                                                                       quantiles = my.probs)
era5_wind_intens_quant_trends_1979_2014_MA <- plot_era5_intensification_quantile_trends(era5_intens_rates_wind_1979_2014_MA, 
                                                                                        "max_wind", 
                                                                                        quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-ERA5, 1979-2014

# Full dataset
crcm5_era5_slp_intens_quant_trends_1979_2014 <- plot_era5_intensification_quantile_trends(crcm5_era5_intens_rates_slp_1979_2014, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
crcm5_era5_wind_intens_quant_trends_1979_2014 <- plot_era5_intensification_quantile_trends(crcm5_era5_intens_rates_wind_1979_2014, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# SON
crcm5_era5_slp_intens_quant_trends_1979_2014_SON <- plot_era5_intensification_quantile_trends(crcm5_era5_intens_rates_slp_1979_2014_SON, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
crcm5_era5_wind_intens_quant_trends_1979_2014_SON <- plot_era5_intensification_quantile_trends(crcm5_era5_intens_rates_wind_1979_2014_SON, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# DJF
crcm5_era5_slp_intens_quant_trends_1979_2014_DJF <- plot_era5_intensification_quantile_trends(crcm5_era5_intens_rates_slp_1979_2014_DJF, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
crcm5_era5_wind_intens_quant_trends_1979_2014_DJF <- plot_era5_intensification_quantile_trends(crcm5_era5_intens_rates_wind_1979_2014_DJF, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# MA
crcm5_era5_slp_intens_quant_trends_1979_2014_MA <- plot_era5_intensification_quantile_trends(crcm5_era5_intens_rates_slp_1979_2014_MA, 
                                                                                       "min_slp", 
                                                                                       quantiles = my.probs)
crcm5_era5_wind_intens_quant_trends_1979_2014_MA <- plot_era5_intensification_quantile_trends(crcm5_era5_intens_rates_wind_1979_2014_MA, 
                                                                                        "max_wind", 
                                                                                        quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-CMIP6 ensemble (pooled pilot models)

# 1950-2014
crcm5_cmip6_ensemble_slp_intens_quant_trends_1950_2014 <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1950_2014,  
                                                                                                  quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_1950_2014 <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1950_2014, 
                                                                                                   quantiles = my.probs)
# SON
crcm5_cmip6_ensemble_slp_intens_quant_trends_1950_2014_SON <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1950_2014_SON,  
                                                                                                        quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_1950_2014_SON <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1950_2014_SON, 
                                                                                                         quantiles = my.probs)
# DJF
crcm5_cmip6_ensemble_slp_intens_quant_trends_1950_2014_DJF <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1950_2014_DJF,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_1950_2014_DJF <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1950_2014_DJF, 
                                                                                                             quantiles = my.probs)
# MA
crcm5_cmip6_ensemble_slp_intens_quant_trends_1950_2014_MA <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1950_2014_MA,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_1950_2014_MA <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1950_2014_MA, 
                                                                                                             quantiles = my.probs)

# 1979-2014
crcm5_cmip6_ensemble_slp_intens_quant_trends_1979_2014 <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1979_2014, 
                                                                                                        quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_1979_2014 <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1979_2014, 
                                                                                                         quantiles = my.probs)
# SON
crcm5_cmip6_ensemble_slp_intens_quant_trends_1979_2014_SON <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1979_2014_SON,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_1979_2014_SON <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1979_2014_SON, 
                                                                                                             quantiles = my.probs)
# DJF
crcm5_cmip6_ensemble_slp_intens_quant_trends_1979_2014_DJF <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1979_2014_DJF,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_1979_2014_DJF <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1979_2014_DJF, 
                                                                                                             quantiles = my.probs)
# MA
crcm5_cmip6_ensemble_slp_intens_quant_trends_1979_2014_MA <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1979_2014_MA,  
                                                                                                           quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_1979_2014_MA <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1979_2014_MA, 
                                                                                                            quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-CMIP6 intensification rate quantiles by model

# 1950-2014
crcm5_cmip6_models_slp_intens_quant_trends_1950_2014 <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1950_2014,  
                                                                                                        quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_1950_2014 <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1950_2014, 
                                                                                                         quantiles = my.probs)
# SON
crcm5_cmip6_models_slp_intens_quant_trends_1950_2014_SON <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1950_2014_SON,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_1950_2014_SON <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1950_2014_SON, 
                                                                                                             quantiles = my.probs)
# DJF
crcm5_cmip6_models_slp_intens_quant_trends_1950_2014_DJF <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1950_2014_DJF,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_1950_2014_DJF <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1950_2014_DJF, 
                                                                                                             quantiles = my.probs)
# MA
crcm5_cmip6_models_slp_intens_quant_trends_1950_2014_MA <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1950_2014_MA,  
                                                                                                           quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_1950_2014_MA <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1950_2014_MA, 
                                                                                                            quantiles = my.probs)

# 1979-2014
crcm5_cmip6_models_slp_intens_quant_trends_1979_2014 <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1979_2014, 
                                                                                                        quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_1979_2014 <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1979_2014, 
                                                                                                         quantiles = my.probs)
# SON
crcm5_cmip6_models_slp_intens_quant_trends_1979_2014_SON <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1979_2014_SON,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_1979_2014_SON <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1979_2014_SON, 
                                                                                                             quantiles = my.probs)
# DJF
crcm5_cmip6_models_slp_intens_quant_trends_1979_2014_DJF <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1979_2014_DJF,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_1979_2014_DJF <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1979_2014_DJF, 
                                                                                                             quantiles = my.probs)
# MA
crcm5_cmip6_models_slp_intens_quant_trends_1979_2014_MA <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_slp_1979_2014_MA,  
                                                                                                           quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_1979_2014_MA <- plot_intensification_quantile_trends(crcm5_cmip6_intens_rates_wind_1979_2014_MA, 
                                                                                                            quantiles = my.probs)

# ==============================================================================

# Grid-based storm track counting, 1979-2014

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CFSR

n_lat_cfsr <- nrow(cfsr_lat)
n_lon_cfsr <- ncol(cfsr_lon)
lon_out <- seq(-100, -10, by = 1)
lat_out <- seq(10, 70, by = 1)   # 1 deg. by 1 deg. grid

cat("\n=== Counting storms in grid: 1979-2014 ===\n")
cfsr_track_counts_1979_2014 <- count_storms_in_grid_era5( # same function as ERA5
  cfsr_strong_noreasters_1979_2014,
  n_lat_cfsr, n_lon_cfsr, cfsr_lat, cfsr_lon
)
cfsr_track_counts_1979_2014 <- mat_to_regular_grid_rect(cfsr_track_counts_1979_2014, 
                                                        lon_out, lat_out)
cfsr_track_counts_1979_2014$z <- cfsr_track_counts_1979_2014$z / 34    # normalize by number of storm seasons

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# ERA5

n_lat_era5 <- nrow(era5_lat)
n_lon_era5 <- ncol(era5_lon)
lon_out <- seq(-90, -40, by = 1)
lat_out <- seq(20, 50, by = 1)   # 1 deg. by 1 deg. grid

cat("\n=== Counting storms in grid: 1979-2014 ===\n")
era5_track_counts_1979_2014 <- count_storms_in_grid_era5(
  era5_strong_noreasters_1979_2014,
  n_lat_era5, n_lon_era5, era5_lat, era5_lon
)
era5_track_counts_1979_2014 <- mat_to_regular_grid_rect(era5_track_counts_1979_2014, 
                                                        lon_out, lat_out)
era5_track_counts_1979_2014$z <- era5_track_counts_1979_2014$z / 34    # normalize by number of storm seasons

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-ERA5 & CRCM5-CMIP6 (by model)

n_lat_crcm5 <- nrow(crcm5_lat)  # 628
n_lon_crcm5 <- ncol(crcm5_lon)  # 655
lon_out <- seq(-110, -7, by = 1)
lat_out <- seq(15, 82, by = 1)    # 1 deg. by 1 deg. grid

cat(sprintf("Using CRCM5 grid: %d x %d cells\n", n_lat_crcm5, n_lon_crcm5))

cat("\n=== Counting storms in grid: 1979-2014 ===\n")
crcm5_track_counts_1979_2014 <- count_storms_in_grid(
  strong_crcm5_noreasters_1979_2014,
  c("canesm5", "cnrm", "noresm2", "mpi", "era5"),
  n_lat_crcm5, n_lon_crcm5, crcm5_lat, crcm5_lon
)

for (model_name in c("canesm5", "cnrm", "mpi", "noresm2", "era5")) {
  crcm5_track_counts_1979_2014[[model_name]] <- mat_to_regular_grid(crcm5_track_counts_1979_2014[[model_name]], 
                                                                    lon_out, lat_out)
  # normalize to change in number of storms per year
  crcm5_track_counts_1979_2014[[model_name]]$z <- crcm5_track_counts_1979_2014[[model_name]]$z / 34
}

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# CRCM5-CMIP6 ensemble calculation using individual model track counts

available_models <- gcm_models[gcm_models %in% names(crcm5_track_counts_1979_2014)]
if (length(available_models) > 1) {
  z_dims <- dim(crcm5_track_counts_1979_2014[[available_models[1]]]$z)
  crcm5_cmip6_ensemble_array <- array(NA, dim = c(z_dims[1], z_dims[2], 
                                                  length(available_models)))
  
  for (i in 1:length(available_models)) {
    crcm5_cmip6_ensemble_array[, , i] <- crcm5_track_counts_1979_2014[[available_models[i]]]$z
  }
  
  crcm5_cmip6_ensemble_mean_z_1979_2014 <- apply(crcm5_cmip6_ensemble_array, 
                                                 c(1, 2), mean, na.rm = TRUE)
  crcm5_cmip6_ensemble_counts_1979_2014 <- crcm5_track_counts_1979_2014[[available_models[1]]]
  crcm5_cmip6_ensemble_counts_1979_2014$z <- crcm5_cmip6_ensemble_mean_z_1979_2014
  cat("crcm5_cmip6_ensemble mean 1979-2014 calculated\n")
}

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# Calculate track count color bar breaks for each dataset, 1979-2014

# CFSR
max_cfsr_track_counts_1979_2014 <- max(abs(range(cfsr_track_counts_1979_2014$z)))
colbar_breaks_cfsr_track_counts_1979_2014 <- seq(0, max_cfsr_track_counts_1979_2014, length.out = 301)
# ERA5
max_era5_track_counts_1979_2014 <- max(abs(range(era5_track_counts_1979_2014$z)))
colbar_breaks_era5_track_counts_1979_2014 <- seq(0, max_era5_track_counts_1979_2014, length.out = 301)
# CRCM5-ERA5
max_crcm5_era5_track_counts_1979_2014 <- max(abs(range(crcm5_track_counts_1979_2014$era5$z, na.rm = TRUE)))
colbar_breaks_crcm5_era5_track_counts_1979_2014 <- seq(0, max_crcm5_era5_track_counts_1979_2014, length.out = 301)
# CRCM5-CMIP6 ensemble
max_crcm5_cmip6_track_counts_1979_2014 <- max(abs(range(crcm5_cmip6_ensemble_counts_1979_2014$z, na.rm = TRUE)))
colbar_breaks_crcm5_cmip6_track_counts_1979_2014 <- seq(0, max_crcm5_cmip6_track_counts_1979_2014, length.out = 301)

# ==============================================================================

# Extract maximum intensity locations & calculate kernel density estimate
#     (location of storm center at the time peak intensity occurs)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# CFSR, min. SLP & max. wind

cfsr_min_slp_locations_1979_2014 <- extract_era5_min_slp_locations(cfsr_strong_noreasters_1979_2014, 
                                                                   cfsr_lon, cfsr_lat)
cfsr_max_wind_locations_1979_2014 <- extract_era5_max_wind_locations(cfsr_strong_noreasters_1979_2014, 
                                                                     wind_variable = "max_wind_in_box", 
                                                                     cfsr_lon, cfsr_lat)

cfsr_min_slp_kde_1979_2014 <- kde2d(cfsr_min_slp_locations_1979_2014$longitude, 
                                    cfsr_min_slp_locations_1979_2014$latitude, 
                                    n = 100, lims = c(-100, -10, 10, 70))
cfsr_max_wind_kde_1979_2014 <- kde2d(cfsr_max_wind_locations_1979_2014$longitude, 
                                     cfsr_max_wind_locations_1979_2014$latitude, 
                                     n = 100, lims = c(-100, -10, 10, 70))

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# ERA5, min. SLP & max. wind

era5_min_slp_locations_1979_2014 <- extract_era5_min_slp_locations(era5_strong_noreasters_1979_2014, 
                                                                   era5_lon, era5_lat)
era5_max_wind_locations_1979_2014 <- extract_era5_max_wind_locations(era5_strong_noreasters_1979_2014, 
                                                                     wind_variable = "max_w_h_750 (m/s)", 
                                                                     era5_lon, era5_lat)

era5_min_slp_kde_1979_2014 <- kde2d(era5_min_slp_locations_1979_2014$longitude, 
                                    era5_min_slp_locations_1979_2014$latitude, 
                                    n = 100, lims = c(-90, -40, 20, 50))
era5_max_wind_kde_1979_2014 <- kde2d(era5_max_wind_locations_1979_2014$longitude, 
                                     era5_max_wind_locations_1979_2014$latitude, 
                                     n = 100, lims = c(-90, -40, 20, 50))

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# CRCM5-ERA5, min. SLP & max. wind

# Make CRCM5-ERA5 list (make data compatible with location extraction function)
crcm5_era5_strong_noreasters_1979_2014 <- list(era5 = strong_crcm5_noreasters_1979_2014$era5)

crcm5_era5_min_slp_locations_1979_2014 <- extract_min_slp_locations(crcm5_era5_strong_noreasters_1979_2014, 
                                                                    crcm5_lon, crcm5_lat)
crcm5_era5_max_wind_locations_1979_2014 <- extract_max_wind_locations(crcm5_era5_strong_noreasters_1979_2014, 
                                                                            crcm5_lon, crcm5_lat)

crcm5_era5_min_slp_kde_1979_2014 <- calculate_max_intensity_kde(crcm5_era5_min_slp_locations_1979_2014, 
                                                           kde_resolution = 100, 
                                                           min(crcm5_lon), min(crcm5_lat), 
                                                           max(crcm5_lon), max(crcm5_lat), 
                                                           multi_member_models = c())
crcm5_era5_max_wind_kde_1979_2014 <- calculate_max_intensity_kde(crcm5_era5_max_wind_locations_1979_2014, 
                                                            kde_resolution = 100, 
                                                            min(crcm5_lon), min(crcm5_lat), 
                                                            max(crcm5_lon), max(crcm5_lat), 
                                                            multi_member_models = c())

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# CRCM5-CMIP6, min. SLP & max. wind

crcm5_cmip6_min_slp_locations_1979_2014 <- extract_min_slp_locations(strong_crcm5_noreasters_1979_2014[gcm_models], 
                                                                     crcm5_lon, crcm5_lat)
crcm5_cmip6_max_wind_locations_1979_2014 <- extract_max_wind_locations(strong_crcm5_noreasters_1979_2014[gcm_models], 
                                                                     crcm5_lon, crcm5_lat)

crcm5_cmip6_min_slp_kde_1979_2014 <- calculate_max_intensity_kde(crcm5_cmip6_min_slp_locations_1979_2014, 
                                                                 kde_resolution = 100, 
                                                                 min(crcm5_lon), min(crcm5_lat), 
                                                                 max(crcm5_lon), max(crcm5_lat), 
                                                                 multi_member_models = c())
crcm5_cmip6_max_wind_kde_1979_2014 <- calculate_max_intensity_kde(crcm5_cmip6_max_wind_locations_1979_2014, 
                                                                 kde_resolution = 100, 
                                                                 min(crcm5_lon), min(crcm5_lat), 
                                                                 max(crcm5_lon), max(crcm5_lat), 
                                                                 multi_member_models = c())

# Min. SLP ensemble kde
available_models <- gcm_models[gcm_models %in% names(crcm5_cmip6_min_slp_kde_1979_2014)]
z_dims <- dim(crcm5_cmip6_min_slp_kde_1979_2014[[available_models[1]]]$z)
ensemble_array <- array(NA, dim = c(z_dims[1], z_dims[2], length(available_models)))
for (i in 1:length(available_models)) {
  ensemble_array[, , i] <- crcm5_cmip6_min_slp_kde_1979_2014[[available_models[i]]]$z
}
crcm5_cmip6_ensemble_min_slp_kde_1979_2014 <- crcm5_cmip6_min_slp_kde_1979_2014[[available_models[1]]]
crcm5_cmip6_ensemble_min_slp_kde_1979_2014$z <- apply(ensemble_array, c(1, 2), mean, na.rm = TRUE)

# Max. wind ensemble kde
available_models <- gcm_models[gcm_models %in% names(crcm5_cmip6_max_wind_kde_1979_2014)]
z_dims <- dim(crcm5_cmip6_max_wind_kde_1979_2014[[available_models[1]]]$z)
ensemble_array <- array(NA, dim = c(z_dims[1], z_dims[2], length(available_models)))
for (i in 1:length(available_models)) {
  ensemble_array[, , i] <- crcm5_cmip6_max_wind_kde_1979_2014[[available_models[i]]]$z
}
crcm5_cmip6_ensemble_max_wind_kde_1979_2014 <- crcm5_cmip6_max_wind_kde_1979_2014[[available_models[1]]]
crcm5_cmip6_ensemble_max_wind_kde_1979_2014$z <- apply(ensemble_array, c(1, 2), mean, na.rm = TRUE)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# Calculate min. SLP & max. wind color bar breaks for each dataset, 1979-2014

# CFSR
max_cfsr_min_slp_1979_2014 <- max(abs(range(cfsr_min_slp_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_cfsr_min_slp_1979_2014 <- seq(0, max_cfsr_min_slp_1979_2014, length.out = 301)
max_cfsr_max_wind_1979_2014 <- max(abs(range(cfsr_max_wind_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_cfsr_max_wind_1979_2014 <- seq(0, max_cfsr_max_wind_1979_2014, length.out = 301)

# ERA5
max_era5_min_slp_1979_2014 <- max(abs(range(era5_min_slp_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_era5_min_slp_1979_2014 <- seq(0, max_era5_min_slp_1979_2014, length.out = 301)
max_era5_max_wind_1979_2014 <- max(abs(range(era5_max_wind_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_era5_max_wind_1979_2014 <- seq(0, max_era5_max_wind_1979_2014, length.out = 301)

# CRCM5-ERA5
max_crcm5_era5_min_slp_1979_2014 <- max(abs(range(crcm5_era5_min_slp_kde_1979_2014$era5$z, na.rm = TRUE)))
colbar_breaks_crcm5_era5_min_slp_1979_2014 <- seq(0, max_crcm5_era5_min_slp_1979_2014, length.out = 301)
max_crcm5_era5_max_wind_1979_2014 <- max(abs(range(crcm5_era5_max_wind_kde_1979_2014$era5$z, na.rm = TRUE)))
colbar_breaks_crcm5_era5_max_wind_1979_2014 <- seq(0, max_crcm5_era5_max_wind_1979_2014, length.out = 301)

# CRCM5-CMIP6
max_crcm5_cmip6_ensemble_min_slp_1979_2014 <- max(abs(range(crcm5_cmip6_ensemble_min_slp_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_crcm5_cmip6_min_slp_1979_2014 <- seq(0, max_crcm5_cmip6_ensemble_min_slp_1979_2014, length.out = 301)
max_crcm5_cmip6_ensemble_max_wind_1979_2014 <- max(abs(range(crcm5_cmip6_ensemble_max_wind_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_crcm5_cmip6_max_wind_1979_2014 <- seq(0, max_crcm5_cmip6_ensemble_max_wind_1979_2014, length.out = 301)

# ==============================================================================

# Extract maximum intensification rate locations & calculate kernel densiy estimate
#   (location of storm center at time maximum rate occurs)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# CFSR

cfsr_slp_intens_locations_1979_2014 <- extract_era5_max_intensification_locations(cfsr_intens_rates_slp_1979_2014, 
                                                                                  cfsr_lon, cfsr_lat)
cfsr_wind_intens_locations_1979_2014 <- extract_era5_max_intensification_locations(cfsr_intens_rates_wind_1979_2014, 
                                                                                  cfsr_lon, cfsr_lat)

cfsr_slp_intens_kde_1979_2014 <- kde2d(cfsr_slp_intens_locations_1979_2014$longitude, 
                                       cfsr_slp_intens_locations_1979_2014$latitude, 
                                       n = 100, lims = c(-100, -10, 10, 70))
cfsr_wind_intens_kde_1979_2014 <- kde2d(cfsr_wind_intens_locations_1979_2014$longitude, 
                                       cfsr_wind_intens_locations_1979_2014$latitude, 
                                       n = 100, lims = c(-100, -10, 10, 70))

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# ERA5

era5_slp_intens_locations_1979_2014 <- extract_era5_max_intensification_locations(era5_intens_rates_slp_1979_2014, 
                                                                                  era5_lon, era5_lat)
era5_wind_intens_locations_1979_2014 <- extract_era5_max_intensification_locations(era5_intens_rates_wind_1979_2014, 
                                                                                   era5_lon, era5_lat)

era5_slp_intens_kde_1979_2014 <- kde2d(era5_slp_intens_locations_1979_2014$longitude, 
                                       era5_slp_intens_locations_1979_2014$latitude, 
                                       n = 100, lims = c(-90, -40, 20, 50))
era5_wind_intens_kde_1979_2014 <- kde2d(era5_wind_intens_locations_1979_2014$longitude, 
                                        era5_wind_intens_locations_1979_2014$latitude, 
                                        n = 100, lims = c(-90, -40, 20, 50))

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# CRCM5-ERA5

# Make CRCM5-ERA5 intens. rates list (make compatible with extraction function)
crcm5_era5_intens_rates_slp_1979_2014 <- list(era5 = crcm5_era5_intens_rates_slp_1979_2014)
crcm5_era5_intens_rates_wind_1979_2014 <- list(era5 = crcm5_era5_intens_rates_wind_1979_2014)

crcm5_era5_slp_intens_locations_1979_2014 <- extract_max_intensification_locations(crcm5_era5_intens_rates_slp_1979_2014, 
                                                                                   crcm5_lon, crcm5_lat)
crcm5_era5_wind_intens_locations_1979_2014 <- extract_max_intensification_locations(crcm5_era5_intens_rates_wind_1979_2014, 
                                                                                   crcm5_lon, crcm5_lat)

crcm5_era5_slp_intens_kde_1979_2014 <- calculate_max_intensity_kde(crcm5_era5_slp_intens_locations_1979_2014, 
                                                                   kde_resolution = 100, 
                                                                   min(crcm5_lon), min(crcm5_lat), 
                                                                   max(crcm5_lon), max(crcm5_lat), 
                                                                   multi_member_models = c())
crcm5_era5_wind_intens_kde_1979_2014 <- calculate_max_intensity_kde(crcm5_era5_wind_intens_locations_1979_2014, 
                                                                   kde_resolution = 100, 
                                                                   min(crcm5_lon), min(crcm5_lat), 
                                                                   max(crcm5_lon), max(crcm5_lat), 
                                                                   multi_member_models = c())

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# CRCM5-CMIP6

crcm5_cmip6_slp_intens_locations_1979_2014 <- extract_max_intensification_locations(crcm5_cmip6_intens_rates_slp_1979_2014, 
                                                                                    crcm5_lon, crcm5_lat)
crcm5_cmip6_wind_intens_locations_1979_2014 <- extract_max_intensification_locations(crcm5_cmip6_intens_rates_wind_1979_2014, 
                                                                                    crcm5_lon, crcm5_lat)

crcm5_cmip6_slp_intens_kde_1979_2014 <- calculate_max_intensity_kde(crcm5_cmip6_slp_intens_locations_1979_2014, 
                                                               kde_resolution = 100, 
                                                               min(crcm5_lon), min(crcm5_lat), 
                                                               max(crcm5_lon), max(crcm5_lat), 
                                                               multi_member_models = c())
crcm5_cmip6_wind_intens_kde_1979_2014 <- calculate_max_intensity_kde(crcm5_cmip6_wind_intens_locations_1979_2014, 
                                                               kde_resolution = 100, 
                                                               min(crcm5_lon), min(crcm5_lat), 
                                                               max(crcm5_lon), max(crcm5_lat), 
                                                               multi_member_models = c())

# Min. SLP rate ensemble kde
available_models <- gcm_models[gcm_models %in% names(crcm5_cmip6_slp_intens_kde_1979_2014)]
z_dims <- dim(crcm5_cmip6_slp_intens_kde_1979_2014[[available_models[1]]]$z)
ensemble_array <- array(NA, dim = c(z_dims[1], z_dims[2], length(available_models)))
for (i in 1:length(available_models)) {
  ensemble_array[, , i] <- crcm5_cmip6_slp_intens_kde_1979_2014[[available_models[i]]]$z
}
crcm5_cmip6_ensemble_slp_intens_kde_1979_2014 <- crcm5_cmip6_slp_intens_kde_1979_2014[[available_models[1]]]
crcm5_cmip6_ensemble_slp_intens_kde_1979_2014$z <- apply(ensemble_array, c(1, 2), mean, na.rm = TRUE)

# Max. wind rate ensemble kde
available_models <- gcm_models[gcm_models %in% names(crcm5_cmip6_wind_intens_kde_1979_2014)]
z_dims <- dim(crcm5_cmip6_wind_intens_kde_1979_2014[[available_models[1]]]$z)
ensemble_array <- array(NA, dim = c(z_dims[1], z_dims[2], length(available_models)))
for (i in 1:length(available_models)) {
  ensemble_array[, , i] <- crcm5_cmip6_wind_intens_kde_1979_2014[[available_models[i]]]$z
}
crcm5_cmip6_ensemble_wind_intens_kde_1979_2014 <- crcm5_cmip6_wind_intens_kde_1979_2014[[available_models[1]]]
crcm5_cmip6_ensemble_wind_intens_kde_1979_2014$z <- apply(ensemble_array, c(1, 2), mean, na.rm = TRUE)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
# Calculate min. SLP & max. wind rate color bar breaks for each dataset, 1979-2014

# CFSR
max_cfsr_slp_intens_1979_2014 <- max(abs(range(cfsr_slp_intens_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_cfsr_slp_intens_1979_2014 <- seq(0, max_cfsr_slp_intens_1979_2014, length.out = 301)
max_cfsr_wind_intens_1979_2014 <- max(abs(range(cfsr_wind_intens_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_cfsr_wind_intens_1979_2014 <- seq(0, max_cfsr_wind_intens_1979_2014, length.out = 301)

# ERA5
max_era5_slp_intens_1979_2014 <- max(abs(range(era5_slp_intens_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_era5_slp_intens_1979_2014 <- seq(0, max_era5_slp_intens_1979_2014, length.out = 301)
max_era5_wind_intens_1979_2014 <- max(abs(range(era5_wind_intens_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_era5_wind_intens_1979_2014 <- seq(0, max_era5_wind_intens_1979_2014, length.out = 301)

# CRCM5-ERA5
max_crcm5_era5_slp_intens_1979_2014 <- max(abs(range(crcm5_era5_slp_intens_kde_1979_2014$era5$z, na.rm = TRUE)))
colbar_breaks_crcm5_era5_slp_intens_1979_2014 <- seq(0, max_crcm5_era5_slp_intens_1979_2014, length.out = 301)
max_crcm5_era5_wind_intens_1979_2014 <- max(abs(range(crcm5_era5_wind_intens_kde_1979_2014$era5$z, na.rm = TRUE)))
colbar_breaks_crcm5_era5_wind_intens_1979_2014 <- seq(0, max_crcm5_era5_wind_intens_1979_2014, length.out = 301)

# CRCM5-CMIP6
max_crcm5_cmip6_ensemble_slp_intens_1979_2014 <- max(abs(range(crcm5_cmip6_ensemble_slp_intens_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_crcm5_cmip6_slp_intens_1979_2014 <- seq(0, max_crcm5_cmip6_ensemble_slp_intens_1979_2014, length.out = 301)
max_crcm5_cmip6_ensemble_wind_intens_1979_2014 <- max(abs(range(crcm5_cmip6_ensemble_wind_intens_kde_1979_2014$z, na.rm = TRUE)))
colbar_breaks_crcm5_cmip6_wind_intens_1979_2014 <- seq(0, max_crcm5_cmip6_ensemble_wind_intens_1979_2014, length.out = 301)

# ==============================================================================

# Historical MPI member intensity analysis ONLY: includes r1i1p1f1, r2i1p1f1,
#   r3i1p1f1, r4i1p1f1, and r5i1p1f1

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Load MPI member data

# Define base directory
mpi_base_dir <- "/Path/to/CRCM5/MPI/data/mpi_historical_unfiltered_tracks"
mpi_members <- c("sim1", "sim2", "sim3", "sim4", "sim5")

# Define members
members <- list(
  sim1 = file.path(mpi_base_dir, "sim1"),
  sim2 = file.path(mpi_base_dir, "sim2"),
  sim3 = file.path(mpi_base_dir, "sim3"),
  sim4 = file.path(mpi_base_dir, "sim4"),
  sim5 = file.path(mpi_base_dir, "sim5")
)

# Initialize storage list
mpi_noreasters <- list()

cat("Loading members...\n")
for (model in names(members)) {
  cat(paste("Loading", model, "...\n"))
  
  model_path <- members[[model]]
  csv_files <- list.files(model_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Create a list to store this model's data
  model_data <- list()
  
  for (csv_file in csv_files) {
    year <- extract_year(csv_file)
    
    if (!is.na(year)) {
      tryCatch({
        df <- read_storm_csv(csv_file)
        model_data[[as.character(year)]] <- df
        cat(paste("  ", year, ":", nrow(df), "storms\n"))
      }, error = function(e) {
        cat(paste("  Error reading", csv_file, ":", e$message, "\n"))
      })
    } else {
      cat(paste("  Warning: Could not extract year from", basename(csv_file), "\n"))
    }
  }
  
  mpi_noreasters[[model]] <- model_data
}

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Filter MPI member data: time, distance, TC removal, USEC box, strong storms, 
#   & time periods

mpi_wind_dirs <- list(
  sim1 = file.path("/Path/to/CRCM5/MPI/wind/data/MPI_wind_data/mpi_historical1"),
  sim2 = file.path("/Path/to/CRCM5/MPI/wind/dataMPI_wind_data/mpi_historical2"),
  sim3 = file.path("/Path/to/CRCM5/MPI/wind/data/MPI_wind_data/mpi_historical3"),
  sim4 = file.path("/Path/to/CRCM5/MPI/wind/data/MPI_wind_data/mpi_historical4"),
  sim5 = file.path("/Path/to/CRCM5/MPI/wind/data/MPI_wind_Data/mpi_historical5")
)

mpi_noreasters <- filter_crcm5_noreasters(mpi_noreasters, crcm5_lat, crcm5_lon, 
                                        wind_dirs = mpi_wind_dirs,
                                        apply_region_filter = FALSE, 
                                        apply_TC_filter = TRUE,
                                        TC_min_lat = 30, 
                                        TC_min_slp = 98000, 
                                        TC_max_wind = 20)

# Filter for just strong nor'easters
mpi_strong_noreasters <- apply_SLP_filter(mpi_noreasters, limit = 98000, 
                                          lim_min = FALSE, era5 = FALSE, 
                                          dataframe = TRUE)

# Filter mpi storms to common time period w/ other datasets
mpi_strong_noreasters_1979_2014 <- filter_years_multi_model(mpi_strong_noreasters, 
                                                            start_year = 1980, # keyed by end year
                                                            end_year = 2014)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Extract intensity values, 1979-2014

mpi_min_slp_1979_2014 <- extract_mpi_intensity_values(mpi_strong_noreasters_1979_2014, 
                                                      metric = "min_slp")
mpi_max_wind_1979_2014 <- extract_mpi_intensity_values(mpi_strong_noreasters_1979_2014, 
                                                       metric = "max_wind")

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Calculate pooled intensity quantiles for each MPI member, 1979-2014

# Sim1
mpi_sim1_min_slp_quants_1979_2014 <- calculate_mpi_intensity_quantiles(mpi_min_slp_1979_2014$sim1, 
                                                                       metric = "min_slp", my.probs)
mpi_sim1_max_wind_quants_1979_2014 <- calculate_mpi_intensity_quantiles(mpi_max_wind_1979_2014$sim1, 
                                                                        metric = "max_wind", my.probs)

# Sim2
mpi_sim2_min_slp_quants_1979_2014 <- calculate_mpi_intensity_quantiles(mpi_min_slp_1979_2014$sim2, 
                                                                       metric = "min_slp", my.probs)
mpi_sim2_max_wind_quants_1979_2014 <- calculate_mpi_intensity_quantiles(mpi_max_wind_1979_2014$sim2, 
                                                                        metric = "max_wind", my.probs)

# Sim3
mpi_sim3_min_slp_quants_1979_2014 <- calculate_mpi_intensity_quantiles(mpi_min_slp_1979_2014$sim3, 
                                                                       metric = "min_slp", my.probs)
mpi_sim3_max_wind_quants_1979_2014 <- calculate_mpi_intensity_quantiles(mpi_max_wind_1979_2014$sim3, 
                                                                        metric = "max_wind", my.probs)

# Sim4
mpi_sim4_min_slp_quants_1979_2014 <- calculate_mpi_intensity_quantiles(mpi_min_slp_1979_2014$sim4, 
                                                                       metric = "min_slp", my.probs)
mpi_sim4_max_wind_quants_1979_2014 <- calculate_mpi_intensity_quantiles(mpi_max_wind_1979_2014$sim4, 
                                                                        metric = "max_wind", my.probs)

# Sim5
mpi_sim5_min_slp_quants_1979_2014 <- calculate_mpi_intensity_quantiles(mpi_min_slp_1979_2014$sim5, 
                                                                       metric = "min_slp", my.probs)
mpi_sim5_max_wind_quants_1979_2014 <- calculate_mpi_intensity_quantiles(mpi_max_wind_1979_2014$sim5, 
                                                                        metric = "max_wind", my.probs)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Calculate intensity quantiles & trends over time for each MPI member, 1979-2014

mpi_slp_quant_trends_1979_2014 <- plot_intensity_quantile_trends(mpi_min_slp_1979_2014, 
                                                                     "min_slp", 
                                                                     quantiles = my.probs)
mpi_wind_quant_trends_1979_2014 <- plot_intensity_quantile_trends(mpi_max_wind_1979_2014, 
                                                                      "max_wind", 
                                                                      quantiles = my.probs)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Calculate intensification rates for each MPI member, 1979-2014

mpi_intens_rates_slp_1979_2014 <- calculate_max_intensification_rate(mpi_strong_noreasters_1979_2014, 
                                                                     metric ="min_slp", 
                                                                     time_window_hours = 24)
mpi_intens_rates_wind_1979_2014 <- calculate_max_intensification_rate(mpi_strong_noreasters_1979_2014, 
                                                                     metric ="max_wind", 
                                                                     time_window_hours = 24)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Calculate intensification rate quantiles & trends over timefor each MPI 
#   member, 1979-2014

mpi_slp_intens_quant_trends_1979_2014 <- plot_intensification_quantile_trends(mpi_intens_rates_slp_1979_2014, 
                                                                              "min_slp", 
                                                                              quantiles = my.probs)
mpi_wind_intens_quant_trends_1979_2014 <- plot_intensification_quantile_trends(mpi_intens_rates_wind_1979_2014, 
                                                                               "max_wind", 
                                                                               quantiles = my.probs)



















