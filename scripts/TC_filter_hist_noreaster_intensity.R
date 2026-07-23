# ==============================================================================

# Analysis of historical nor'easter intensity characteristics in CFSR, ERA5, 
#   CRCM5, and CRCM5-CMIP6: Tropical cyclone filter development & testing

# ==============================================================================

# Associated publication: Weaver, M.M. et al. (2026). Historical nor'easter
#                         intensity trends in reanalysis and regional climate
#                         simulations: A cross-dataset perspective.

# This file contains all code necessary to recreate Table 1 & Fig. 1, S3-6, & S10
#   in the aforementioned publication. Run this script independently of
#   hist_noreaster_intensity_analysis.R.

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
library(mapdata)

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
#   travel distance, **KEEP** tropical cyclones, U.S. east coast proximity, 
#   strong storms only

# Define wind_1979_2014 data directory (needed for tropical cyclone filtering)
cfsr_wind_dir <- file.path(base_dir, "Path/to/CFSR/wind/data/CFSR_wind_data")

cfsr_noreasters_wTC <- filter_cfsr_noreasters(cfsr_noreasters, cfsr_lat, cfsr_lon, 
                                          wind_dirs = cfsr_wind_dir,
                                          apply_region_filter = FALSE, 
                                          apply_TC_filter = FALSE,
                                          TC_min_lat = 30, TC_min_slp = 98000, 
                                          TC_max_wind = 20)

# Filter for just strong nor'easters (lifetime min. sea-level pressure ≤ 98000 Pa)
cfsr_strong_noreasters_wTC <- apply_SLP_filter(cfsr_noreasters_wTC, limit = 98000, 
                                           lim_min = FALSE, era5 = TRUE, dataframe = TRUE)
                                        # era5 = TRUE indicates non-nested data structure

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Apply post-processing filters to CFSR nor'easters: 24 hr lifetime, 1000 km 
#   travel distance, **REMOVE** tropical cyclones, U.S. east coast proximity, 
#   strong storms only

cfsr_noreasters_noTC <- filter_cfsr_noreasters(cfsr_noreasters, cfsr_lat, cfsr_lon, 
                                                  wind_dirs = cfsr_wind_dir,
                                                  apply_region_filter = FALSE, 
                                                  apply_TC_filter = TRUE,
                                                  TC_min_lat = 30, TC_min_slp = 98000, 
                                                  TC_max_wind = 20)

# Filter for just strong nor'easters (lifetime min. sea-level pressure ≤ 98000 Pa)
cfsr_strong_noreasters_noTC <- apply_SLP_filter(cfsr_noreasters_noTC, limit = 98000, 
                                               lim_min = FALSE, era5 = TRUE, dataframe = TRUE)
                                            # era5 = TRUE indicates non-nested data structure

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
# Apply post-processing to ERA5 nor'easters: **KEEP** tropical cyclones, strong
#   storms only
#   ***since ERA5 nor'easters are from an existing tracked cyclone dataset
#      (Chen et al., 2025), the 24 hr lifetime, 1000 km travel distance, and
#      U.S. east coast post-processing filters have already been applied***

# ERA5 nor'easter dataset pre-TC filtering:
era5_noreasters_wTC <- era5_noreasters

# Filter for just strong nor'easters (lifetime min. sea-level pressure ≤ 98000 Pa)
#   (tropical cyclones still included)
era5_strong_noreasters_wTC <- apply_SLP_filter(era5_noreasters, limit = 98000, 
                                           lim_min = FALSE, era5 = TRUE, dataframe = TRUE)


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Apply post-processing to ERA5 nor'easters: **REMOVE** tropical cyclones, strong
#   storms only

# Apply tropical cyclone filter to full ERA5 nor'easter dataset (no era5_wind_dir directory 
#   needed, ERA5 wind data is integrated into storm track data (Chen et al., 2025))
era5_noreasters_noTC <- TC_filter_ERA5(era5_noreasters, era5_lon, era5_lat, 
                                              min_lat = 30, min_slp = 98000, max_wind = 20)

# Apply tropical cyclone filter to strong ERA5 nor'easters (no era5_wind_dir directory 
#   needed, ERA5 wind data is integrated into storm track data (Chen et al., 2025))
era5_strong_noreasters_noTC <- TC_filter_ERA5(era5_strong_noreasters_wTC, era5_lon, era5_lat, 
                                         min_lat = 30, min_slp = 98000, max_wind = 20)

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
#   travel distance, **KEEP** tropical cyclones, U.S. east coast proximity, 
#   intensity subsets, common time period(s), seasons

# Define wind_1979_2014 data directories
crcm5_wind_dirs <- list(
  canesm5 = file.path(base_dir, "path/to/wind/data/CanESM5_wind_data/canesm5_historical"),
  cnrm = file.path(base_dir, "path/to/wind/data/CNRM_wind_data/cnrm_historical"),
  noresm2 = file.path(base_dir, "path/to/wind/data/NorESM2_wind_data/noresm2_historical"),
  mpi = file.path(base_dir, "path/to/wind/data/MPI_wind_data/mpi_historical1"),
  era5 = file.path(base_dir, "path/to/wind/data/CRCM5-ERA5_wind_data/era5_v1-r1")
)

# Apply main post-processing function (USEC, time, distance, hurricanes, genesis, NO target region)
crcm5_noreasters_wTC <- filter_crcm5_noreasters(crcm5_noreasters, crcm5_lat, crcm5_lon, 
                                                wind_dirs = crcm5_wind_dirs,
                                                apply_region_filter = FALSE, 
                                                apply_TC_filter = FALSE,
                                                TC_min_lat = 30, TC_min_slp = 98000, 
                                                TC_max_wind = 20)

# Filter for just strong nor'easters (lifetime min. sea-level pressure ≤ 98000 Pa)
strong_crcm5_noreasters_wTC <- apply_SLP_filter(crcm5_noreasters_wTC, 
                                                limit = 98000, lim_min = FALSE, 
                                                era5 = FALSE, dataframe = TRUE)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Apply post-processing filters to CRCM5 nor'easters: 24 hr lifetime, 1000 km 
#   travel distance, **KEEP** tropical cyclones, U.S. east coast proximity, 
#   intensity subsets, common time period(s), seasons

# Apply main post-processing function (USEC, time, distance, hurricanes, genesis, NO target region)
crcm5_noreasters_noTC <- filter_crcm5_noreasters(crcm5_noreasters, crcm5_lat, crcm5_lon, 
                                                 wind_dirs = crcm5_wind_dirs,
                                                 apply_region_filter = FALSE, 
                                                 apply_TC_filter = TRUE,
                                                 TC_min_lat = 30, TC_min_slp = 98000, 
                                                 TC_max_wind = 20)

# Filter for just strong nor'easters (lifetime min. sea-level pressure ≤ 98000 Pa)
strong_crcm5_noreasters_noTC <- apply_SLP_filter(crcm5_noreasters_noTC, 
                                                 limit = 98000, lim_min = FALSE, 
                                                 era5 = FALSE, dataframe = TRUE)

# ==============================================================================

# Minimum sea-level pressure data extraction & quantile trend calculation

my.probs = seq(from = 0, to = 1, by = 0.01)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Extract min. SLP values from all datasets, w/ & w/o TCs

# CFSR
cfsr_min_slp_wTC <- process_era5_for_intensity(cfsr_strong_noreasters_wTC, 
                                               metric = "min_slp")
cfsr_min_slp_noTC <- process_era5_for_intensity(cfsr_strong_noreasters_noTC, 
                                                metric = "min_slp")

# ERA5
era5_min_slp_wTC <- process_era5_for_intensity(era5_strong_noreasters_wTC, 
                                               metric = "min_slp")
era5_min_slp_noTC <- process_era5_for_intensity(era5_strong_noreasters_noTC, 
                                                metric = "min_slp")

# CRCM5-ERA5
crcm5_era5_min_slp_wTC <- extract_intensity_values(strong_crcm5_noreasters_wTC["era5"], 
                                                   metric = "min_slp")
crcm5_era5_min_slp_wTC <- crcm5_era5_min_slp_wTC$era5
crcm5_era5_min_slp_noTC <- extract_intensity_values(strong_crcm5_noreasters_noTC["era5"], 
                                                    metric = "min_slp")
crcm5_era5_min_slp_noTC <- crcm5_era5_min_slp_noTC$era5

# CRCM5
crcm5_cmip6_min_slp_wTC <- extract_intensity_values(strong_crcm5_noreasters_wTC[gcm_models], 
                                                    metric = "min_slp")
crcm5_cmip6_min_slp_noTC <- extract_intensity_values(strong_crcm5_noreasters_noTC[gcm_models], 
                                                     metric = "min_slp")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Calculate min. SLP quantiles for all datasets, full time period available for
#   each dataset

# CFSR
cfsr_slp_quant_trends_wTC <- plot_era5_intensity_quantile_trends(cfsr_min_slp_wTC, 
                                                                 metric = "min_slp", 
                                                                 quantiles = my.probs)
cfsr_slp_quant_trends_noTC <- plot_era5_intensity_quantile_trends(cfsr_min_slp_noTC, 
                                                                  metric = "min_slp", 
                                                                  quantiles = my.probs)

# ERA5
era5_slp_quant_trends_wTC <- plot_era5_intensity_quantile_trends(era5_min_slp_wTC, 
                                                                 metric = "min_slp", 
                                                                 quantiles = my.probs)
era5_slp_quant_trends_noTC <- plot_era5_intensity_quantile_trends(era5_min_slp_noTC, 
                                                                  metric = "min_slp", 
                                                                  quantiles = my.probs)

# CRCM5-ERA5
crcm5_era5_slp_quant_trends_wTC <- plot_era5_intensity_quantile_trends(crcm5_era5_min_slp_wTC, 
                                                                       metric = "min_slp", 
                                                                       quantiles = my.probs)
crcm5_era5_slp_quant_trends_noTC <- plot_era5_intensity_quantile_trends(crcm5_era5_min_slp_noTC, 
                                                                        metric = "min_slp", 
                                                                        quantiles = my.probs)

# CRCM5
crcm5_cmip6_slp_quant_trends_wTC <- plot_ensemble_quantile_trends(crcm5_cmip6_min_slp_wTC, 
                                                                  metric = "min_slp", 
                                                                  quantiles = my.probs)
crcm5_cmip6_slp_quant_trends_noTC <- plot_ensemble_quantile_trends(crcm5_cmip6_min_slp_noTC, 
                                                                   metric = "min_slp", 
                                                                   quantiles = my.probs)

# ==============================================================================

# Maximum surface wind speed data extraction & quantile trend calculation

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Extract max. wind values from all datasets, w/ & w/o TCs

# CFSR
cfsr_max_wind_wTC <- process_era5_for_intensity(cfsr_strong_noreasters_wTC, 
                                                metric = "max_wind")
cfsr_max_wind_noTC <- process_era5_for_intensity(cfsr_strong_noreasters_noTC, 
                                                 metric = "max_wind")

# ERA5
era5_max_wind_wTC <- process_era5_for_intensity(era5_strong_noreasters_wTC, 
                                                metric = "max_wind")
era5_max_wind_noTC <- process_era5_for_intensity(era5_strong_noreasters_noTC, 
                                                 metric = "max_wind")

# CRCM5-ERA5
crcm5_era5_max_wind_wTC <- extract_intensity_values(strong_crcm5_noreasters_wTC["era5"], 
                                                    metric = "max_wind")
crcm5_era5_max_wind_wTC <- crcm5_era5_max_wind_wTC$era5
crcm5_era5_max_wind_noTC <- extract_intensity_values(strong_crcm5_noreasters_noTC["era5"], 
                                                     metric = "max_wind")
crcm5_era5_max_wind_noTC <- crcm5_era5_max_wind_noTC$era5

# CRCM5
crcm5_cmip6_max_wind_wTC <- extract_intensity_values(strong_crcm5_noreasters_wTC[gcm_models], 
                                                     metric = "max_wind")
crcm5_cmip6_max_wind_noTC <- extract_intensity_values(strong_crcm5_noreasters_noTC[gcm_models], 
                                                      metric = "max_wind")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Calculate max. wind quantiles for all datasets, full time period available for
#   each dataset

# CFSR
cfsr_wind_quant_trends_wTC <- plot_era5_intensity_quantile_trends(cfsr_max_wind_wTC, 
                                                                  metric = "max_wind", 
                                                                  quantiles = my.probs)
cfsr_wind_quant_trends_noTC <- plot_era5_intensity_quantile_trends(cfsr_max_wind_noTC, 
                                                                   metric = "max_wind", 
                                                                   quantiles = my.probs)

# ERA5
era5_wind_quant_trends_wTC <- plot_era5_intensity_quantile_trends(era5_max_wind_wTC, 
                                                                  metric = "max_wind", 
                                                                  quantiles = my.probs)
era5_wind_quant_trends_noTC <- plot_era5_intensity_quantile_trends(era5_max_wind_noTC, 
                                                                   metric = "max_wind", 
                                                                   quantiles = my.probs)

# CRCM5-ERA5
crcm5_era5_wind_quant_trends_wTC <- plot_era5_intensity_quantile_trends(crcm5_era5_max_wind_wTC, 
                                                                        metric = "max_wind", 
                                                                        quantiles = my.probs)
crcm5_era5_wind_quant_trends_noTC <- plot_era5_intensity_quantile_trends(crcm5_era5_max_wind_noTC, 
                                                                         metric = "max_wind", 
                                                                         quantiles = my.probs)

# CRCM5
crcm5_cmip6_wind_quant_trends_wTC <- plot_ensemble_quantile_trends(crcm5_cmip6_max_wind_wTC, 
                                                                   metric = "max_wind", 
                                                                   quantiles = my.probs)
crcm5_cmip6_wind_quant_trends_noTC <- plot_ensemble_quantile_trends(crcm5_cmip6_max_wind_noTC, 
                                                                    metric = "max_wind", 
                                                                    quantiles = my.probs)

# ==============================================================================

# Minimum sea-level pressure intensification rate & quantile trend calculation

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Calculate min. SLP rates for all datasets, w/ & w/o TCs

# CFSR
cfsr_slp_intens_rates_wTC <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_wTC, 
                                                                      metric ="min_slp", 
                                                                      time_window_hours = 24, 
                                                                      timestep_hours = 6)
cfsr_slp_intens_rates_noTC <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_noTC, 
                                                                       metric ="min_slp", 
                                                                       time_window_hours = 24, 
                                                                       timestep_hours = 6)

# ERA5
era5_slp_intens_rates_wTC <- calculate_era5_max_intensification_rate(era5_strong_noreasters_wTC, 
                                                                      metric ="min_slp", 
                                                                      time_window_hours = 24)
era5_slp_intens_rates_noTC <- calculate_era5_max_intensification_rate(era5_strong_noreasters_noTC, 
                                                                       metric ="min_slp", 
                                                                       time_window_hours = 24)

# CRCM5-ERA5
crcm5_era5_slp_intens_rates_wTC <- calculate_max_intensification_rate(strong_crcm5_noreasters_wTC["era5"], 
                                                                       metric ="min_slp", 
                                                                       time_window_hours = 24)
crcm5_era5_slp_intens_rates_wTC <- crcm5_era5_slp_intens_rates_wTC$era5
crcm5_era5_slp_intens_rates_noTC <- calculate_max_intensification_rate(strong_crcm5_noreasters_noTC["era5"], 
                                                                        metric ="min_slp", 
                                                                        time_window_hours = 24)
crcm5_era5_slp_intens_rates_noTC <- crcm5_era5_slp_intens_rates_noTC$era5

# CRCM5
crcm5_cmip6_slp_intens_rates_wTC <- calculate_max_intensification_rate(strong_crcm5_noreasters_wTC[gcm_models], 
                                                                        metric ="min_slp", 
                                                                        time_window_hours = 24)
crcm5_cmip6_slp_intens_rates_noTC <- calculate_max_intensification_rate(strong_crcm5_noreasters_noTC[gcm_models], 
                                                                         metric ="min_slp", 
                                                                         time_window_hours = 24)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Calculate min. SLP rate quantiles for all datasets, full time period available
#   for each dataset

# CFSR
cfsr_slp_intens_quant_trends_wTC <- plot_era5_intensification_quantile_trends(cfsr_slp_intens_rates_wTC, 
                                                                       metric = "min_slp", 
                                                                       quantiles = my.probs)
cfsr_slp_intens_quant_trends_noTC <- plot_era5_intensification_quantile_trends(cfsr_slp_intens_rates_noTC, 
                                                                        metric = "min_slp", 
                                                                        quantiles = my.probs)

# ERA5
era5_slp_intens_quant_trends_wTC <- plot_era5_intensification_quantile_trends(era5_slp_intens_rates_wTC, 
                                                                       metric = "min_slp", 
                                                                       quantiles = my.probs)
era5_slp_intens_quant_trends_noTC <- plot_era5_intensification_quantile_trends(era5_slp_intens_rates_noTC, 
                                                                        metric = "min_slp", 
                                                                        quantiles = my.probs)

# CRCM5-ERA5
crcm5_era5_slp_intens_quant_trends_wTC <- plot_era5_intensification_quantile_trends(crcm5_era5_slp_intens_rates_wTC, 
                                                                             metric = "min_slp", 
                                                                             quantiles = my.probs)
crcm5_era5_slp_intens_quant_trends_noTC <- plot_era5_intensification_quantile_trends(crcm5_era5_slp_intens_rates_noTC, 
                                                                              metric = "min_slp", 
                                                                              quantiles = my.probs)

# CRCM5
crcm5_cmip6_slp_intens_quant_trends_wTC <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_slp_intens_rates_wTC, 
                                                                                  quantiles = my.probs)
crcm5_cmip6_slp_intens_quant_trends_noTC <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_slp_intens_rates_noTC, 
                                                                                   quantiles = my.probs)

# ==============================================================================

# Maximum surface wind speed intensification rate & quantile trend calculation

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Calculate max. wind rates for all datasets, w/ & w/o TCs

# CFSR
cfsr_wind_intens_rates_wTC <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_wTC, 
                                                                      metric ="max_wind", 
                                                                      time_window_hours = 24, 
                                                                      timestep_hours = 6)
cfsr_wind_intens_rates_noTC <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_noTC, 
                                                                       metric ="max_wind", 
                                                                       time_window_hours = 24, 
                                                                       timestep_hours = 6)

# ERA5
era5_wind_intens_rates_wTC <- calculate_era5_max_intensification_rate(era5_strong_noreasters_wTC, 
                                                                      metric ="max_wind", 
                                                                      time_window_hours = 24)
era5_wind_intens_rates_noTC <- calculate_era5_max_intensification_rate(era5_strong_noreasters_noTC, 
                                                                       metric ="max_wind", 
                                                                       time_window_hours = 24)

# CRCM5-ERA5
crcm5_era5_wind_intens_rates_wTC <- calculate_max_intensification_rate(strong_crcm5_noreasters_wTC["era5"], 
                                                                       metric ="max_wind", 
                                                                       time_window_hours = 24)
crcm5_era5_wind_intens_rates_wTC <- crcm5_era5_wind_intens_rates_wTC$era5
crcm5_era5_wind_intens_rates_noTC <- calculate_max_intensification_rate(strong_crcm5_noreasters_noTC["era5"], 
                                                                        metric ="max_wind", 
                                                                        time_window_hours = 24)
crcm5_era5_wind_intens_rates_noTC <- crcm5_era5_wind_intens_rates_noTC$era5

# CRCM5
crcm5_cmip6_wind_intens_rates_wTC <- calculate_max_intensification_rate(strong_crcm5_noreasters_wTC[gcm_models], 
                                                                        metric ="max_wind", 
                                                                        time_window_hours = 24)
crcm5_cmip6_wind_intens_rates_noTC <- calculate_max_intensification_rate(strong_crcm5_noreasters_noTC[gcm_models], 
                                                                         metric ="max_wind", 
                                                                         time_window_hours = 24)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Calculate max. wind rate quantiles for all datasets, full time period available
#   for each dataset

# CFSR
cfsr_wind_intens_quant_trends_wTC <- plot_era5_intensification_quantile_trends(cfsr_wind_intens_rates_wTC, 
                                                                              metric = "max_wind", 
                                                                              quantiles = my.probs)
cfsr_wind_intens_quant_trends_noTC <- plot_era5_intensification_quantile_trends(cfsr_wind_intens_rates_noTC, 
                                                                               metric = "max_wind", 
                                                                               quantiles = my.probs)

# ERA5
era5_wind_intens_quant_trends_wTC <- plot_era5_intensification_quantile_trends(era5_wind_intens_rates_wTC, 
                                                                              metric = "max_wind", 
                                                                              quantiles = my.probs)
era5_wind_intens_quant_trends_noTC <- plot_era5_intensification_quantile_trends(era5_wind_intens_rates_noTC, 
                                                                               metric = "max_wind", 
                                                                               quantiles = my.probs)

# CRCM5-ERA5
crcm5_era5_wind_intens_quant_trends_wTC <- plot_era5_intensification_quantile_trends(crcm5_era5_wind_intens_rates_wTC, 
                                                                                    metric = "max_wind", 
                                                                                    quantiles = my.probs)
crcm5_era5_wind_intens_quant_trends_noTC <- plot_era5_intensification_quantile_trends(crcm5_era5_wind_intens_rates_noTC, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)

# CRCM5
crcm5_cmip6_wind_intens_quant_trends_wTC <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_wind_intens_rates_wTC, 
                                                                                         quantiles = my.probs)
crcm5_cmip6_wind_intens_quant_trends_noTC <- plot_ensemble_intensification_quantile_trends(crcm5_cmip6_wind_intens_rates_noTC, 
                                                                                          quantiles = my.probs)

# ==============================================================================

# FIGURE 1

# Percentile vs. p-values intensity fig. with & w/o TCs
#   - min. SLP in top row, max. wind in bottom row
#   - each column is a different dataset

pdf("/Path/for/output/figures/fig_1_TC_intensity_p_values.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 7, 2))

# CFSR min. SLP
plot(rev(my.probs), cfsr_slp_quant_trends_wTC$p_values, type = "l", col = ("#E69F00"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), cfsr_slp_quant_trends_noTC$p_values, lwd = 3, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.15), xpd = "True", legend = c("with TCs", "without TCs"), 
       col = "#E69F00", lty = c(2, 1), lwd = 3, cex = 2.5, ncol = 2, bty = "n")
mtext("Min. SLP", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.75, family = 'Helvetica', font = 2, cex = 2)
mtext("CFSR, 1979-2025", 3, 4, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "a", cex = 3.5, font = 4, col = "black", 
     family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# ERA5 min. SLP
plot(rev(my.probs), era5_slp_quant_trends_wTC$p_values, type = "l", col = ("#00A9A5"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_noTC$p_values, lwd = 3, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.15), xpd = "True", legend = c("with TCs", "without TCs"), 
       col = "#00A9A5", lty = c(2, 1), lwd = 3, cex = 2.5, ncol = 2, bty = "n")
mtext("ERA5, 1940-2024", 3, 4, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "b", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CRCM5-ERA5 min. SLP
plot(rev(my.probs), crcm5_era5_slp_quant_trends_wTC$p_values, type = "l", col = ("#CC79A7"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), crcm5_era5_slp_quant_trends_noTC$p_values, lwd = 3, col = "#CC79A7")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.15), xpd = "True", legend = c("with TCs", "without TCs"), 
       col = "#CC79A7", lty = c(2, 1), lwd = 3, cex = 2.5, ncol = 2, bty = "n")
mtext("CRCM5-ERA5, 1979-2020", 3, 4, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "c", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CRCM5-CMIP6 min. SLP
plot(rev(my.probs), crcm5_cmip6_slp_quant_trends_wTC$p_values, type = "l", col = ("#882255"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), crcm5_cmip6_slp_quant_trends_noTC$p_values, lwd = 3, col = "#882255")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.15), xpd = "True", legend = c("with TCs", "without TCs"), 
       col = "#882255", lty = c(2, 1), lwd = 3, cex = 2.5, ncol = 2, bty = "n")
mtext("CRCM5-CMIP6, 1950-2014", 3, 4, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "d", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CFSR max. wind
plot(my.probs, cfsr_wind_quant_trends_wTC$p_values, type = "l", col = ("#E69F00"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, cfsr_wind_quant_trends_noTC$p_values, lwd = 3, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.75, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "e", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# ERA5 max. wind
plot(my.probs, era5_wind_quant_trends_wTC$p_values, type = "l", col = ("#00A9A5"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_noTC$p_values, lwd = 3, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "f", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CRCM5-ERA5 max. wind
plot(my.probs, crcm5_era5_wind_quant_trends_wTC$p_values, type = "l", col = ("#CC79A7"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_era5_wind_quant_trends_noTC$p_values, lwd = 3, col = "#CC79A7")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "g", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CRCM5-CMIP6 max. wind
plot(my.probs, crcm5_cmip6_wind_quant_trends_wTC$p_values, type = "l", col = ("#882255"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_cmip6_wind_quant_trends_noTC$p_values, lwd = 3, col = "#882255")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "h", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

dev.off()

# ==============================================================================

# FIGURE S10

# CRCM5-CMIP6. plot quantile vs. p-values intensification fig.
#   - min. SLP in top row, max. wind in bottom row
#   - each column is a different dataset

pdf("/Path/for/output/figures/fig_s10_TC_intensification_p_values.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 7, 2))

# CFSR min. SLP
plot(my.probs, cfsr_slp_intens_quant_trends_wTC$p_values, type = "l", 
     col = ("#E69F00"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, cfsr_slp_intens_quant_trends_noTC$p_values, lwd = 3, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.15), xpd = "True", legend = c("with TCs", "without TCs"), 
       col = "#E69F00", lty = c(2, 1), lwd = 3, cex = 2.5, ncol = 2, bty = "n")
mtext("Min. SLP Rate", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.75, family = 'Helvetica', font = 2, cex = 2)
mtext("CFSR, 1979-2025", 3, 4, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "a", cex = 3.25, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# ERA5 min. SLP
plot(my.probs, era5_slp_intens_quant_trends_wTC$p_values, type = "l", 
     col = ("#00A9A5"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_noTC$p_values, lwd = 3, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.15), xpd = "True", legend = c("with TCs", "without TCs"), 
       col = "#00A9A5", lty = c(2, 1), lwd = 3, cex = 2.5, ncol = 2, bty = "n")
mtext("ERA5, 1940-2024", 3, 4, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "b", cex = 3.25, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CRCM5-ERA5 min. SLP
plot(my.probs, crcm5_era5_slp_intens_quant_trends_wTC$p_values, type = "l", 
     col = ("#CC79A7"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_era5_slp_intens_quant_trends_noTC$p_values, lwd = 3, col = "#CC79A7")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.15), xpd = "True", legend = c("with TCs", "without TCs"), 
       col = "#CC79A7", lty = c(2, 1), lwd = 3, cex = 2.5, ncol = 2, bty = "n")
mtext("CRCM5-ERA5, 1979-2020", 3, 4, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "c", cex = 3.25, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CRCM5-CMIP6 min. SLP
plot(my.probs, crcm5_cmip6_slp_intens_quant_trends_wTC$p_values, type = "l", 
     col = ("#882255"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_cmip6_slp_intens_quant_trends_noTC$p_values, lwd = 3, col = "#882255")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.15), xpd = "True", legend = c("with TCs", "without TCs"), 
       col = "#882255", lty = c(2, 1), lwd = 3, cex = 2.5, ncol = 2, bty = "n")
mtext("CRCM5-CMIP6, 1950-2014", 3, 4, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "d", cex = 3.25, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CFSR max. wind
plot(my.probs, cfsr_wind_intens_quant_trends_wTC$p_values, type = "l", 
     col = ("#E69F00"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, cfsr_wind_intens_quant_trends_noTC$p_values, lwd = 3, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind Rate", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.75, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "e", cex = 3.25, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# ERA5 max. wind
plot(my.probs, era5_wind_intens_quant_trends_wTC$p_values, type = "l", 
     col = ("#00A9A5"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_noTC$p_values, lwd = 3, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.04, labels = "f", cex = 3.25, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CRCM5-ERA5 max. wind
plot(my.probs, crcm5_era5_wind_intens_quant_trends_wTC$p_values, type = "l", 
     col = ("#CC79A7"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_era5_wind_intens_quant_trends_noTC$p_values, lwd = 3, col = "#CC79A7")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "g", cex = 3.25, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CRCM5-CMIP6 max. wind
plot(my.probs, crcm5_cmip6_wind_intens_quant_trends_wTC$p_values, type = "l", 
     col = ("#882255"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_cmip6_wind_intens_quant_trends_noTC$p_values, lwd = 3, col = "#882255")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "h", cex = 3.25, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

dev.off()

# ==============================================================================

# TC filter testing: using known historical cyclones tracked in CFSR & ERA5

# The following code tests the TC filter functionality by determining whether
#   key autumn nor'easters and tropical cyclones were retained/discarded 
#   correctly in the ERA5 cyclone dataset

# The filter is only tested on & applied to cyclones that occur in September,
#   October, or November (the overlap between Atlantic hurricane season and the 
#   typical nor'easter season) to avoid filtering exceptionally strong nor'easters
#   that occur in winter and spring out of the tracked cyclone datasets

# ==============================================================================

# CFSR NOR'EASTERS

# Nor'easter IDs manually identified by matching dates, genesis locations, 
#   storm lifetimes, and maximum intensity values between tracked cyclone 
#   dataset and documented characteristics for each storm
cfsr_famous_noreasters <- list(
# "Great Appalachian Storm (1950)" N/A since CFSR data starts in 1979
  "Perfect Storm (1991)" = "10291991-06-92-54",
  "Nor'Ida (2009)" = "11122009-00-44-76",
  "Halloween nor'easter (2011)" = "10292011-12-51-68",
  "November 2012 nor'easter (2012)" = "11072012-00-54-72"
)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Search for famous nor'easters in CFSR cyclone dataset both BEFORE and AFTER
#   the TC filter is applied

# Search for each storm in cyclone dataset BEFORE TC filter 
#   (with potential TCs included)
print("--- BEFORE TC FILTER (cfsr_noreasters_wTC) ---")
for (storm_name in names(cfsr_famous_noreasters)) {
  storm_id <- cfsr_famous_noreasters[[storm_name]]
  result <- find_storm_by_id(cfsr_noreasters_wTC, storm_id)
  
  print(paste("\nSearching for:", storm_name, "-", storm_id))
  if (!is.null(result)) {
    print(paste("✓ Found in year:", result$year))
    # ^^^ cyclones are keyed by end year of the nor'easter season (i.e., a fall 1991 storm is keyed to the 1992 season)
    print(paste("  Min pressure:", result$storm_data$`min_p (Pa)`, "Pa"))
    print(paste("  Max wind:", max(unlist(result$storm_data$max_wind_in_box), 
                                   na.rm = TRUE), "m/s"))
    
    # Extract genesis position
    cx_p_vals <- extract_coordinates_robust(result$storm_data$cx_p)
    cy_p_vals <- extract_coordinates_robust(result$storm_data$cy_p)
    cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
    cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
    genesis_coords <- lat_lon(cx_p_first, cy_p_first, cfsr_lon, cfsr_lat)
    
    print(paste("  Genesis latitude:", round(genesis_coords[2], 2), "°N"))
    print(paste("  Genesis longitude:", round(genesis_coords[1], 2), "°E"))
    
  } else {
    print("✗ NOT FOUND in pre-filtered data")
  }
}

# Search for each storm in cyclone dataset AFTER TC filter  (TCs removed)
print("\n--- AFTER TC FILTER (cfsr_noreasters_noTC) ---")
for (storm_name in names(cfsr_famous_noreasters)) {
  storm_id <- cfsr_famous_noreasters[[storm_name]]
  result <- find_storm_by_id(cfsr_noreasters_noTC, storm_id)
  
  print(paste("\nSearching for:", storm_name, "-", storm_id))
  if (!is.null(result)) {
    print(paste("✓ RETAINED - Found in year:", result$year))
    print(paste("  Min pressure:", result$storm_data$`min_p (Pa)`, "Pa"))
    print(paste("  Max wind:", max(unlist(result$storm_data$max_wind_in_box), na.rm = TRUE), "m/s"))
    
    # Extract genesis position
    cx_p_vals <- extract_coordinates_robust(result$storm_data$cx_p)
    cy_p_vals <- extract_coordinates_robust(result$storm_data$cy_p)
    cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
    cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
    genesis_coords <- lat_lon(cx_p_first, cy_p_first, cfsr_lon, cfsr_lat)
    
    print(paste("  Genesis latitude:", round(genesis_coords[2], 2), "°N"))
    print(paste("  Genesis longitude:", round(genesis_coords[1], 2), "°E"))
    
    # Check against TC criteria
    print("\n  TC Filter Criteria Check:")
    print(paste("    Genesis south of 30°N?", genesis_coords[2] < 30, 
                ifelse(genesis_coords[2] < 30, "(would flag as TC)", "(OK)")))
    print(paste("    Min pressure <= 98000 Pa?", result$storm_data$`min_p (Pa)` <= 98000,
                ifelse(result$storm_data$`min_p (Pa)` <= 98000, "(would flag as TC)", "(OK)")))
    print(paste("    Max wind >= 20 m/s?", max(unlist(result$storm_data$max_wind_in_box), na.rm = TRUE) >= 20,
                ifelse(max(unlist(result$storm_data$max_wind_in_box), na.rm = TRUE) >= 20, "(would flag as TC)", "(OK)")))
    
    meets_all_criteria <- (genesis_coords[2] < 30) & 
      (result$storm_data$`min_p (Pa)` <= 98000) & 
      (max(unlist(result$storm_data$max_wind_in_box), na.rm = TRUE) >= 20)
    print(paste("    Meets ALL 3 TC criteria?", meets_all_criteria,
                ifelse(meets_all_criteria, "⚠️ (SHOULD HAVE BEEN REMOVED UNLESS IT OCCURS OUTSIDE OF SON)", "(Correctly retained)")))
    
  } else {
    print("✗ REMOVED by TC filter")
    # Try to find it in the original to see why
    result_orig <- find_storm_by_id(cfsr_noreasters_wTC, storm_id)
    if (!is.null(result_orig)) {
      print("  Checking why it was removed:")
      
      cx_p_vals <- extract_coordinates_robust(result_orig$storm_data$cx_p)
      cy_p_vals <- extract_coordinates_robust(result_orig$storm_data$cy_p)
      cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
      cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
      genesis_coords <- lat_lon(cx_p_first, cy_p_first, cfsr_lon, cfsr_lat)
      
      print(paste("    Genesis lat:", round(genesis_coords[2], 2), "°N (< 30°N?", genesis_coords[2] < 30, ")"))
      print(paste("    Min pressure:", result_orig$storm_data$`min_p (Pa)`, "Pa (<= 98000?", 
                  result_orig$storm_data$`min_p (Pa)` <= 98000, ")"))
      print(paste("    Max wind:", result_orig$storm_data$max_wind_in_box, "m/s (>= 20?", 
                  result_orig$storm_data$max_wind_in_box >= 20, ")"))
    }
  }
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Plot CFSR nor'easter tracks to verify that each ID corresponds to the
#   correct storm

# Open PDF for plotting
pdf("/Path/for/output/figures/fig_s3_noreaster_track_verification_cfsr.pdf", 
    width = 14, height = 16)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 2), oma = c(1, 1, 1, 1))

# Loop through each noreaster and plot
storm_count <- 0
for (storm_name in names(cfsr_famous_noreasters)) {
  storm_id <- cfsr_famous_noreasters[[storm_name]]
  
  # Find the storm in the data
  found <- FALSE
  for (year in names(cfsr_noreasters_wTC)) {
    year_data <- cfsr_noreasters_wTC[[year]]
    if (is.data.frame(year_data)) {
      match_idx <- which(year_data$ID == storm_id)
      if (length(match_idx) > 0) {
        storm_data <- year_data[match_idx, ]
        
        cat(paste("Plotting:", storm_name, "from year", year, "\n"))
        
        # Plot the storm
        plot_cfsr_storm_track(storm_data, storm_name, cfsr_lon, cfsr_lat)
        
        storm_count <- storm_count + 1
        found <- TRUE
        break
      }
    }
  }
  
  if (!found) {
    cat(paste("WARNING: Storm not found:", storm_name, "-", storm_id, "\n"))
    
    # Create empty plot with message
    plot(1, type = "n", xlim = c(-100, -10), ylim = c(10, 70),
         xlab = "", ylab = "", main = storm_name, asp = 1.3)
    map("worldHires", add = TRUE, col = "gray90", fill = TRUE, border = "gray60")
    text(-55, 40, paste("Storm not found:\n", storm_id), 
         cex = 1.2, col = "red", font = 2)
  }
}

dev.off()

cat("\nNor'easter track plots saved to: noreaster_track_verification_cfsr.pdf\n")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# OVERALL FILTER STATISTICS: CFSR NOR'EASTERS
print("\n=== OVERALL FILTER STATISTICS ===")

total_before_TC <- sum(sapply(cfsr_noreasters_wTC, nrow))
total_after_TC <- sum(sapply(cfsr_noreasters_noTC, nrow))

print(paste("Total storms before TC filter:", total_before_TC))
print(paste("Total storms after TC filter:", total_after_TC))
print(paste("Storms removed by TC filter:", total_before_TC - total_after_TC))
print(paste("Percent removed by TC filter:", 
            round(100 * (total_before_TC - total_after_TC) / total_before_TC, 2), "%"))

# ==============================================================================

# ERA5 NOR'EASTERS

# Nor'easter IDs manually identified by matching dates, genesis locations, 
#   storm lifetimes, and maximum intensity values between tracked cyclone 
#   dataset and documented characteristics for each storm
era5_famous_noreasters <- list(
  "Great Appalachian Storm (1950)" = "11251950-00-40-56",
  "Perfect Storm (1991)" = "10291991-00-131-33",
  "Nor'Ida (2009)" = "11072009-12-17-120",
  "Halloween nor'easter (2011)" = "10292011-12-60-55",
  "November 2012 nor'easter (2012)" = "11062012-18-59-67"
)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Search for famous nor'easters in ERA5 cyclone dataset both BEFORE and AFTER
#   the TC filter is applied

# Search for each storm in cyclone dataset BEFORE TC filter 
#   (with potential TCs included)
print("--- BEFORE TC FILTER (era5_noreasters_wTC) ---")
for (storm_name in names(era5_famous_noreasters)) {
  storm_id <- era5_famous_noreasters[[storm_name]]
  result <- find_storm_by_id(era5_noreasters_wTC, storm_id)
  
  print(paste("\nSearching for:", storm_name, "-", storm_id))
  if (!is.null(result)) {
    print(paste("✓ Found in year:", result$year))
    # ^^^ cyclones are keyed by start year of the nor'easter season (i.e., a fall 1991 storm is keyed to the 1991 season)
    print(paste("  Min pressure:", result$storm_data$`min_p_h (Pa)`, "Pa"))
    print(paste("  Max wind:", result$storm_data$`max_w_all_h_750 (m/s)`, "m/s"))
    
    # Extract genesis position
    cx_p_vals <- extract_coordinates_robust(result$storm_data$cx_p)
    cy_p_vals <- extract_coordinates_robust(result$storm_data$cy_p)
    cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
    cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
    genesis_coords <- lat_lon(cx_p_first, cy_p_first, era5_lon, era5_lat)
    
    print(paste("  Genesis latitude:", round(genesis_coords[2], 2), "°N"))
    print(paste("  Genesis longitude:", round(genesis_coords[1], 2), "°E"))
    
  } else {
    print("✗ NOT FOUND in pre-filtered data")
  }
}

# Search for each storm in cyclone dataset AFTER TC filter  (TCs removed)
print("\n--- AFTER TC FILTER (era5_noreasters_noTC) ---")
for (storm_name in names(era5_famous_noreasters)) {
  storm_id <- era5_famous_noreasters[[storm_name]]
  result <- find_storm_by_id(era5_noreasters_noTC, storm_id)
  
  print(paste("\nSearching for:", storm_name, "-", storm_id))
  if (!is.null(result)) {
    print(paste("✓ RETAINED - Found in year:", result$year))
    print(paste("  Min pressure:", result$storm_data$`min_p_h (Pa)`, "Pa"))
    print(paste("  Max wind:", result$storm_data$`max_w_all_h_750 (m/s)`, "m/s"))
    
    # Extract genesis position
    cx_p_vals <- extract_coordinates_robust(result$storm_data$cx_p)
    cy_p_vals <- extract_coordinates_robust(result$storm_data$cy_p)
    cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
    cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
    genesis_coords <- lat_lon(cx_p_first, cy_p_first, era5_lon, era5_lat)
    
    print(paste("  Genesis latitude:", round(genesis_coords[2], 2), "°N"))
    print(paste("  Genesis longitude:", round(genesis_coords[1], 2), "°E"))
    
    # Check against TC criteria
    print("\n  TC Filter Criteria Check:")
    print(paste("    Genesis south of 30°N?", genesis_coords[2] < 30, 
                ifelse(genesis_coords[2] < 30, "(would flag as TC)", "(OK)")))
    print(paste("    Min pressure <= 98000 Pa?", result$storm_data$`min_p_h (Pa)` <= 98000,
                ifelse(result$storm_data$`min_p_h (Pa)` <= 98000, "(would flag as TC)", "(OK)")))
    print(paste("    Max wind >= 20 m/s?", result$storm_data$`max_w_all_h_750 (m/s)` >= 20,
                ifelse(result$storm_data$`max_w_all_h_750 (m/s)` >= 20, "(would flag as TC)", "(OK)")))
    
    meets_all_criteria <- (genesis_coords[2] < 30) & 
      (result$storm_data$`min_p_h (Pa)` <= 98000) & 
      (result$storm_data$`max_w_all_h_750 (m/s)` >= 20)
    print(paste("    Meets ALL 3 TC criteria?", meets_all_criteria,
                ifelse(meets_all_criteria, "⚠️ (SHOULD HAVE BEEN REMOVED UNLESS IT OCCURS OUTSIDE OF SON)", "(Correctly retained)")))
    
  } else {
    print("✗ REMOVED by TC filter")
    # Try to find it in the original to see why
    result_orig <- find_storm_by_id(era5_noreasters_wTC, storm_id)
    if (!is.null(result_orig)) {
      print("  Checking why it was removed:")
      
      cx_p_vals <- extract_coordinates_robust(result_orig$storm_data$cx_p)
      cy_p_vals <- extract_coordinates_robust(result_orig$storm_data$cy_p)
      cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
      cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
      genesis_coords <- lat_lon(cx_p_first, cy_p_first, era5_lon, era5_lat)
      
      print(paste("    Genesis lat:", round(genesis_coords[2], 2), "°N (< 30°N?", genesis_coords[2] < 30, ")"))
      print(paste("    Min pressure:", result_orig$storm_data$`min_p_h (Pa)`, "Pa (<= 98000?", 
                  result_orig$storm_data$`min_p_h (Pa)` <= 98000, ")"))
      print(paste("    Max wind:", result_orig$storm_data$`max_w_all_h_750 (m/s)`, "m/s (>= 20?", 
                  result_orig$storm_data$`max_w_all_h_750 (m/s)` >= 20, ")"))
    }
  }
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Plot ERA5 nor'easter tracks to verify that each ID corresponds to the
#   correct storm

# Open PDF for plotting
pdf("/Path/for/output/figures/fig_s4_noreaster_track_verification_era5.pdf", 
    width = 14, height = 16)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 2), oma = c(1, 1, 1, 1))

# Loop through each noreaster and plot
storm_count <- 0
for (storm_name in names(era5_famous_noreasters)) {
  storm_id <- era5_famous_noreasters[[storm_name]]
  
  # Find the storm in the data
  found <- FALSE
  for (year in names(era5_noreasters_wTC)) {
    year_data <- era5_noreasters_wTC[[year]]
    if (is.data.frame(year_data)) {
      match_idx <- which(year_data$ID == storm_id)
      if (length(match_idx) > 0) {
        storm_data <- year_data[match_idx, ]
        
        cat(paste("Plotting:", storm_name, "from year", year, "\n"))
        
        # Plot the storm
        plot_era5_storm_track(storm_data, storm_name, era5_lon, era5_lat)
        
        storm_count <- storm_count + 1
        found <- TRUE
        break
      }
    }
  }
  
  if (!found) {
    cat(paste("WARNING: Storm not found:", storm_name, "-", storm_id, "\n"))
    
    # Create empty plot with message
    plot(1, type = "n", xlim = c(-100, -10), ylim = c(10, 70),
         xlab = "", ylab = "", main = storm_name, asp = 1.3)
    map("worldHires", add = TRUE, col = "gray90", fill = TRUE, border = "gray60")
    text(-55, 40, paste("Storm not found:\n", storm_id), 
         cex = 1.2, col = "red", font = 2)
  }
}

dev.off()

cat("\nNor'easter track plots saved to: noreaster_track_verification_era5.pdf\n")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# OVERALL FILTER STATISTICS: ERA5 NOR'EASTERS
print("\n=== OVERALL FILTER STATISTICS ===")

total_before_TC <- sum(sapply(era5_noreasters_wTC, nrow))
total_after_TC <- sum(sapply(era5_noreasters_noTC, nrow))

print(paste("Total storms before TC filter:", total_before_TC))
print(paste("Total storms after TC filter:", total_after_TC))
print(paste("Storms removed by TC filter:", total_before_TC - total_after_TC))
print(paste("Percent removed by TC filter:", 
            round(100 * (total_before_TC - total_after_TC) / total_before_TC, 2), "%"))

# ==============================================================================

# CFSR TROPICAL CYCLONES

# Tropical cyclone IDs manually identified by matching dates, genesis locations, 
#   storm lifetimes, and maximum intensity values between tracked cyclone 
#   dataset and documented characteristics for each storm
cfsr_notable_tropical_cyclones <- list(
  "Hurricane Frederic (1979)" = "09081979-18-39-96", # Cat 4, merged w/ ETC ~9/14
  "Hurricane Gloria (1985)" = "09181985-06-132-110", # Cat 4, extratropical transition ~9/28
  "Hurricane Opal (1995)" = "09261995-06-29-105", # Cat 1, extratropical transition ~10/5
  "Hurricane Ivan (2004)" = "09042004-06-127-120", # Cat 3, extratropical transition ~9/18
  "Hurricane Sandy (2012)" = "10222012-06-45-113" # Cat 1 / merged w/ nor'easter
)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Search for notable tropical cyclones in CFSR cyclone dataset both BEFORE and AFTER
#   the TC filter is applied

# Search for each storm in cyclone dataset BEFORE TC filter 
#   (with potential TCs included)
print("--- BEFORE TC FILTER (cfsr_noreasters_wTC) ---")
for (storm_name in names(cfsr_notable_tropical_cyclones)) {
  storm_id <- cfsr_notable_tropical_cyclones[[storm_name]]
  result <- find_storm_by_id(cfsr_noreasters_wTC, storm_id)
  
  print(paste("\nSearching for:", storm_name, "-", storm_id))
  if (!is.null(result)) {
    print(paste("✓ Found in year:", result$year))
    print(paste("  Min pressure:", result$storm_data$`min_p (Pa)`, "Pa"))
    print(paste("  Max wind:", max(unlist(result$storm_data$max_wind_in_box), na.rm = TRUE), "m/s"))
    
    # Extract genesis position
    cx_p_vals <- extract_coordinates_robust(result$storm_data$cx_p)
    cy_p_vals <- extract_coordinates_robust(result$storm_data$cy_p)
    cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
    cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
    genesis_coords <- lat_lon(cx_p_first, cy_p_first, cfsr_lon, cfsr_lat)
    
    print(paste("  Genesis latitude:", round(genesis_coords[2], 2), "°N"))
    print(paste("  Genesis longitude:", round(genesis_coords[1], 2), "°E"))
    
  } else {
    print("✗ NOT FOUND in pre-filtered data")
  }
}

# Search for each storm in cyclone dataset AFTER TC filter (TCs removed)
print("\n--- AFTER TC FILTER (cfsr_noreasters_noTC) ---")
for (storm_name in names(cfsr_notable_tropical_cyclones)) {
  storm_id <- cfsr_notable_tropical_cyclones[[storm_name]]
  result <- find_storm_by_id(cfsr_noreasters_noTC, storm_id)
  
  print(paste("\nSearching for:", storm_name, "-", storm_id))
  if (!is.null(result)) {
    print(paste("✓ RETAINED - Found in year:", result$year))
    print(paste("  Min pressure:", result$storm_data$`min_p (Pa)`, "Pa"))
    print(paste("  Max wind:", max(unlist(result$storm_data$max_wind_in_box), na.rm = TRUE), "m/s"))
    
    # Extract genesis position
    cx_p_vals <- extract_coordinates_robust(result$storm_data$cx_p)
    cy_p_vals <- extract_coordinates_robust(result$storm_data$cy_p)
    cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
    cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
    genesis_coords <- lat_lon(cx_p_first, cy_p_first, cfsr_lon, cfsr_lat)
    
    print(paste("  Genesis latitude:", round(genesis_coords[2], 2), "°N"))
    print(paste("  Genesis longitude:", round(genesis_coords[1], 2), "°E"))
    
    # Check against TC criteria
    print("\n  TC Filter Criteria Check:")
    print(paste("    Genesis south of 30°N?", genesis_coords[2] < 30, 
                ifelse(genesis_coords[2] < 30, "(would flag as TC)", "(OK)")))
    print(paste("    Min pressure <= 98000 Pa?", result$storm_data$`min_p (Pa)` <= 98000,
                ifelse(result$storm_data$`min_p (Pa)` <= 98000, "(would flag as TC)", "(OK)")))
    print(paste("    Max wind >= 20 m/s?", max(unlist(result$storm_data$max_wind_in_box), na.rm = TRUE) >= 20,
                ifelse(max(unlist(result$storm_data$max_wind_in_box), na.rm = TRUE) >= 20, "(would flag as TC)", "(OK)")))
    
    meets_all_criteria <- (genesis_coords[2] < 30) & 
      (result$storm_data$`min_p (Pa)` <= 98000) & 
      (max(unlist(result$storm_data$max_wind_in_box), na.rm = TRUE) >= 20)
    print(paste("    Meets ALL 3 TC criteria?", meets_all_criteria,
                ifelse(meets_all_criteria, "⚠️ (SHOULD HAVE BEEN REMOVED UNLESS IT OCCURS OUTSIDE OF SON)", "(Correctly retained)")))
    
  } else {
    print("✗ REMOVED by TC filter")
    # Try to find it in the original to see why
    result_orig <- find_storm_by_id(cfsr_noreasters_wTC, storm_id)
    if (!is.null(result_orig)) {
      print("  Checking why it was removed:")
      
      cx_p_vals <- extract_coordinates_robust(result_orig$storm_data$cx_p)
      cy_p_vals <- extract_coordinates_robust(result_orig$storm_data$cy_p)
      cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
      cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
      genesis_coords <- lat_lon(cx_p_first, cy_p_first, cfsr_lon, cfsr_lat)
      
      print(paste("    Genesis lat:", round(genesis_coords[2], 2), "°N (< 30°N?", genesis_coords[2] < 30, ")"))
      print(paste("    Min pressure:", result_orig$storm_data$`min_p (Pa)`, "Pa (<= 98000?", 
                  result_orig$storm_data$`min_p (Pa)` <= 98000, ")"))
      print(paste("    Max wind:", max(unlist(result_orig$storm_data$max_wind_in_box), na.rm = TRUE), "m/s (>= 20?", 
                  max(unlist(result_orig$storm_data$max_wind_in_box), na.rm = TRUE) >= 20, ")"))
    }
  }
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Plot CFSR tropical cyclone tracks to verify that each ID corresponds to the
#   correct storm

pdf("/Path/for/output/figures/fig_s5_TC_track_verification_cfsr.pdf", 
    width = 14, height = 16)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 2), oma = c(1, 1, 1, 1))

# Loop through each tropical cyclone and plot
storm_count <- 0
for (storm_name in names(cfsr_notable_tropical_cyclones)) {
  storm_id <- cfsr_notable_tropical_cyclones[[storm_name]]
  
  # Find the storm in the data
  found <- FALSE
  for (year in names(cfsr_noreasters_wTC)) {
    year_data <- cfsr_noreasters_wTC[[year]]
    if (is.data.frame(year_data)) {
      match_idx <- which(year_data$ID == storm_id)
      if (length(match_idx) > 0) {
        storm_data <- year_data[match_idx, ]
        
        cat(paste("Plotting:", storm_name, "from year", year, "\n"))
        
        # Plot the storm
        plot_cfsr_storm_track(storm_data, storm_name, cfsr_lon, cfsr_lat)
        
        storm_count <- storm_count + 1
        found <- TRUE
        break
      }
    }
  }
  
  if (!found) {
    cat(paste("WARNING: Storm not found:", storm_name, "-", storm_id, "\n"))
    
    # Create empty plot with message
    plot(1, type = "n", xlim = c(-100, -10), ylim = c(10, 70),
         xlab = "", ylab = "", main = storm_name, asp = 1.3)
    map("worldHires", add = TRUE, col = "gray90", fill = TRUE, border = "gray60")
    text(-60, 40, paste("Storm not found:\n", storm_id), 
         cex = 1.2, col = "red", font = 2)
  }
}

dev.off()

cat("\nTropical cyclone track plots saved to: TC_track_verification_cfsr.pdf\n")

# ==============================================================================

# ERA5 TROPICAL CYCLONES

# Tropical cyclone IDs manually identified by matching dates, genesis locations, 
#   storm lifetimes, and maximum intensity values between tracked cyclone 
#   dataset and documented characteristics for each storm
era5_notable_tropical_cyclones <- list(
  "Hurricane Frederic (1979)" = "09061979-06-72-120", # Cat 4, merged w/ ETC ~9/14
  "Hurricane Gloria (1985)" = "09221985-18-123-120", # Cat 4, extratropical transition ~9/28
  "Hurricane Opal (1995)" = "09261995-00-14-120", # Cat 1, extratropical transition ~10/5
  "Hurricane Ivan (2004)" = "09102004-18-51-120", # Cat 3, extratropical transition ~9/18
  "Hurricane Sandy (2012)" = "10222012-12-65-120" # Cat 1 / merged w/ nor'easter
)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Search for notable tropical cyclones in ERA5 cyclone dataset both BEFORE and AFTER
#   the TC filter is applied

# Search for each storm in cyclone dataset BEFORE TC filter 
#   (with potential TCs included)
print("--- BEFORE TC FILTER (era5_noreasters_wTC) ---")
for (storm_name in names(era5_notable_tropical_cyclones)) {
  storm_id <- era5_notable_tropical_cyclones[[storm_name]]
  result <- find_storm_by_id(era5_noreasters_wTC, storm_id)
  
  print(paste("\nSearching for:", storm_name, "-", storm_id))
  if (!is.null(result)) {
    print(paste("✓ Found in year:", result$year))
    print(paste("  Min pressure:", result$storm_data$`min_p_h (Pa)`, "Pa"))
    print(paste("  Max wind:", result$storm_data$`max_w_all_h_750 (m/s)`, "m/s"))
    
    # Extract genesis position
    cx_p_vals <- extract_coordinates_robust(result$storm_data$cx_p)
    cy_p_vals <- extract_coordinates_robust(result$storm_data$cy_p)
    cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
    cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
    genesis_coords <- lat_lon(cx_p_first, cy_p_first, era5_lon, era5_lat)
    
    print(paste("  Genesis latitude:", round(genesis_coords[2], 2), "°N"))
    print(paste("  Genesis longitude:", round(genesis_coords[1], 2), "°E"))
    
  } else {
    print("✗ NOT FOUND in pre-filtered data")
  }
}

# Search for each storm in cyclone dataset AFTER TC filter (TCs removed)
print("\n--- AFTER TC FILTER (era5_noreasters_noTC) ---")
for (storm_name in names(era5_notable_tropical_cyclones)) {
  storm_id <- era5_notable_tropical_cyclones[[storm_name]]
  result <- find_storm_by_id(era5_noreasters_noTC, storm_id)
  
  print(paste("\nSearching for:", storm_name, "-", storm_id))
  if (!is.null(result)) {
    print(paste("✓ RETAINED - Found in year:", result$year))
    print(paste("  Min pressure:", result$storm_data$`min_p_h (Pa)`, "Pa"))
    print(paste("  Max wind:", result$storm_data$`max_w_all_h_750 (m/s)`, "m/s"))
    
    # Extract genesis position
    cx_p_vals <- extract_coordinates_robust(result$storm_data$cx_p)
    cy_p_vals <- extract_coordinates_robust(result$storm_data$cy_p)
    cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
    cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
    genesis_coords <- lat_lon(cx_p_first, cy_p_first, era5_lon, era5_lat)
    
    print(paste("  Genesis latitude:", round(genesis_coords[2], 2), "°N"))
    print(paste("  Genesis longitude:", round(genesis_coords[1], 2), "°E"))
    
    # Check against TC criteria
    print("\n  TC Filter Criteria Check:")
    print(paste("    Genesis south of 30°N?", genesis_coords[2] < 30, 
                ifelse(genesis_coords[2] < 30, "(would flag as TC)", "(OK)")))
    print(paste("    Min pressure <= 98000 Pa?", result$storm_data$`min_p_h (Pa)` <= 98000,
                ifelse(result$storm_data$`min_p_h (Pa)` <= 98000, "(would flag as TC)", "(OK)")))
    print(paste("    Max wind >= 20 m/s?", result$storm_data$`max_w_all_h_750 (m/s)` >= 20,
                ifelse(result$storm_data$`max_w_all_h_750 (m/s)` >= 20, "(would flag as TC)", "(OK)")))
    
    meets_all_criteria <- (genesis_coords[2] < 30) & 
      (result$storm_data$`min_p_h (Pa)` <= 98000) & 
      (result$storm_data$`max_w_all_h_750 (m/s)` >= 20)
    print(paste("    Meets ALL 3 TC criteria?", meets_all_criteria,
                ifelse(meets_all_criteria, "⚠️ (SHOULD HAVE BEEN REMOVED UNLESS IT OCCURS OUTSIDE OF SON)", "(Correctly retained)")))
    
  } else {
    print("✗ REMOVED by TC filter")
    # Try to find it in the original to see why
    result_orig <- find_storm_by_id(era5_noreasters_wTC, storm_id)
    if (!is.null(result_orig)) {
      print("  Checking why it was removed:")
      
      cx_p_vals <- extract_coordinates_robust(result_orig$storm_data$cx_p)
      cy_p_vals <- extract_coordinates_robust(result_orig$storm_data$cy_p)
      cx_p_first <- fix_corrupted_values(cx_p_vals[1])[1]
      cy_p_first <- fix_corrupted_values(cy_p_vals[1])[1]
      genesis_coords <- lat_lon(cx_p_first, cy_p_first, era5_lon, era5_lat)
      
      print(paste("    Genesis lat:", round(genesis_coords[2], 2), "°N (< 30°N?", genesis_coords[2] < 30, ")"))
      print(paste("    Min pressure:", result_orig$storm_data$`min_p_h (Pa)`, "Pa (<= 98000?", 
                  result_orig$storm_data$`min_p_h (Pa)` <= 98000, ")"))
      print(paste("    Max wind:", result_orig$storm_data$`max_w_all_h_750 (m/s)`, "m/s (>= 20?", 
                  result_orig$storm_data$`max_w_all_h_750 (m/s)` >= 20, ")"))
    }
  }
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Plot ERA5 tropical cyclone tracks to verify that each ID corresponds to the
#   correct storm

pdf("/Path/for/output/figures/fig_s6_TC_track_verification_era5.pdf", 
    width = 14, height = 16)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 2), oma = c(1, 1, 1, 1))

# Loop through each tropical cyclone and plot
storm_count <- 0
for (storm_name in names(era5_notable_tropical_cyclones)) {
  storm_id <- era5_notable_tropical_cyclones[[storm_name]]
  
  # Find the storm in the data
  found <- FALSE
  for (year in names(era5_noreasters_wTC)) {
    year_data <- era5_noreasters_wTC[[year]]
    if (is.data.frame(year_data)) {
      match_idx <- which(year_data$ID == storm_id)
      if (length(match_idx) > 0) {
        storm_data <- year_data[match_idx, ]
        
        cat(paste("Plotting:", storm_name, "from year", year, "\n"))
        
        # Plot the storm
        plot_era5_storm_track(storm_data, storm_name, era5_lon, era5_lat)
        
        storm_count <- storm_count + 1
        found <- TRUE
        break
      }
    }
  }
  
  if (!found) {
    cat(paste("WARNING: Storm not found:", storm_name, "-", storm_id, "\n"))
    
    # Create empty plot with message
    plot(1, type = "n", xlim = c(-100, -10), ylim = c(10, 70),
         xlab = "", ylab = "", main = storm_name, asp = 1.3)
    map("worldHires", add = TRUE, col = "gray90", fill = TRUE, border = "gray60")
    text(-60, 40, paste("Storm not found:\n", storm_id), 
         cex = 1.2, col = "red", font = 2)
  }
}

dev.off()

cat("\nTropical cyclone track plots saved to: TC_track_verification_era5.pdf\n")









