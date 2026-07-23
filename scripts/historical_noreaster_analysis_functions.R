# ==============================================================================

# Collection of functions used in analyses of historical nor'easter intensity
#   characteristics in CFSR, ERA5, CRCM5-ERA5, and CRCM5-CMIP6 datasets

# ==============================================================================

# Associated publication: Weaver, M.M. et al. (2026). Historical nor'easter
#                         intensity trends in reanalysis and regional climate
#                         simulations: A cross-dataset perspective.

# This file contains all functions necessary to run hist_noreaster_intensity_analysis.R,
#   which contains code to recreate the main text figures, tables, and 
#   supplementary figures in the aforementioned publication.

# ***Functions are listed in the order of their use in 
#    hist_noreaster_intensity_analysis.R, followed by order of use in
#    TC_filter_hist_noreaster_intensity.R***

# ==============================================================================

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# Function to extract year from filename (storm tracks)
extract_year <- function(filepath) {
  filename <- basename(filepath)
  # Handle both {YYYY} and YYYY patterns
  year_match <- str_extract(filename, "\\{?(\\d{4})\\}?")
  year_match <- str_remove_all(year_match, "[{}]")
  return(as.numeric(year_match))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# Function to properly read storm CSV files
read_storm_csv <- function(filepath) {
  # Read all lines as text
  lines <- readLines(filepath, warn = FALSE)
  
  # Remove empty lines
  lines <- lines[nzchar(lines)]
  
  if (length(lines) < 2) {
    warning(paste("File", basename(filepath), "has insufficient data"))
    return(data.frame())
  }
  
  # The header is the first line
  header <- lines[1]
  
  # Reconstruct proper data rows
  reconstructed_rows <- c()
  current_row <- ""
  
  for (i in 2:length(lines)) {
    line <- lines[i]
    
    # Check if this line starts a new storm (begins with a number and comma)
    if (grepl("^[0-9]+,", line)) {
      # If we have a current row being built, save it
      if (nzchar(current_row)) {
        reconstructed_rows <- c(reconstructed_rows, current_row)
      }
      # Start new row
      current_row <- line
    } else {
      # This is a continuation line - append to current row
      current_row <- paste0(current_row, " ", trimws(line))
    }
  }
  
  # Don't forget the last row
  if (nzchar(current_row)) {
    reconstructed_rows <- c(reconstructed_rows, current_row)
  }
  
  # Now create a proper CSV string
  csv_content <- paste(c(header, reconstructed_rows), collapse = "\n")
  
  # Read using fread on the reconstructed content
  tryCatch({
    df <- fread(text = csv_content, 
                stringsAsFactors = FALSE, 
                na.strings = c("", "NA", "NULL", "-999", "9999"))
    
    # CRITICAL: Ensure we actually have a data frame
    if (!is.data.frame(df)) {
      # Convert to data frame if it's not already
      df <- as.data.frame(df)
    }
    
    # Remove entirely empty rows
    if (nrow(df) > 0) {
      df <- df[!apply(is.na(df), 1, all), ]
    }
    
    # Verify we have the expected structure
    if (nrow(df) > 0 && ncol(df) > 1) {
      return(df)
    } else {
      warning(paste("File", basename(filepath), "resulted in invalid data frame structure"))
      return(data.frame())
    }
    
  }, error = function(e) {
    warning(paste("Error parsing reconstructed CSV for", basename(filepath), ":", e$message))
    return(data.frame())
  })
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# Create CFSR grid (0.5 degree resolution)
create_cfsr_grid <- function() {
  lon_seq <- seq(-100, -10, by = 0.5)
  lat_seq <- rev(seq(10, 70, by = 0.5))  # 70N, 69.5N, ..., 10N (flipped so that
                                         # lats are indexed correctly later)
  # Create coordinate matrices
  # lon_grid should vary by column, lat_grid should vary by row
  lon_grid <- matrix(rep(lon_seq, each = length(lat_seq)), 
                     nrow = length(lat_seq), ncol = length(lon_seq))
  lat_grid <- matrix(rep(lat_seq, times = length(lon_seq)), 
                     nrow = length(lat_seq), ncol = length(lon_seq))
  
  return(list(
    lon = lon_grid,
    lat = lat_grid,
    lon_seq = lon_seq,
    lat_seq = lat_seq
  ))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Create ERA5 grid (0.25 degree resolution)
create_era5_grid <- function() {
  # -90 to -40 longitude (match Chen et al., 2025)
  lon_seq <- seq(-90, -40, by = 0.25)
  # 20 to 50 latitude (match Chen et al., 2025)
  lat_seq <- rev(seq(20, 50, by = 0.25))  # 50N, 49.75N, ..., 20N (flipped so 
                                          # lats are indexed correctly later)
  # Create coordinate matrices
  # lon_grid should vary by column, lat_grid should vary by row
  lon_grid <- matrix(rep(lon_seq, each = length(lat_seq)), 
                     nrow = length(lat_seq), ncol = length(lon_seq))
  lat_grid <- matrix(rep(lat_seq, times = length(lon_seq)), 
                     nrow = length(lat_seq), ncol = length(lon_seq))
  
  return(list(
    lon = lon_grid,
    lat = lat_grid,
    lon_seq = lon_seq,
    lat_seq = lat_seq
  ))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# count how many times a storm passes through the USEC box
pass_USEC <- function(xps_raw, yps_raw, xcs_raw, ycs_raw, lat_grid, lon_grid) {
  passed <- 0
  
  xps <- extract_coordinates_robust(xps_raw)
  yps <- extract_coordinates_robust(yps_raw)
  xcs <- extract_coordinates_robust(xcs_raw)
  ycs <- extract_coordinates_robust(ycs_raw)
  
  len <- min(length(xps), length(yps), length(xcs), length(ycs))
  if (len == 0) return(0)
  
  # Bounding box for U.S. East Coast
  usne_N = 45 # deg lat
  usne_S = 33 # deg lat
  usne_E = -70 # deg lon
  usne_W = -80 # deg lon
  
  for (i in seq_len(len)) {
    if (is.na(xps[i]) || is.na(yps[i]) || is.na(xcs[i]) || is.na(ycs[i])) next
    
    lonp_latp <- lat_lon(xps[i], yps[i], lon_grid, lat_grid)
    lonc_latc <- lat_lon(xcs[i], ycs[i], lon_grid, lat_grid)
    
    lonp <- lonp_latp[1]
    latp <- lonp_latp[2]
    lonc <- lonc_latc[1]
    latc <- lonc_latc[2]
    
    if (any(is.na(c(lonp, latp, lonc, latc)))) next
    
    p_in <- lonp >= usne_W && lonp <= usne_E &&
      latp >= usne_S && latp <= usne_N
    
    c_in <- lonc >= usne_W && lonc <= usne_E &&
      latc >= usne_S && latc <= usne_N
    
    if (p_in || c_in) {
      passed <- passed + 1
    }
  }
  
  return(passed)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# identify storms that pass through the USEC box
filter_USEC <- function(df_track, lat_grid, lon_grid) {
  passed_counts <- integer(nrow(df_track))
  
  for (i in seq_len(nrow(df_track))) {
    xps_raw <- df_track$cx_p[i]
    yps_raw <- df_track$cy_p[i]
    xcs_raw <- df_track$cx_c[i]
    ycs_raw <- df_track$cy_c[i]
    
    passed_counts[i] <- pass_USEC(xps_raw, yps_raw, xcs_raw, ycs_raw, lat_grid, lon_grid)
  }
  
  df_track$passed_USEC <- passed_counts
  return(df_track)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# keep only storms that pass through USEC at least once
filter_USECPassedOnly <- function(df_track) {
  df_track <- df_track[df_track$passed_USEC > 0, ]
  rownames(df_track) <- NULL
  return(df_track)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# remove storms that last less than 24 h
filter_time <- function(df_out) {
  df_out <- df_out[df_out[["lifetime (h)"]] >= 24, ]
  rownames(df_out) <- NULL
  return(df_out)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# remove storms that don't travel at least 1000 km
filter_dist <- function(df_out) {
  df_out <- df_out[
    df_out[["tot_dist_p (km)"]] >= 1000 & 
      df_out[["tot_dist_c (km)"]] >= 1000, ]
  rownames(df_out) <- NULL
  return(df_out)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to clean and standardize datetime/storm_id columns
clean_datetime <- function(dt) {
  # Handle different storm_id formats
  dt[, storm_id_clean := fifelse(
    grepl("^\\d{4}-\\d{2}-\\d{2}", storm_id), 
    storm_id,
    as.character(as.POSIXct(as.numeric(storm_id), origin = "1899-12-30"))
  )]
  
  # Handle different datetime formats
  dt[, datetime_clean := fifelse(
    grepl("^\\d{4}-\\d{2}-\\d{2}", datetime),
    as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S"),
    as.POSIXct(as.numeric(datetime), origin = "1899-12-30")
  )]
  
  return(dt)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to read and process a single wind data CSV file
read_wind_file <- function(file_path, model_name = NULL) {
  # Extract metadata from filename and path
  filename <- basename(file_path)
  dir_path <- dirname(file_path)
  
  # Parse filename
  parts <- strsplit(filename, "_")[[1]]
  end_year <- as.numeric(parts[1])
  season <- paste0(end_year - 1, "-", end_year)
  
  # Extract model and member info from directory structure
  path_parts <- strsplit(dir_path, "/")[[1]]
  
  # If model_name is provided, use it directly
  if (!is.null(model_name)) {
    base_model <- toupper(model_name)
    member <- "1"  # Default for single-member models
    model_full <- base_model
  } else {
    # Find the model directory (ends with _wind_data)
    model_dir_idx <- which(grepl("_wind_data$", path_parts))
    
    if (length(model_dir_idx) == 0) {
      stop("Could not find model directory in path: ", dir_path)
    }
    
    # Take the last match if multiple found
    model_dir_idx <- model_dir_idx[length(model_dir_idx)]
    model_dir <- path_parts[model_dir_idx]
    base_model <- gsub("_wind_data$", "", model_dir)
    
    # Handle special case for CRCM5-ERA5
    if (length(base_model) == 1 && base_model == "CRCM5-ERA5") {
      base_model <- "ERA5"
    }
    
    # Extract member information
    if (model_dir_idx < length(path_parts)) {
      member_dir <- path_parts[model_dir_idx + 1]
    } else {
      stop("Could not find member directory in path: ", dir_path)
    }
    
    if (base_model %in% c("MPI", "ERA5")) {
      # Multi-member models
      if (base_model == "MPI") {
        member <- gsub("historical", "", member_dir)
      } else { # ERA5
        member <- gsub("v1-r", "", member_dir)
      }
      model_full <- paste0(base_model, "_", member)
    } else {
      # Single-member models
      member <- "1"
      model_full <- base_model
    }
  }
  
  # Read the file
  dt <- fread(file_path)
  
  # Clean datetime columns
  dt <- clean_datetime(dt)
  
  # Add metadata columns
  dt[, `:=`(
    model = base_model,
    model_member = model_full,
    member = member,
    season = season,
    season_start_year = end_year - 1,
    season_end_year = end_year
  )]
  
  return(dt)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to merge CFSR wind data into existing storm data structure
merge_wind_into_cfsr_storms <- function(noreasters, wind_dirs) {
  
  require(data.table)
  
  cat("Merging wind data into storm tracks...\n")
  
  model_name <- "CFSR"  # fixed label since this function only ever handles CFSR data
  
  wind_dir <- wind_dirs
  
  if (!dir.exists(wind_dir)) {
    cat(paste("  Warning: Wind directory does not exist:", wind_dir, "\n"))
    return(noreasters)
  }
  
  # Get all wind CSV files
  wind_files <- list.files(wind_dir, pattern = "\\.csv$", 
                           full.names = TRUE, recursive = TRUE)
  
  if (length(wind_files) == 0) {
    cat(paste("  Warning: No wind CSV files found in:", wind_dir, "\n"))
    return(noreasters)
  }
  
  # Read all wind data
  cat(paste("  Reading", length(wind_files), "wind files for", model_name, "...\n"))
  
  wind_data_list <- list()
  for (i in seq_along(wind_files)) {
    if (i %% 10 == 0) {
      cat(paste("    Processing file", i, "of", length(wind_files), "\n"))
    }
    tryCatch({
      wind_data_list[[i]] <- read_wind_file(wind_files[i], model_name = model_name)
    }, error = function(e) {
      cat(paste("    ERROR in file:", basename(wind_files[i]), "\n"))
      cat(paste("    Error message:", e$message, "\n"))
    })
  }
  
  # Remove NULL entries
  wind_data_list <- wind_data_list[!sapply(wind_data_list, is.null)]
  
  if (length(wind_data_list) == 0) {
    cat(paste("  Error: No wind files could be read successfully\n"))
    return(noreasters)
  }
  
  cat(paste("  Successfully read", length(wind_data_list), "files. Now combining...\n"))
  
  model_wind_data <- data.table::rbindlist(wind_data_list, fill = TRUE)
  
  if (!is.data.table(model_wind_data)) {
    model_wind_data <- as.data.table(model_wind_data)
  }
  
  cat(paste("  Combined into", nrow(model_wind_data), "wind track points\n"))
  
  # Aggregate wind data by storm (one row per storm with list columns)
  cat("  Aggregating wind data by storm...\n")
  wind_aggregated_all <- model_wind_data[, .(
    u_wind = list(u_wind),
    v_wind = list(v_wind),
    wind_speed = list(wind_speed_center),
    wind_direction = list(wind_direction_center),
    max_wind_in_box = list(max_wind_in_box),
    max_wind_x_idx = list(max_wind_lon_idx),
    max_wind_y_idx = list(max_wind_lat_idx),
    rlat_idx = list(lat_idx),
    rlon_idx = list(lon_idx),
    storm_radius_km = list(storm_radius_km),
    grid_radius_cells = list(grid_radius_cells),
    datetime = list(datetime),
    season_start_year = first(season_start_year),
    season_end_year = first(season_end_year)
  ), by = storm_id]
  
  cat(paste("  Aggregated to", nrow(wind_aggregated_all), "unique storms\n"))
  
  # Process each year
  for (year in names(noreasters)) {
    storm_df <- noreasters[[year]]
    
    if (is.null(storm_df) || nrow(storm_df) == 0) next
    
    cat(paste("  Processing year", year, "...\n"))
    cat(paste("    Original storm count:", nrow(storm_df), "\n"))
    
    noreaster_ids <- storm_df$ID
    
    year_num <- as.numeric(year)
    year_wind <- wind_aggregated_all[
      storm_id %in% noreaster_ids &
        (season_start_year == year_num | season_end_year == year_num)
    ]
    
    if (nrow(year_wind) == 0) {
      cat(paste("    Warning: No matching wind data for year", year, "\n"))
      next
    }
    
    cat(paste("    Found", nrow(year_wind), "matching storms with wind data\n"))
    
    storm_dt <- as.data.table(storm_df)
    
    merged_dt <- merge(
      storm_dt,
      year_wind,
      by.x = "ID",
      by.y = "storm_id",
      all.x = TRUE,
      all.y = FALSE,
      sort = FALSE
    )
    
    if (nrow(merged_dt) != nrow(storm_df)) {
      cat(paste("    WARNING: Row count changed from", nrow(storm_df), 
                "to", nrow(merged_dt), "\n"))
      cat(paste("    This suggests duplicate IDs. Investigating...\n"))
      
      dup_ids <- storm_df$ID[duplicated(storm_df$ID)]
      if (length(dup_ids) > 0) {
        cat(paste("    Found", length(dup_ids), "duplicate IDs in storm data\n"))
      }
      
      dup_wind_ids <- year_wind$storm_id[duplicated(year_wind$storm_id)]
      if (length(dup_wind_ids) > 0) {
        cat(paste("    Found", length(dup_wind_ids), "duplicate IDs in wind data\n"))
      }
    }
    
    if (!is.data.table(storm_df)) {
      merged_dt <- as.data.frame(merged_dt)
    }
    
    noreasters[[year]] <- merged_dt
    
    n_total <- nrow(merged_dt)
    n_matched <- sum(!is.na(merged_dt$wind_speed))
    cat(paste("    Final: Year", year, ":", n_matched, "of", n_total, 
              "storms have wind data\n"))
  }
  
  cat("Wind data merge complete!\n")
  return(noreasters)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Filter TCs out of nor'easter dataset via genesis, min. SLP, & max. wind parameters
TC_filter_CFSR <- function(noreasters, lon_grid, lat_grid, min_lat = 30, min_slp = 98000, max_wind = 20) {
  
  cat("Filtering tropical cyclones from nor'easter dataset...\n")
  cat("Criteria: Genesis south of", min_lat, "°N, min pressure <=", min_slp, "Pa, max wind >=", max_wind, "m/s\n\n")
  
  filtered_noreasters <- list()
  
  # Process each year
  for (year in names(noreasters)) {
    cat(paste("Processing year:", year, "...\n"))
    
    # Get the dataframe for this year
    storm_df <- noreasters[[year]]
    
    # Check if year has data
    if (is.null(storm_df) || nrow(storm_df) == 0) {
      cat("  No storms found for this year\n")
      filtered_noreasters[[year]] <- storm_df
      next
    }
    
    original_count <- nrow(storm_df)
    
    # Extract month from storm ID (MMDDYYYY-HH-XX-XX)
    storm_id <- storm_df$ID
    storm_month <- substr(storm_id, 1, 2)
    
    # Identify Sep–Oct–Nov storms
    is_SON <- storm_month %in% c("09", "10", "11")
    n_SON <- sum(is_SON, na.rm = TRUE)
    
    # Calculate genesis latitude for each storm
    genesis_lat <- numeric(original_count)
    
    for (i in 1:original_count) {
      # Extract first position (genesis point)
      cx_p_first <- extract_coordinates_robust(storm_df[i, "cx_p"])[1]
      cy_p_first <- extract_coordinates_robust(storm_df[i, "cy_p"])[1]
      
      # Skip if coordinates are invalid
      if (is.na(cx_p_first) || is.na(cy_p_first) || 
          length(cx_p_first) == 0 || length(cy_p_first) == 0) {
        genesis_lat[i] <- NA
        next
      }
      
      # Apply corruption fix
      cx_p_first <- fix_corrupted_values(cx_p_first)[1]
      cy_p_first <- fix_corrupted_values(cy_p_first)[1]
      
      # Convert to lat/lon
      coords <- lat_lon(cx_p_first, cy_p_first, lon_grid, lat_grid)
      genesis_lat[i] <- coords[2]  # latitude is second element
    }
    
    # Get minimum pressure for each storm
    min_pressure <- storm_df$`min_p (Pa)`
    
    # Get maximum wind in box for each storm
    max_wind_per_storm <- sapply(storm_df$max_wind_in_box, function(wind_list) {
      if (is.null(wind_list) || length(wind_list) == 0) {
        return(NA)
      }
      max(wind_list, na.rm = TRUE)
    })
    
    # Define TC criteria (ONLY applied to Sep–Oct–Nov storms)
    is_TC <- is_SON &                          # Only SON storms evaluated
      (genesis_lat < min_lat) &                # Forms south of min_lat
      (min_pressure <= min_slp) &              # Reaches min pressure threshold
      (max_wind_per_storm >= max_wind)          # Reaches max wind threshold
    
    # Treat NAs as non-TC
    is_TC[is.na(is_TC)] <- FALSE
    
    # Filter out TCs - keep the dataframe structure
    filtered_df <- storm_df[!is_TC, ]
    
    # Reset row names
    rownames(filtered_df) <- NULL
    
    # Store filtered dataframe
    filtered_noreasters[[year]] <- filtered_df
    
    # Report statistics
    n_filtered <- sum(is_TC)
    n_kept <- sum(!is_TC)
    
    cat(paste(
      "  Year", year, ":",
      n_kept, "storms kept out of", original_count,
      "(removed", n_filtered, "TCs)\n"
    ))
    
    cat(paste(
      "    →", n_SON, "SON storms evaluated,",
      n_filtered, "removed\n"
    ))
  }
  
  cat("\nTC filtering complete!\n")
  return(filtered_noreasters)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Filter storms to those that passed through analysis region
#   ***not used in the present analysis, but useful if a more restricted
#      analysis region is desired***
filter_analysis_region <- function(df_track, lat_grid, lon_grid,
                                   lon_min = -90, lon_max = -40, 
                                   lat_min = 20, lat_max = 50) {
  keep_storms <- logical(nrow(df_track))
  
  for (i in seq_len(nrow(df_track))) {
    cx <- extract_coordinates_robust(df_track$cx_p[i])
    cy <- extract_coordinates_robust(df_track$cy_p[i])
    
    if (length(cx) > 0 && length(cy) > 0) {
      len <- min(length(cx), length(cy))
      
      # Convert to lat/lon for all points in track
      lons <- numeric(len)
      lats <- numeric(len)
      
      for (j in seq_len(len)) {
        if (!is.na(cx[j]) && !is.na(cy[j])) {
          lon_lat <- lat_lon(cx[j], cy[j], lon_grid, lat_grid)
          lons[j] <- lon_lat[1]
          lats[j] <- lon_lat[2]
        } else {
          lons[j] <- NA
          lats[j] <- NA
        }
      }
      
      # Keep storm if ANY point is in analysis region
      in_region <- (lons >= lon_min & lons <= lon_max &
                      lats >= lat_min & lats <= lat_max)
      
      keep_storms[i] <- any(in_region, na.rm = TRUE)
    }
  }
  
  df_filtered <- df_track[keep_storms, ]
  rownames(df_filtered) <- NULL
  return(df_filtered)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# Master filtering function for CFSR storm track data

# Structure of cfsr_noreasters:
#   List of 46 data.frames keyed by year
#   e.g. cfsr_noreasters$`1980` is a data.frame of all storms
#   tracked in the 1979-80 season

filter_cfsr_noreasters <- function(noreasters, lat_grid, lon_grid,
                                   wind_dirs = NULL, apply_TC_filter = FALSE,
                                   TC_min_lat = 30, TC_min_slp = 98000,
                                   TC_max_wind = 20, apply_region_filter = FALSE,
                                   lon_min = -90, lon_max = -40,
                                   lat_min = 20, lat_max = 50) {
  
  filtered_data <- list()
  
  for (year in names(noreasters)) {
    
    df <- noreasters[[year]]
    
    # Skip if not a valid data.frame or empty
    if (!is.data.frame(df) || nrow(df) == 0) {
      cat("Skipping year", year, "- empty or invalid data\n")
      next
    }
    
    # Apply core filters (same as CRCM5 version, functions unchanged)
    df <- df |>
      filter_USEC(lat_grid, lon_grid) |>
      filter_USECPassedOnly() |>
      filter_time() |>
      filter_dist()
    
    if (nrow(df) == 0) next
    
    filtered_data[[year]] <- df
  }
  
  # Merge wind data if wind directories provided
  # (same as CRCM5 version -- happens after core filters, before TC/region filters)
  if (!is.null(wind_dirs)) {
    cat("\nMerging wind data...\n")
    filtered_data <- merge_wind_into_cfsr_storms(filtered_data, wind_dirs)
  }
  
  # Apply TC filter if requested and wind data available
  if (apply_TC_filter) {
    if (!is.null(wind_dirs)) {
      cat("\nApplying TC filter...\n")
      filtered_data <- TC_filter_CFSR(
        filtered_data,
        lon_grid,
        lat_grid,
        min_lat  = TC_min_lat,
        min_slp  = TC_min_slp,
        max_wind = TC_max_wind
      )
    } else {
      cat("\nWarning: Cannot apply TC filter - wind data not available\n")
      cat("Provide wind_dirs to enable TC filtering\n")
    }
  }
  
  # Apply genesis/analysis region filter if requested
  if (apply_region_filter) {
    cat("\nApplying genesis region filter...\n")
    for (year in names(filtered_data)) {
      df <- filtered_data[[year]]
      if (!is.data.frame(df) || nrow(df) == 0) next
      df <- filter_analysis_region(df, lat_grid, lon_grid,
                                   lon_min, lon_max, lat_min, lat_max)
      if (nrow(df) > 0) {
        filtered_data[[year]] <- df
      } else {
        filtered_data[[year]] <- NULL
      }
    }
  }
  
  return(filtered_data)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Takes dataframe/list of storm tracks and drops storms that don't reach a min SLP
SLP_filter <- function(df_out, limit = 98000, lim_min = TRUE, era5 = FALSE, dataframe = TRUE) {
  
  if (dataframe == FALSE) {
    # List of vectors structure (ERA5 or CRCM5)
    # Determine which pressure variable to use
    pressure_var <- if (era5) "min_p_h (Pa)" else "min_p (Pa)"
    
    # Create a logical vector for filtering
    if (lim_min == TRUE) {
      # retain weak nor'easters
      keep <- df_out[[pressure_var]] > limit
    } else {
      # retain strong nor'easters
      keep <- df_out[[pressure_var]] <= limit
    }
    
    # Apply filter to all variables in the list
    filtered_list <- list()
    for (var_name in names(df_out)) {
      filtered_list[[var_name]] <- df_out[[var_name]][keep]
    }
    
    return(filtered_list)
    
  } else {
    # Dataframe structure (ERA5, CFSR, or CRCM5)
    # Determine which pressure variable to use
    if ("min_p_h (Pa)" %in% colnames(df)) {
      pressure_var <- "min_p_h (Pa)"
    } else {
      pressure_var <- "min_p (Pa)"
    }
    
    if (lim_min == TRUE) {
      # retain weak nor'easters
      df_out <- df_out[df_out[[pressure_var]] > limit, ]
    } else {
      # retain strong nor'easters
      df_out <- df_out[df_out[[pressure_var]] <= limit, ]
    }
    
    rownames(df_out) <- NULL
    return(df_out)
  }
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Wrapper function to apply SLP_filter to any data structure 
#   (CFSR, ERA5, CRCM5-ERA5, or CRCM5-CMIP6)
apply_SLP_filter <- function(noreasters, limit = 98000, lim_min = TRUE, 
                             era5 = FALSE, dataframe = TRUE) {
  filtered_data <- list()
  
  cat("Starting SLP filter with limit =", limit, "and lim_min =", lim_min, "\n")
  
  if (dataframe == FALSE) {
    # List of vectors structure (ERA5 or CRCM5)
    cat("Processing list of vectors structure\n")
    cat("Number of years:", length(names(noreasters)), "\n")
    filtered_list_data <- list()
    
    # Determine which pressure variable to use
    pressure_var <- if (era5) "min_p_h (Pa)" else "min_p (Pa)"
    
    for (year in names(noreasters)) {
      cat("  Processing year:", year, "\n")
      year_data <- noreasters[[year]]
      
      # Check if the pressure column exists - use names() for list structure
      if (!pressure_var %in% names(year_data)) {
        cat("    No", pressure_var, "column - keeping unchanged\n")
        filtered_list_data[[year]] <- year_data
        next
      }
      
      cat("    min_p range:", range(year_data[[pressure_var]], na.rm = TRUE), "\n")
      cat("    Number of NAs in min_p:", sum(is.na(year_data[[pressure_var]])), "\n")
      
      # Apply the SLP filter
      list_filtered <- SLP_filter(year_data, limit = limit, lim_min = lim_min, era5 = era5, dataframe = dataframe)
      
      # Check length instead of nrow since it returns a list
      if (length(list_filtered$ID) > 0) {
        filtered_list_data[[year]] <- list_filtered
      }
    }
    
    cat("\nFinal result has", length(filtered_list_data), "years\n")
    return(filtered_list_data)
    
  } else if (dataframe == TRUE && era5 == TRUE) {
    # ERA5 dataframe structure: noreasters$year = dataframe
    cat("Processing ERA5 or CFSR dataframe structure\n")
    cat("Number of years:", length(names(noreasters)), "\n")
    filtered_era5_data <- list()
    
    for (year in names(noreasters)) {
      cat("  Processing year:", year, "\n")
      df <- noreasters[[year]]
      
      # Skip if not a data frame or empty
      if (!is.data.frame(df) || nrow(df) == 0) {
        cat("    Skipping - not a data frame or empty\n")
        next
      }
      
      cat("    Rows:", nrow(df), "Columns:", ncol(df), "\n")
      
      if ("min_p_h (Pa)" %in% colnames(df)) {
        pressure_var <- "min_p_h (Pa)"
      } else {
        pressure_var <- "min_p (Pa)"
      }
      
      # Check if the min_p column exists
      if (!pressure_var %in% colnames(df)) {
        cat("    No", pressure_var, "column - keeping unchanged\n")
        filtered_era5_data[[year]] <- df
        next
      }
      
      cat("    min_p range:", range(df[[pressure_var]], na.rm = TRUE), "\n")
      cat("    Number of NAs in min_p:", sum(is.na(df[[pressure_var]])), "\n")
      
      # Troubleshoot CFSR issue (IDs reordered to begin in Jan., should be Sept.)
      cat("    First 3 IDs before filter:", head(df$ID, 3), "\n")
      df_filtered <- SLP_filter(df, limit = limit, lim_min = lim_min, era5 = era5, dataframe = dataframe)
      cat("    First 3 IDs after filter:", head(df_filtered$ID, 3), "\n")
      
      # Apply the SLP filter
      df_filtered <- SLP_filter(df, limit = limit, lim_min = lim_min, era5 = era5, dataframe = dataframe)
      
      cat("    Filtered rows:", nrow(df_filtered), "\n")
      
      if (nrow(df_filtered) > 0) {
        filtered_era5_data[[year]] <- df_filtered
      }
    }
    
    cat("\nFinal result has", length(filtered_era5_data), "years\n")
    return(filtered_era5_data)
    
  } else {
    # CRCM5 dataframe structure: noreasters$model$member/year$dataframe
    cat("Processing CRCM5 dataframe structure\n")
    pressure_var <- "min_p (Pa)"
    
    cat("Number of models:", length(names(noreasters)), "\n")
    
    for (model in names(noreasters)) {
      cat("\nProcessing model:", model, "\n")
      model_data <- noreasters[[model]]
      filtered_model_data <- list()
      
      for (member_or_year in names(model_data)) {
        cat("  Processing member/year:", member_or_year, "\n")
        member_data <- model_data[[member_or_year]]
        
        if (is.list(member_data) && !is.data.frame(member_data)) {
          # This should be ensemble data with years as keys
          cat("    Ensemble data detected\n")
          filtered_member_data <- list()
          
          for (year in names(member_data)) {
            cat("    Processing year:", year, "\n")
            df <- member_data[[year]]
            
            cat("      Data type:", class(df), "\n")
            cat("      Is data frame?", is.data.frame(df), "\n")
            
            # Skip if not a data frame or empty
            if (!is.data.frame(df) || nrow(df) == 0) {
              cat("      Skipping - not a data frame or empty\n")
              next
            }
            
            cat("      Rows:", nrow(df), "Columns:", ncol(df), "\n")
            cat("      Column names:", paste(colnames(df), collapse = ", "), "\n")
            
            # Check if the min_p column exists
            if (!pressure_var %in% colnames(df)) {
              cat("      No", pressure_var, "column - keeping unchanged\n")
              filtered_member_data[[year]] <- df
              next
            }
            
            cat("      min_p range:", range(df[[pressure_var]], na.rm = TRUE), "\n")
            cat("      Number of NAs in min_p:", sum(is.na(df[[pressure_var]])), "\n")
            
            # Apply the SLP filter
            df_filtered <- SLP_filter(df, limit = limit, lim_min = lim_min, era5 = era5, dataframe = dataframe)
            
            cat("      Filtered rows:", nrow(df_filtered), "\n")
            
            if (nrow(df_filtered) > 0) {
              filtered_member_data[[year]] <- df_filtered
            }
          }
          
          if (length(filtered_member_data) > 0) {
            filtered_model_data[[member_or_year]] <- filtered_member_data
          }
        } else {
          # This should be single-member data (direct data frame)
          cat("    Single-member data detected\n")
          df <- member_data
          
          cat("    Data type:", class(df), "\n")
          cat("    Is data frame?", is.data.frame(df), "\n")
          
          # Skip if not a data frame or empty
          if (!is.data.frame(df) || nrow(df) == 0) {
            cat("    Skipping - not a data frame or empty\n")
            next
          }
          
          cat("    Rows:", nrow(df), "Columns:", ncol(df), "\n")
          cat("    Column names:", paste(colnames(df), collapse = ", "), "\n")
          
          # Check if the min_p column exists
          if (!pressure_var %in% colnames(df)) {
            cat("    No", pressure_var, "column - keeping unchanged\n")
            filtered_model_data[[member_or_year]] <- df
            next
          }
          
          cat("    min_p range:", range(df[[pressure_var]], na.rm = TRUE), "\n")
          cat("    Number of NAs in min_p:", sum(is.na(df[[pressure_var]])), "\n")
          
          # Apply the SLP filter
          df_filtered <- SLP_filter(df, limit = limit, lim_min = lim_min, era5 = era5, dataframe = dataframe)
          
          cat("    Filtered rows:", nrow(df_filtered), "\n")
          
          if (nrow(df_filtered) > 0) {
            filtered_model_data[[member_or_year]] <- df_filtered
          }
        }
      }
      
      cat("  Model", model, "has", length(filtered_model_data), "filtered datasets\n")
      
      if (length(filtered_model_data) > 0) {
        filtered_data[[model]] <- filtered_model_data
      }
    }
    
    cat("\nFinal result has", length(filtered_data), "models\n")
    return(filtered_data)
  }
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# extract correct grid cell values from storm csv cx_p and cy_p columns
extract_coordinates_robust <- function(coord_string) {
  coord_string <- as.character(coord_string)
  
  if (isTRUE(is.na(coord_string)) || isTRUE(coord_string == "") || isTRUE(coord_string == "NA")) {
    return(numeric(0))
  }
  
  tryCatch({
    # Handle numpy format: [np.int64(394), np.int64(382), ...]
    if (grepl("np\\.int64", coord_string)) {
      # Extract numbers from np.int64(...) format
      numbers <- regmatches(coord_string, gregexpr("np\\.int64\\((\\d+)\\)", coord_string))[[1]]
      coords <- as.numeric(gsub("np\\.int64\\((\\d+)\\)", "\\1", numbers))
    } else if (grepl("\\[.*\\]", coord_string)) {
      # Handle regular list format [1, 2, 3, ...]
      # Remove brackets and split by comma
      clean_string <- gsub("\\[|\\]", "", coord_string)
      coords <- as.numeric(unlist(strsplit(clean_string, ",")))
      coords <- coords[!is.na(coords)]
    } else {
      # Handle simple comma-separated format
      coords <- as.numeric(unlist(strsplit(gsub("[^0-9,.-]", "", coord_string), ",")))
      coords <- coords[!is.na(coords)]
    }
    
    # Convert from 0-based to 1-based indexing
    return(coords + 1)
    
  }, error = function(e) {
    return(numeric(0))
  })
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# covert grid cell values to lat/lon coords
lat_lon <- function(cx_index, cy_index, lon_grid, lat_grid) {
  if (cx_index < 1 || cx_index > ncol(lon_grid) || 
      cy_index < 1 || cy_index > nrow(lat_grid)) {
    return(c(NA, NA))
  }
  
  cx_int <- as.integer(round(cx_index))
  cy_int <- as.integer(round(cy_index))
  
  longitude <- lon_grid[cy_int, cx_int]
  latitude <- lat_grid[cy_int, cx_int]
  
  return(c(longitude, latitude))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Fix corrupted cx_p & cy_p values by extracting digits after "1.64"
fix_corrupted_values <- function(vals) {
  if (length(vals) == 0) return(vals)
  
  # Convert to character to work with string manipulation
  vals_char <- as.character(vals)
  
  # Extract everything after "1.64"
  fixed_vals <- numeric(length(vals))
  
  for (i in seq_along(vals_char)) {
    val_str <- vals_char[i]
    
    # Find position after "1.64"
    pos <- regexpr("1\\.64", val_str)
    if (pos > 0) {
      # Extract everything after "1.64"
      remaining <- substr(val_str, pos + 4, nchar(val_str))
      # Convert back to numeric
      fixed_vals[i] <- as.numeric(remaining)
    } else {
      # If "1.64" pattern not found, keep original value
      fixed_vals[i] <- vals[i]
    }
  }
  
  return(fixed_vals)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Filter TCs out of nor'easter dataset via genesis, min. SLP, & max. wind parameters
#   (ERA5 data structure)
TC_filter_ERA5 <- function(noreasters, lon_grid, lat_grid, min_lat = 30, 
                           min_slp = 98000, max_wind = 20) {
  
  cat("Filtering tropical cyclones from nor'easter dataset...\n")
  cat("Criteria: Genesis south of", min_lat, "°N, min pressure <=", min_slp, "Pa, max wind >=", max_wind, "m/s\n\n")
  
  filtered_noreasters <- list()
  
  # Process each year
  for (year in names(noreasters)) {
    cat(paste("Processing year:", year, "...\n"))
    
    # Get the dataframe for this year
    storm_df <- noreasters[[year]]
    
    # Check if year has data
    if (is.null(storm_df) || nrow(storm_df) == 0) {
      cat("  No storms found for this year\n")
      filtered_noreasters[[year]] <- storm_df
      next
    }
    
    original_count <- nrow(storm_df)
    
    # Extract month from storm ID (MMDDYYYY-HH-XX-XX)
    storm_id <- storm_df$ID
    storm_month <- substr(storm_id, 1, 2)
    
    # Identify Sep–Oct–Nov storms
    is_SON <- storm_month %in% c("09", "10", "11")
    n_SON <- sum(is_SON, na.rm = TRUE)
    
    # Calculate genesis latitude for each storm
    genesis_lat <- numeric(original_count)
    
    for (i in 1:original_count) {
      # Extract first position (genesis point)
      cx_p_first <- extract_coordinates_robust(storm_df[i, "cx_p"])[1]
      cy_p_first <- extract_coordinates_robust(storm_df[i, "cy_p"])[1]
      
      # Skip if coordinates are invalid
      if (is.na(cx_p_first) || is.na(cy_p_first) || 
          length(cx_p_first) == 0 || length(cy_p_first) == 0) {
        genesis_lat[i] <- NA
        next
      }
      
      # Apply corruption fix
      cx_p_first <- fix_corrupted_values(cx_p_first)[1]
      cy_p_first <- fix_corrupted_values(cy_p_first)[1]
      
      # Convert to lat/lon
      coords <- lat_lon(cx_p_first, cy_p_first, lon_grid, lat_grid)
      genesis_lat[i] <- coords[2]  # latitude is second element
    }
    
    # Get minimum pressure for each storm
    min_pressure <- storm_df$`min_p_h (Pa)`
    
    # Get maximum wind for each storm
    max_wind_per_storm <- storm_df$`max_w_all_h_750 (m/s)`
    
    # Define TC criteria (ONLY applied to Sep–Oct–Nov storms)
    is_TC <- is_SON &                          # Only SON storms evaluated
      (genesis_lat < min_lat) &                # Forms south of min_lat
      (min_pressure <= min_slp) &              # Reaches min pressure threshold
      (max_wind_per_storm >= max_wind)          # Reaches max wind threshold
    
    # Treat NAs as non-TC
    is_TC[is.na(is_TC)] <- FALSE
    
    # Filter out TCs - keep the dataframe structure
    filtered_df <- storm_df[!is_TC, ]
    
    # Reset row names
    rownames(filtered_df) <- NULL
    
    # Store filtered dataframe
    filtered_noreasters[[year]] <- filtered_df
    
    # Report statistics
    n_filtered <- sum(is_TC)
    n_kept <- sum(!is_TC)
    
    cat(paste(
      "  Year", year, ":",
      n_kept, "storms kept out of", original_count,
      "(removed", n_filtered, "TCs)\n"
    ))
    
    cat(paste(
      "    →", n_SON, "SON storms evaluated,",
      n_filtered, "removed\n"
    ))
  }
  
  cat("\nTC filtering complete!\n")
  return(filtered_noreasters)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to filter storms based on year range (ERA5 data structure)
filter_years_era5 <- function(storm_data, start_year = 1980, end_year = 2014) {
  # Create filtered version of storm_data
  filtered_data <- list()
  
  cat(paste("Filtering ERA5 data to years", start_year, "to", end_year, "\n"))
  
  # Get available years and filter to desired range
  available_years <- names(storm_data)
  years_to_keep <- available_years[as.numeric(available_years) >= start_year & 
                                     as.numeric(available_years) <= end_year]
  
  storms_kept <- 0
  storms_total <- 0
  
  # Copy data for years in the specified range
  for (year in years_to_keep) {
    if (year %in% available_years) {
      year_data <- storm_data[[year]]
      filtered_data[[year]] <- year_data
      storms_kept <- storms_kept + nrow(year_data)
      cat(paste("Year", year, ": kept", nrow(year_data), "storms\n"))
    }
  }
  
  # Count total storms in original data
  for (year in available_years) {
    storms_total <- storms_total + nrow(storm_data[[year]])
  }
  
  cat(paste("Total: kept", length(years_to_keep), "years,", 
            storms_kept, "storms out of", storms_total, "total\n"))
  
  return(filtered_data)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Filter a flat year-keyed ERA5 or CFSR list by season
# era5_data structure: era5_noreasters[["1980"]] = dataframe
filter_season_era5 <- function(era5_data, season) {
  season_months <- switch(season,
                          SON = c(9, 10, 11),
                          DJF = c(12,  1,  2),
                          MA  = c(3,   4),
                          stop("Unknown season '", season, "'. Use 'SON', 'DJF', or 'MA'.")
  )
  lapply(era5_data, function(year_df) {
    filter_df_by_months(year_df, season_months)  # reuses existing helper
  })
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Filter a single year-dataframe to storms whose genesis month is in `months`
filter_df_by_months <- function(df, months) {
  if (is.null(df) || nrow(df) == 0) return(df)
  storm_months <- sapply(df$ID, extract_month)   # <-- your existing function
  df[storm_months %in% months, ]
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to extract month from storm ID
extract_month <- function(id) {
  if (grepl("^\\d{8}-", id)) {
    # Format: MMDDYYYY-HH-XX-XX (e.g., "09141980-18-66-23")
    return(as.integer(substr(id, 1, 2)))
  } else if (grepl("^\\d{4}-\\d{2}-\\d{2}", id)) {
    # Format: YYYY-MM-DD HH:MM:SS-XXX-XXX (e.g., "1979-09-04 06:00:00-435-145")
    return(as.integer(substr(id, 6, 7)))
  } else {
    return(NA)
  }
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to merge wind data into existing CRCM5 storm data structure
merge_wind_into_storms <- function(noreasters, wind_dirs) {
  
  require(data.table)
  
  cat("Merging wind data into storm tracks...\n")
  
  # Process each model in noreasters
  for (model in names(noreasters)) {
    cat(paste("Processing model:", model, "...\n"))
    
    # Check if wind directory is provided for this model
    if (!model %in% names(wind_dirs)) {
      cat(paste("  Warning: No wind directory specified for model", model, "\n"))
      next
    }
    
    wind_dir <- wind_dirs[[model]]
    
    if (!dir.exists(wind_dir)) {
      cat(paste("  Warning: Wind directory does not exist:", wind_dir, "\n"))
      next
    }
    
    # Get all wind CSV files for this model
    wind_files <- list.files(wind_dir, pattern = "\\.csv$", 
                             full.names = TRUE, recursive = TRUE)
    
    if (length(wind_files) == 0) {
      cat(paste("  Warning: No wind CSV files found in:", wind_dir, "\n"))
      next
    }
    
    # Read all wind data for this model
    cat(paste("  Reading", length(wind_files), "wind files for", model, "...\n"))
    
    wind_data_list <- list()
    for (i in seq_along(wind_files)) {
      if (i %% 10 == 0) {
        cat(paste("    Processing file", i, "of", length(wind_files), "\n"))
      }
      tryCatch({
        wind_data_list[[i]] <- read_wind_file(wind_files[i], model_name = model)
      }, error = function(e) {
        cat(paste("    ERROR in file:", basename(wind_files[i]), "\n"))
        cat(paste("    Error message:", e$message, "\n"))
      })
    }
    
    # Remove NULL entries
    wind_data_list <- wind_data_list[!sapply(wind_data_list, is.null)]
    
    if (length(wind_data_list) == 0) {
      cat(paste("  Error: No wind files could be read successfully\n"))
      next
    }
    
    cat(paste("  Successfully read", length(wind_data_list), "files. Now combining...\n"))
    
    model_wind_data <- data.table::rbindlist(wind_data_list, fill = TRUE)
    
    if (!is.data.table(model_wind_data)) {
      model_wind_data <- as.data.table(model_wind_data)
    }
    
    cat(paste("  Combined into", nrow(model_wind_data), "wind track points\n"))
    
    # Aggregate wind data by storm (one row per storm with list columns)
    cat("  Aggregating wind data by storm...\n")
    wind_aggregated_all <- model_wind_data[, .(
      u_wind = list(u_wind),
      v_wind = list(v_wind),
      wind_speed = list(wind_speed),
      wind_direction = list(wind_direction),
      max_wind_in_box = list(max_wind_in_box),
      max_wind_x_idx = list(max_wind_x_idx),
      max_wind_y_idx = list(max_wind_y_idx),
      rlat_idx = list(rlat_idx),
      rlon_idx = list(rlon_idx),
      storm_radius_km = list(storm_radius_km),
      grid_radius_cells = list(grid_radius_cells),
      datetime = list(datetime),
      season_start_year = first(season_start_year),
      season_end_year = first(season_end_year)
    ), by = storm_id]
    
    cat(paste("  Aggregated to", nrow(wind_aggregated_all), "unique storms\n"))
    
    # Process each year
    for (year in names(noreasters[[model]])) {
      storm_df <- noreasters[[model]][[year]]
      
      if (is.null(storm_df) || nrow(storm_df) == 0) next
      
      cat(paste("  Processing year", year, "...\n"))
      cat(paste("    Original storm count:", nrow(storm_df), "\n"))
      
      # Get the storm IDs that exist in the nor'easter data for this year
      noreaster_ids <- storm_df$ID
      
      # Filter wind data to only storms that exist in nor'easter data
      year_num <- as.numeric(year)
      year_wind <- wind_aggregated_all[
        storm_id %in% noreaster_ids &
          (season_start_year == year_num | season_end_year == year_num)
      ]
      
      if (nrow(year_wind) == 0) {
        cat(paste("    Warning: No matching wind data for year", year, "\n"))
        next
      }
      
      cat(paste("    Found", nrow(year_wind), "matching storms with wind data\n"))
      
      # Convert storm_df to data.table for easier manipulation
      storm_dt <- as.data.table(storm_df)
      
      # Merge by matching storm IDs
      # Using data.table merge which is more reliable
      merged_dt <- merge(
        storm_dt,
        year_wind,
        by.x = "ID",
        by.y = "storm_id",
        all.x = TRUE,  # Keep all nor'easter storms
        all.y = FALSE  # Don't add storms that aren't in nor'easter data
      )
      
      # Verify row count hasn't changed
      if (nrow(merged_dt) != nrow(storm_df)) {
        cat(paste("    WARNING: Row count changed from", nrow(storm_df), 
                  "to", nrow(merged_dt), "\n"))
        cat(paste("    This suggests duplicate IDs. Investigating...\n"))
        
        # Check for duplicate IDs
        dup_ids <- storm_df$ID[duplicated(storm_df$ID)]
        if (length(dup_ids) > 0) {
          cat(paste("    Found", length(dup_ids), "duplicate IDs in storm data\n"))
        }
        
        dup_wind_ids <- year_wind$storm_id[duplicated(year_wind$storm_id)]
        if (length(dup_wind_ids) > 0) {
          cat(paste("    Found", length(dup_wind_ids), "duplicate IDs in wind data\n"))
        }
      }
      
      # Convert back to data.frame if original was data.frame
      if (!is.data.table(storm_df)) {
        merged_dt <- as.data.frame(merged_dt)
      }
      
      # Store back into the structure
      noreasters[[model]][[year]] <- merged_dt
      
      # Report merge statistics
      n_total <- nrow(merged_dt)
      n_matched <- sum(!is.na(merged_dt$wind_speed))
      cat(paste("    Final: Year", year, ":", n_matched, "of", n_total, 
                "storms have wind data\n"))
    }
  }
  
  cat("Wind data merge complete!\n")
  return(noreasters)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Filter TCs out of CRCM5 nor'easter dataset via genesis, min. SLP, & max. wind parameters
TC_filter <- function(noreasters, lon_grid, lat_grid, min_lat = 30, min_slp = 98000, max_wind = 20) {
  
  cat("Filtering tropical cyclones from nor'easter dataset...\n")
  cat("Criteria: Genesis south of", min_lat, "°N, min pressure <=", min_slp, "Pa, max wind >=", max_wind, "m/s\n\n")
  
  filtered_noreasters <- list()
  
  # Process each model
  for (model in names(noreasters)) {
    cat(paste("Processing model:", model, "...\n"))
    filtered_noreasters[[model]] <- list()
    
    # Process each year
    for (year in names(noreasters[[model]])) {
      storm_df <- noreasters[[model]][[year]]
      
      if (is.null(storm_df) || nrow(storm_df) == 0) {
        filtered_noreasters[[model]][[year]] <- storm_df
        next
      }
      
      original_count <- nrow(storm_df)
      
      # Extract month from storm ID (MMDDYYYY-HH-XX-XX)
      storm_id <- storm_df$ID
      storm_month <- sapply(storm_id, extract_month)
      # Identify Sep–Oct–Nov storms
      is_SON <- storm_month %in% c(9, 10, 11)
      n_SON <- sum(is_SON, na.rm = TRUE)
      
      # Calculate genesis latitude for each storm
      genesis_lat <- numeric(nrow(storm_df))
      
      for (i in 1:nrow(storm_df)) {
        # Extract first position (genesis point)
        cx_p_first <- extract_coordinates_robust(storm_df[i, "cx_p"])[1]
        cy_p_first <- extract_coordinates_robust(storm_df[i, "cy_p"])[1]
        
        # Skip if coordinates are invalid
        if (is.na(cx_p_first) || is.na(cy_p_first) || 
            length(cx_p_first) == 0 || length(cy_p_first) == 0) {
          genesis_lat[i] <- NA
          next
        }
        
        # Apply corruption fix
        cx_p_first <- fix_corrupted_values(cx_p_first)[1]
        cy_p_first <- fix_corrupted_values(cy_p_first)[1]
        
        # Convert to lat/lon
        coords <- lat_lon(cx_p_first, cy_p_first, lon_grid, lat_grid)
        genesis_lat[i] <- coords[2]  # latitude is second element
      }
      
      # Get maximum wind in box for each storm
      max_wind_per_storm <- sapply(storm_df$max_wind_in_box, function(wind_list) {
        if (is.null(wind_list) || length(wind_list) == 0) {
          return(NA)
        }
        max(wind_list, na.rm = TRUE)
      })
      
      # Get minimum pressure for each storm
      min_pressure <- storm_df$`min_p (Pa)`
      
      # Define TC criteria (ONLY applied to Sep–Oct–Nov storms)
      is_TC <- is_SON &                          # Only SON storms evaluated
        (genesis_lat < min_lat) &                # Forms south of min_lat
        (min_pressure <= min_slp) &              # Reaches min pressure threshold
        (max_wind_per_storm >= max_wind)          # Reaches max wind threshold
      
      # Treat NAs as non-TC
      is_TC[is.na(is_TC)] <- FALSE
      
      # Filter out TCs (keep only non-TCs)
      filtered_df <- storm_df[!is_TC, ]
      
      # Reset row names
      rownames(filtered_df) <- NULL
      
      # Store filtered data
      filtered_noreasters[[model]][[year]] <- filtered_df
      
      # Report statistics
      n_filtered <- sum(is_TC)
      n_kept <- sum(!is_TC)
      
      cat(paste(
        "  Year", year, ":",
        n_kept, "storms kept out of", original_count,
        "(removed", n_filtered, "TCs)\n"
      ))
      
      cat(paste(
        "    →", n_SON, "SON storms evaluated,",
        n_filtered, "removed\n"
      ))
    }
  }
  
  cat("\nTC filtering complete!\n")
  return(filtered_noreasters)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# master function for filtering CRCM5 nor'easters
filter_crcm5_noreasters <- function(noreasters, lat_grid, lon_grid,
                                  wind_dirs = NULL,
                                  apply_region_filter = FALSE,
                                  apply_TC_filter = FALSE,
                                  TC_min_lat = 30,
                                  TC_min_slp = 98000,
                                  TC_max_wind = 20,
                                  lon_min = -90, lon_max = -40,
                                  lat_min = 20, lat_max = 50) {
  filtered_data <- list()
  
  for (model in names(noreasters)) {
    model_data <- noreasters[[model]]
    filtered_model_data <- list()
    
    for (member_or_year in names(model_data)) {
      member_data <- model_data[[member_or_year]]
      filtered_member_data <- list()
      
      if (is.list(member_data) && !is.data.frame(member_data)) {
        # This should be ensemble data with years as keys
        for (year in names(member_data)) {
          df <- member_data[[year]]
          
          # CRITICAL: Handle the case where data is stored as separate vectors
          if (!is.data.frame(df)) {
            # Convert list of vectors to data frame
            if (is.list(df) && length(df) > 0) {
              # Check if we have named elements that look like columns
              if (!is.null(names(df)) && any(c("ID", "cx_p", "cy_p") %in% names(df))) {
                # Convert to data frame
                max_length <- max(sapply(df, length))
                df_list <- lapply(df, function(x) {
                  if (length(x) < max_length) {
                    c(x, rep(NA, max_length - length(x)))
                  } else {
                    x[1:max_length]
                  }
                })
                df <- as.data.frame(df_list, stringsAsFactors = FALSE)
              } else {
                next  # Skip if we can't convert to proper format
              }
            } else {
              next  # Skip if not valid data
            }
          }
          
          # Safely skip if not a valid data frame or empty
          if (!is.data.frame(df) || nrow(df) == 0) next
          
          # Apply filters with updated functions
          df <- df |>
            filter_USEC(lat_grid, lon_grid) |>
            filter_USECPassedOnly() |>
            filter_time() |>
            filter_dist()
          
          if (nrow(df) > 0) {
            filtered_member_data[[year]] <- df
          }
        }
        
        if (length(filtered_member_data) > 0) {
          filtered_model_data[[member_or_year]] <- filtered_member_data
        }
      } else {
        # This should be single-member data
        df <- member_data
        
        # CRITICAL: Handle the case where data is stored as separate vectors
        if (!is.data.frame(df)) {
          if (is.list(df) && length(df) > 0) {
            # Check if we have named elements that look like columns
            if (!is.null(names(df)) && any(c("ID", "cx_p", "cy_p") %in% names(df))) {
              # Convert to data frame
              max_length <- max(sapply(df, length))
              df_list <- lapply(df, function(x) {
                if (length(x) < max_length) {
                  c(x, rep(NA, max_length - length(x)))
                } else {
                  x[1:max_length]
                }
              })
              df <- as.data.frame(df_list, stringsAsFactors = FALSE)
            } else {
              next  # Skip if we can't convert to proper format
            }
          } else {
            next  # Skip if not valid data
          }
        }
        
        if (!is.data.frame(df) || nrow(df) == 0) next
        
        # Apply filters with updated functions
        df <- df |>
          filter_USEC(lat_grid, lon_grid) |>
          filter_USECPassedOnly() |>
          filter_time() |>
          filter_dist()
        
        if (nrow(df) > 0) {
          filtered_model_data[[member_or_year]] <- df
        }
      }
    }
    
    if (length(filtered_model_data) > 0) {
      filtered_data[[model]] <- filtered_model_data
    }
  }
  
  # Step: Merge wind data (if wind directories provided)
  # This happens AFTER USEC/time/distance filters but BEFORE TC and genesis filters
  if (!is.null(wind_dirs)) {
    cat("\nMerging wind data...\n")
    filtered_data <- merge_wind_into_storms(filtered_data, wind_dirs)
  }
  
  # Step: Apply TC filter (if requested and wind data available)
  if (apply_TC_filter) {
    if (!is.null(wind_dirs)) {
      cat("\nApplying TC filter...\n")
      filtered_data <- TC_filter(filtered_data, lon_grid, lat_grid,
                                 min_lat = TC_min_lat,
                                 min_slp = TC_min_slp, 
                                 max_wind = TC_max_wind)
    } else {
      cat("\nWarning: Cannot apply TC filter - wind data not available\n")
      cat("Please provide wind_dirs parameter to enable TC filtering\n")
    }
  }
  
  # Step: Apply analysis region filter (genesis filter)
  if (apply_region_filter) {
    cat("\nApplying genesis region filter...\n")
    
    for (model in names(filtered_data)) {
      model_data <- filtered_data[[model]]
      
      for (member_or_year in names(model_data)) {
        member_data <- model_data[[member_or_year]]
        
        if (is.list(member_data) && !is.data.frame(member_data)) {
          # Ensemble data
          for (year in names(member_data)) {
            df <- member_data[[year]]
            if (is.data.frame(df) && nrow(df) > 0) {
              df <- filter_analysis_region(df, lat_grid, lon_grid,
                                           lon_min, lon_max, lat_min, lat_max)
              if (nrow(df) > 0) {
                filtered_data[[model]][[member_or_year]][[year]] <- df
              } else {
                filtered_data[[model]][[member_or_year]][[year]] <- NULL
              }
            }
          }
        } else {
          # Single-member data
          df <- member_data
          if (is.data.frame(df) && nrow(df) > 0) {
            df <- filter_analysis_region(df, lat_grid, lon_grid,
                                         lon_min, lon_max, lat_min, lat_max)
            if (nrow(df) > 0) {
              filtered_data[[model]][[member_or_year]] <- df
            } else {
              filtered_data[[model]][[member_or_year]] <- NULL
            }
          }
        }
      }
    }
  }
  
  return(filtered_data)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to filter storms based on year range for multi-model CRCM5 data structure
filter_years_multi_model <- function(storm_data, start_year = 1980, end_year = 2014) {
  # Create filtered version of storm_data
  filtered_data <- list()
  
  cat(paste("Filtering data to years", start_year, "to", end_year, "\n"))
  
  # Process each model
  for (model in names(storm_data)) {
    cat(paste("Processing model:", model, "\n"))
    
    # Handle single-member models (canesm5, cnrm, noresm2)
    filtered_data[[model]] <- list()
    
    # Filter years for this model
    available_years <- names(storm_data[[model]])
    years_to_keep <- available_years[as.numeric(available_years) >= start_year & 
                                       as.numeric(available_years) <= end_year]
    
    storms_kept <- 0
    storms_total <- 0
    
    for (year in years_to_keep) {
      if (year %in% available_years) {
        year_data <- storm_data[[model]][[year]]
        filtered_data[[model]][[year]] <- year_data
        storms_kept <- storms_kept + nrow(year_data)
      }
    }
    
    # Count total storms in original data for this model
    for (year in available_years) {
      storms_total <- storms_total + nrow(storm_data[[model]][[year]])
    }
    
    cat(paste("  Model", model, ": kept", length(years_to_keep), "years,", 
              storms_kept, "storms out of", storms_total, "total\n"))
  }
  
  return(filtered_data)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Filter an entire multi-model, multi-year list structure by season:
#     "SON" | "DJF" | "MA"
filter_season_multi_model <- function(noreasters_list, season) {
  season_months <- switch(season,
                          SON = c(9, 10, 11),
                          DJF = c(12,  1,  2),
                          MA  = c(3,   4),
                          stop("Unknown season '", season, "'. Use 'SON', 'DJF', or 'MA'.")
  )
  lapply(noreasters_list, function(model_data) {
    lapply(model_data, function(year_df) {
      filter_df_by_months(year_df, season_months)
    })
  })
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function for extracting intensity values from ERA5 & CFSR nor'easters
process_era5_for_intensity <- function(model_data, metric) {
  model_results <- list()
  for (year in names(model_data)) {  # Use names() to get the year names
    storm_tracks <- model_data[[year]]
    year_storms <- list()
    
    # Process each storm individually
    for (i in 1:nrow(storm_tracks)) {
      storm_data <- list()
      
      if (metric == "min_slp") {
        if ("min_p_h (Pa)" %in% names(storm_tracks)) {
          min_slp <- storm_tracks[i, "min_p_h (Pa)"]
        } else if ("min_p (Pa)" %in% names(storm_tracks)) {
          min_slp <- storm_tracks[i, "min_p (Pa)"]
        } else if ("min_p" %in% names(storm_tracks)) {
          min_slp <- storm_tracks[[i, "min_p_h (Pa)"]]
        } else {
          min_slp <- NA
        }
        storm_data[["min_slp"]] <- min_slp
      } else if (metric == "max_wind") {
        # Extract max wind - check for ERA5 variable first, then fallback
        if ("max_w_all_h_750 (m/s)" %in% names(storm_tracks)) {
          max_wind <- storm_tracks[i, "max_w_all_h_750 (m/s)"]
        } else if ("max_wind_in_box" %in% names(storm_tracks)) {
          max_wind <- max(unlist(storm_tracks[i, "max_wind_in_box"]), na.rm = TRUE)
        } else {
          max_wind <- NA
        }
        storm_data[["max_wind"]] <- max_wind
      }
      
      # Store individual storm track
      year_storms[[paste0("storm_", i)]] <- storm_data
    }
    
    model_results[[year]] <- year_storms
  }
  
  return(model_results)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function for extracting specified intensity values from CRCM5 nor'easters:
extract_intensity_values <- function(storm_data_1980_2014, metric) {
  results <- list()
  
  single_member_models <- c("canesm5", "cnrm", "noresm2", "mpi", "era5")
  #multi_member_models <- c("era5", "mpi")  # Update if needed
  
  for (model in single_member_models) {
    if (model %in% names(storm_data_1980_2014)) {
      cat(paste("Processing", model, "...\n"))
      results[[model]] <- list()
      
      for (year_str in names(storm_data_1980_2014[[model]])) {
        year <- as.numeric(year_str)
        #if (year >= 1951 && year <= 2014) {
        year_data <- storm_data_1980_2014[[model]][[year_str]]
        
        if (metric == "min_slp") {
          if (is.data.frame(year_data) && nrow(year_data) > 0) {
            year_storms <- list()
            for (i in 1:nrow(year_data)) {
              min_slp <- year_data[i, "min_p (Pa)"]
              year_storms[[paste0("storm_", i)]] <- list(min_slp = year_data[i, "min_p (Pa)"])
            }
            results[[model]][[year_str]] <- year_storms
          }
        }
        if (metric == "max_wind") {
          if (is.data.frame(year_data) && nrow(year_data) > 0) {
            year_storms <- list()
            for (i in 1:nrow(year_data)) {
              max_wind <- max(year_data[[i, "max_wind_in_box"]])
              year_storms[[paste0("storm_", i)]] <- list(max_wind = max_wind)
            }
            results[[model]][[year_str]] <- year_storms
          }
        }
      }
    }
  }
  
  return(results)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate CFSR intensity quantiles for a given time period
calculate_cfsr_intensity_quantiles <- function(model_data, metric, probs = my.probs) {
  all_intensity_values <- numeric()
  
  for (year in names(model_data)) {
    year_storms <- model_data[[year]]
    
    for (storm_name in names(year_storms)) {
      storm <- year_storms[[storm_name]]
      
      if (metric == "min_slp") {
        val <- storm[["min_slp"]][["min_p (Pa)"]]
      } else if (metric == "max_wind") {
        val <- storm[["max_wind"]]
      }
      
      # Handle data.table/data.frame just in case
      if (inherits(val, "data.table") || inherits(val, "data.frame")) {
        val <- as.numeric(val[[1]])
      }
      
      if (!is.null(val) && length(val) > 0 && !is.na(val)) {
        all_intensity_values <- c(all_intensity_values, val)
      }
    }
  }
  
  print(paste("Model: CFSR"))
  print(paste("Metric:", metric))
  print(paste("Length:", length(all_intensity_values)))
  print(paste("Class of values:", class(all_intensity_values)))
  print(head(all_intensity_values))
  
  if (length(all_intensity_values) > 0) {
    quantile_results <- quantile(all_intensity_values, probs = probs, na.rm = TRUE)
  } else {
    quantile_results <- NA
  }
  
  return(quantile_results)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to calculate ERA5 intensity quantiles for a given time period
calculate_era5_intensity_quantiles <- function(model_data, metric, probs = my.probs) {
  # Initialize an empty numeric vector to collect all intensity values
  all_intensity_values <- numeric()  
  
  # Process the raw ERA5 data for all years and storms
  for (year in names(model_data)) {
    year_storms <- model_data[[year]]
    
    # Process each storm individually
    for (storm_name in names(year_storms)) {
      storm <- year_storms[[storm_name]]
      
      # Access the appropriate metric value for each storm
      if (metric == "min_slp") {
        val <- storm[["min_slp"]]
        if (is.null(val)) val <- storm[["min_p (Pa)"]]  # fallback
      } else if (metric == "max_wind") {
        val <- storm[["max_wind"]]                        # new field name
        if (is.null(val)) val <- storm[["max_w_all_h_750 (m/s)"]]  # fallback to old
        if (is.null(val)) val <- storm[["max_wind_in_box"]]         # original fallback
      }
      
      # Check if val is a data.table or data.frame and extract the numeric value
      if (inherits(val, "data.table") || inherits(val, "data.frame")) {
        val <- as.numeric(val[[1]])  # Extract the first column as numeric
      }
      
      # Check if the value is valid (not NULL or NA)
      if (!is.null(val) && !is.na(val)) {
        all_intensity_values <- c(all_intensity_values, val)  # Append the value
      }
    }
  }
  
  # Output check
  print(paste("Model: ERA5"))
  print(paste("Metric:", metric))
  print(paste("Length:", length(all_intensity_values)))
  print(paste("Class of values:", class(all_intensity_values)))
  print(head(all_intensity_values))
  
  # Calculate the quantiles for all collected intensity values
  if (length(all_intensity_values) > 0) {
    quantile_results <- quantile(all_intensity_values, probs = probs, na.rm = TRUE)
  } else {
    quantile_results <- NA  # If no valid data, return NA
  }
  
  return(quantile_results)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to calculate CRCM5-CMIP6 ensemble intensity quantiles for a given time period
calculate_ensemble_intensity_quantiles <- function(results, metric, probs = my.probs) {
  all_values <- numeric()
  
  # Loop through all models and collect all values
  for (model in c("canesm5", "cnrm", "noresm2", "mpi")) {
    model_data <- results[[model]]
    
    # Loop through years
    for (year in names(model_data)) {
      storms <- model_data[[year]]
      if (is.null(storms) || length(storms) == 0) next
      
      # Loop through storms
      for (storm in storms) {
        if (metric == "min_slp") {
          val <- storm[["min_slp"]]
          if (is.null(val)) val <- storm[["min_p (Pa)"]]
        } else if (metric == "max_wind") {
          val <- storm[["max_wind"]]
        }
        
        if (!is.null(val) && !is.na(val)) {
          all_values <- c(all_values, as.numeric(val))
        }
      }
    }
  }
  
  # Calculate quantiles from combined values
  if (length(all_values) > 0) {
    return(quantile(all_values, probs = probs, na.rm = TRUE))
  } else {
    return(NA)
  }
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate (& plot) specific quantiles over time for ERA5 & CFSR with ONE-TAILED tests
#     (H1 = increasing intensity, either decreasing min. SLP or increasing max. wind)
plot_era5_intensity_quantile_trends <- function(intensity_data, metric, 
                                                quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9), 
                                                file_path = "file_name") {
  all_years <- sort(as.numeric(names(intensity_data)))
  
  pdf(file_path, family = "Helvetica-Narrow", width = 6, height = 4)
  par(mar = c(4, 4, 2, 1), oma = c(2, 2, 2, 1))
  
  yearly_quants <- matrix(NA, nrow = length(all_years), ncol = length(quantiles))
  rownames(yearly_quants) <- all_years
  colnames(yearly_quants) <- paste0("Q", quantiles)
  
  for (i in seq_along(all_years)) {
    year_str <- as.character(all_years[i])
    if (year_str %in% names(intensity_data)) {
      year_storms <- intensity_data[[year_str]]
      values <- c()
      
      for (storm_name in names(year_storms)) {
        storm <- year_storms[[storm_name]]
        
        if (metric == "min_slp") {
          val <- storm[["min_slp"]]
        } else if (metric == "max_wind") {
          val <- storm[["max_wind"]]
        }
        
        # FIXED: guard against length-zero before is.finite()
        if (!is.null(val) && length(val) > 0) {
          if (is.data.frame(val) || is.matrix(val)) {
            val <- as.numeric(val[[1]])
          } else if (is.list(val)) {
            val <- as.numeric(unlist(val))
          } else {
            val <- as.numeric(val)
          }
          
          val <- val[is.finite(val)]
          if (length(val) > 0) {
            values <- c(values, val)
          }
        }
      }
      
      if (length(values) > 0) {
        yearly_quants[i, ] <- quantile(values, probs = quantiles, na.rm = TRUE)
      }
    }
  }
  
  y_label <- ifelse(metric == "min_slp", "Minimum Pressure (Pa)", "Maximum Wind Speed (m/s)")
  
  plot(all_years, yearly_quants[, 1], type = "n",
       ylim = range(yearly_quants, na.rm = TRUE),
       xlab = "", ylab = "")
  
  mtext("Quantile Trends: ERA5", 3, 1, adj = 0, family = "Helvetica", font = 2, cex = 1)
  mtext(y_label, 2, 2.5, family = "Helvetica", font = 2, cex = 1)
  mtext("Year", 1, 2.5, family = "Helvetica", font = 2, cex = 1)
  
  colors <- viridis(length(quantiles))
  p_values <- numeric(length(quantiles))
  
  for (j in seq_along(quantiles)) {
    lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
    
    valid_idx <- !is.na(yearly_quants[, j])
    if (sum(valid_idx) > 2) {
      trend_data <- data.frame(
        year = all_years[valid_idx],
        quant = yearly_quants[valid_idx, j]
      )
      lm_fit <- lm(quant ~ year, data = trend_data)
      abline(lm_fit, col = colors[j], lty = 2, lwd = 1.5)
      
      # ONE-SIDED p-value (unchanged from previous version, just preserved here)
      slope        <- summary(lm_fit)$coefficients[2, 1]
      p_two_tailed <- summary(lm_fit)$coefficients[2, 4]
      if (metric == "max_wind") {
        p_values[j] <- if (slope > 0) p_two_tailed / 2 else 1 - (p_two_tailed / 2)
      } else if (metric == "min_slp") {
        p_values[j] <- if (slope < 0) p_two_tailed / 2 else 1 - (p_two_tailed / 2)
      }
    } else {
      p_values[j] <- NA
    }
  }
  
  legend_text <- paste0(quantiles * 100, "th %ile (p=", 
                        ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
  legend("topright", inset = c(-0.05, -0.225), xpd = "True",
         legend = legend_text, ncol = 2,
         col = colors, lwd = 2, bty = "n", cex = 0.7,
         title = "Quantiles (solid=data, dash=trend)")
  
  dev.off()
  par(mfrow = c(1, 1))
  
  return(list(
    quantiles = yearly_quants,
    p_values = p_values
  ))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate & plot CRCM5-CMIP6 ensemble quantiles (pooling all models together)
#     (one-tailed test, H1 = increasing intensity)
plot_ensemble_quantile_trends <- function(intensity_data, metric, 
                                          quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.97, 0.99), 
                                          file_path = NULL) {
  
  all_years <- sort(unique(unlist(lapply(intensity_data, function(m) as.numeric(names(m))))))
  
  yearly_quants <- matrix(NA, nrow = length(all_years), ncol = length(quantiles))
  rownames(yearly_quants) <- all_years
  colnames(yearly_quants) <- paste0("Q", quantiles)
  
  # For each year, pool all storms from all models
  for (i in seq_along(all_years)) {
    year_str <- as.character(all_years[i])
    values <- c()
    
    # Loop through all models
    for (model in names(intensity_data)) {
      if (year_str %in% names(intensity_data[[model]])) {
        storms <- intensity_data[[model]][[year_str]]
        
        # Extract values from each storm
        for (storm in storms) {
          if (metric == "min_slp") {
            if (!is.null(storm$min_slp)) {
              val <- storm$min_slp
              if (is.data.frame(val) || is.list(val)) {
                if ("min_p (Pa)" %in% names(val)) {
                  val <- val[["min_p (Pa)"]]
                } else if (length(val) == 1) {
                  val <- unlist(val)
                }
              }
              values <- c(values, as.numeric(val))
            }
          } else if (metric == "max_wind") {
            if (!is.null(storm$max_wind)) {
              val <- storm$max_wind
              if (is.data.frame(val) || is.list(val)) {
                if ("max_wind_in_box" %in% names(val)) {
                  val <- val[["max_wind_in_box"]]
                } else if (length(val) == 1) {
                  val <- unlist(val)
                }
              }
              values <- c(values, as.numeric(val))
            }
          } else if (metric == "intensification") {
            if (!is.null(storm$max_intensification)) {
              values <- c(values, storm$max_intensification)
            }
          }
        }
      }
    }
    
    # Calculate quantiles for this year across all models
    values <- values[is.finite(values)]
    if (length(values) > 0) {
      yearly_quants[i, ] <- quantile(values, probs = quantiles)
    }
  }
  
  # Plot
  if (!is.null(file_path)) {
    pdf(file_path, family = "Helvetica-Narrow", width = 10, height = 6)
  }
  
  ylabel <- if (metric == "min_slp") "Minimum Pressure (Pa)" else if (metric == "max_wind") "Maximum Wind Speed (m/s)" else "Max 24hr Pressure Drop (Pa)"
  
  plot(all_years, yearly_quants[, 1], type = "n",
       ylim = range(yearly_quants, na.rm = TRUE),
       main = "",
       xlab = "", ylab = "")
  
  colors <- viridis(length(quantiles))
  p_values <- numeric(length(quantiles))
  
  for (j in seq_along(quantiles)) {
    # Plot quantile line
    lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
    
    # Add trend line
    valid_idx <- !is.na(yearly_quants[, j])
    if (sum(valid_idx) > 2) {
      trend_data <- data.frame(
        year = all_years[valid_idx],
        quant = yearly_quants[valid_idx, j]
      )
      lm_fit <- lm(quant ~ year, data = trend_data)
      
      # Draw trend line (dashed)
      abline(lm_fit, col = colors[j], lty = 2, lwd = 1.5)
      
      # One-sided p-value: H1 is that intensity increases over time.
      # For min_slp, increasing intensity = decreasing pressure (slope < 0).
      # For max_wind, increasing intensity = increasing wind speed (slope > 0).
      slope <- summary(lm_fit)$coefficients[2, 1]
      p_two_sided <- summary(lm_fit)$coefficients[2, 4]
      if (metric == "min_slp") {
        p_values[j] <- if (slope < 0) p_two_sided / 2 else 1 - p_two_sided / 2
      } else {
        p_values[j] <- if (slope > 0) p_two_sided / 2 else 1 - p_two_sided / 2
      }
    } else {
      p_values[j] <- NA
    }
  }
  
  # Create legend with p-values
  legend_text <- paste0(quantiles * 100, "th %ile (p=", 
                        ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
  
  legend("topright", inset = c(0, -0.175), xpd = "True",
         legend = legend_text, ncol = 2,
         col = colors, lwd = 2, bty = "n", cex = 0.7,
         title = "Quantiles (solid=data, dash=trend)")
  
  if (!is.null(file_path)) {
    dev.off()
  }
  
  return(list(quantiles = yearly_quants, p_values = p_values))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate & plot specific CRCM5 quantiles over time BY MODEL (not pooled ensemble)
#     (one-tailed test, H1 = increasing intensity)
plot_intensity_quantile_trends <- function(intensity_data, metric, quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9), file_path = "file_name") {
  all_years <- sort(unique(unlist(lapply(intensity_data, function(m) as.numeric(names(m))))))
  
  # Initialize lists to collect results for all models
  all_model_quants <- list()
  all_model_pvals <- list()
  
  pdf(file_path, family = "Helvetica-Narrow", width = 6, height = 12)
  par(mfrow = c(4, 1), mar = c(4, 4, 2, 1), oma = c(2, 2, 2, 1))
  
  for (model in names(intensity_data)) {
    yearly_quants <- matrix(NA, nrow = length(all_years), ncol = length(quantiles))
    rownames(yearly_quants) <- all_years
    colnames(yearly_quants) <- paste0("Q", quantiles)
    
    for (i in seq_along(all_years)) {
      year_str <- as.character(all_years[i])
      if (year_str %in% names(intensity_data[[model]])) {
        storms <- intensity_data[[model]][[year_str]]
        values <- c()
        
        for (storm in storms) {
          if (metric == "min_slp") {
            if (!is.null(storm$min_slp)) {
              val <- storm$min_slp
              if (is.data.frame(val) || is.list(val)) {
                if ("min_p (Pa)" %in% names(val)) {
                  val <- val[["min_p (Pa)"]]
                } else if (length(val) == 1) {
                  val <- unlist(val)
                }
              }
              values <- c(values, as.numeric(val))
            }
          } else if (metric == "max_wind") {
            if (!is.null(storm$max_wind)) {
              val <- storm$max_wind
              if (is.data.frame(val) || is.list(val)) {
                if ("max_wind_in_box" %in% names(val)) {
                  val <- val[["max_wind_in_box"]]
                } else if (length(val) == 1) {
                  val <- unlist(val)
                }
              }
              values <- c(values, as.numeric(val))
            }
          }
        }
        
        values <- values[is.finite(values)]
        if (length(values) > 0) {
          yearly_quants[i, ] <- quantile(values, probs = quantiles)
        }
      }
    }
    
    # Save this model's quantiles
    all_model_quants[[model]] <- yearly_quants
    
    plot(all_years, yearly_quants[, 1], type = "n",
         ylim = range(yearly_quants, na.rm = TRUE),
         xlab = "", ylab = "")
    
    mtext(paste("Quantile Trends:", model), 3, 1, adj = 0, family = "Helvetica", font = 2, cex = 1)
    mtext("Minimum Pressure (Pa)", 2, 2.5, family = "Helvetica", font = 2, cex = 1)
    mtext("Year", 1, 2.5, family = "Helvetica", font = 2, cex = 1)
    
    colors <- viridis(length(quantiles))
    p_values <- numeric(length(quantiles))
    
    for (j in seq_along(quantiles)) {
      lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
      
      valid_idx <- !is.na(yearly_quants[, j])
      if (sum(valid_idx) > 2) {
        trend_data <- data.frame(
          year = all_years[valid_idx],
          quant = yearly_quants[valid_idx, j]
        )
        lm_fit <- lm(quant ~ year, data = trend_data)
        abline(lm_fit, col = colors[j], lty = 2, lwd = 1.5)
        
        slope        <- summary(lm_fit)$coefficients[2, 1]
        p_two_sided  <- summary(lm_fit)$coefficients[2, 4]
        if (metric == "max_wind") {
          p_values[j] <- if (slope > 0) p_two_sided / 2 else 1 - p_two_sided / 2
        } else {
          p_values[j] <- if (slope < 0) p_two_sided / 2 else 1 - p_two_sided / 2
        }
      } else {
        p_values[j] <- NA
      }
    }
    
    # Save this model's p-values
    all_model_pvals[[model]] <- p_values
    
    legend_text <- paste0(quantiles * 100, "th %ile (p=", 
                          ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
    legend("topright", inset = c(-0.05, -0.225), xpd = "True",
           legend = legend_text, ncol = 2,
           col = colors, lwd = 2, bty = "n", cex = 0.7,
           title = "Quantiles (solid=data, dash=trend)")
  }
  
  dev.off()
  par(mfrow = c(1, 1))
  
  return(list(
    yearly_quants = all_model_quants,
    p_values      = all_model_pvals
  ))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate maximum 24-hour intensification for CFSR & ERA5 data
calculate_era5_max_intensification_rate <- function(storm_data, metric = "min_slp", 
                                                    time_window_hours = 24,
                                                    timestep_hours = 1) {
  intensification_rates <- list()
  
  # sample_interval: how many native timesteps to skip to get 6-hourly resolution
  # timestep_hours=1 (ERA5 hourly)    -> sample every 6th value
  # timestep_hours=6 (CFSR 6-hourly) -> sample every 1st value (no subsampling)
  sample_interval <- max(1, as.integer(6 / timestep_hours))
  lag_timesteps_6hourly <- time_window_hours / 6
  
  if (lag_timesteps_6hourly != round(lag_timesteps_6hourly)) {
    stop("time_window_hours must be a multiple of 6 (e.g., 6, 12, 18, 24, 30, etc.)")
  }
  lag_timesteps_6hourly <- as.integer(lag_timesteps_6hourly)
  
  cat(paste("Processing - Time window:", time_window_hours, "hours (sampled at 6-hourly resolution)...\n"))
  
  for (year_str in names(storm_data)) {
    year <- as.numeric(year_str)
    year_data <- storm_data[[year_str]]
    
    if (!is.null(year_data) && nrow(year_data) > 0) {
      year_rates <- list()
      
      for (i in 1:nrow(year_data)) {
        storm_id <- year_data[i, "ID"]
        
        cx_coords <- extract_coordinates_robust(year_data[[i, "cx_p"]])
        cy_coords <- extract_coordinates_robust(year_data[[i, "cy_p"]])
        
        if (metric == "min_slp") {
          if ("c_p_h (Pa)" %in% names(year_data)) {
            values_hourly <- extract_numpy_values(year_data[[i, "c_p_h (Pa)"]])
          } else if ("c_p (Pa)" %in% names(year_data)) {
            values_hourly <- extract_numpy_values(year_data[[i, "c_p (Pa)"]])
          } else {
            values_hourly <- NULL
          }
          invert <- TRUE
        } else if (metric == "max_wind") {
          if ("max_w_h_750 (m/s)" %in% names(year_data)) {
            values_hourly <- extract_numpy_values(year_data[[i, "max_w_h_750 (m/s)"]])
          } else if ("max_w_h (m/s)" %in% names(year_data)) {
            values_hourly <- extract_numpy_values(year_data[[i, "max_w_h (m/s)"]])
          } else if ("max_wind_in_box" %in% names(year_data)) {
            values_hourly <- unlist(year_data[[i, "max_wind_in_box"]])
          } else {
            values_hourly <- NULL
          }
          invert <- FALSE
        }
        
        if (!is.null(values_hourly) && length(values_hourly) > 0) {
          six_hourly_indices <- seq(1, length(values_hourly), by = sample_interval)
          values <- values_hourly[six_hourly_indices]
          
          if (length(cx_coords) == length(values_hourly)) {
            cx_coords <- cx_coords[six_hourly_indices]
          }
          if (length(cy_coords) == length(values_hourly)) {
            cy_coords <- cy_coords[six_hourly_indices]
          }
        } else {
          values <- NULL
        }
        
        if (!is.null(values) && length(values) >= (lag_timesteps_6hourly + 1)) {
          changes <- diff(values, lag = lag_timesteps_6hourly)
          
          if (invert) {
            intensification_changes <- changes[changes < 0]
            
            if (length(intensification_changes) > 0) {
              max_intensification <- min(intensification_changes)
              max_intensification_value <- -max_intensification
              max_intens_idx <- which.min(changes) + lag_timesteps_6hourly
            } else {
              max_intensification_value <- 0
              max_intens_idx <- NA
            }
            
            weakening_changes <- changes[changes > 0]
            if (length(weakening_changes) > 0) {
              max_weakening <- max(weakening_changes)
            } else {
              max_weakening <- 0
            }
            
          } else {
            intensification_changes <- changes[changes > 0]
            
            if (length(intensification_changes) > 0) {
              max_intensification <- max(intensification_changes)
              max_intensification_value <- max_intensification
              max_intens_idx <- which.max(changes) + lag_timesteps_6hourly
            } else {
              max_intensification_value <- 0
              max_intens_idx <- NA
            }
            
            weakening_changes <- changes[changes < 0]
            if (length(weakening_changes) > 0) {
              max_weakening <- -min(weakening_changes)
            } else {
              max_weakening <- 0
            }
          }
          
          if (!is.na(max_intens_idx)) {
            value_at_peak_intens <- values[max_intens_idx]
            
            if (max_intens_idx <= length(cx_coords) && max_intens_idx <= length(cy_coords)) {
              cx_at_max_intens <- cx_coords[max_intens_idx]
              cy_at_max_intens <- cy_coords[max_intens_idx]
            } else {
              cx_at_max_intens <- NA
              cy_at_max_intens <- NA
            }
          } else {
            value_at_peak_intens <- NA
            cx_at_max_intens <- NA
            cy_at_max_intens <- NA
          }
          
          if (metric == "min_slp") {
            peak_value <- min(values, na.rm = TRUE)
          } else {
            peak_value <- max(values, na.rm = TRUE)
          }
          
          year_rates[[paste0("storm_", i)]] <- list(
            storm_id = storm_id,
            time_window_hours = time_window_hours,
            max_intensification = max_intensification_value,
            max_weakening = max_weakening,
            value_at_max_intens = value_at_peak_intens,
            peak_value = peak_value,
            n_timesteps = length(values),
            cx_at_max_intens = cx_at_max_intens,
            cy_at_max_intens = cy_at_max_intens
          )
        }
      }
      
      if (length(year_rates) > 0) {
        intensification_rates[[year_str]] <- year_rates
      }
    }
  }
  
  return(intensification_rates)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Extract pressure/wind values from numpy string format (not specific to any particular variable)
extract_numpy_values <- function(value_string) {
  value_string <- as.character(value_string)
  
  if (isTRUE(is.na(value_string)) || isTRUE(value_string == "") || isTRUE(value_string == "NA")) {
    return(numeric(0))
  }
  
  tryCatch({
    # Handle numpy float32 format: [np.float32(100608.0), np.float32(100544.0), ...]
    if (grepl("np\\.float32", value_string)) {
      # Extract numbers from np.float32(...) format
      numbers <- regmatches(value_string, gregexpr("np\\.float32\\(([0-9.eE+-]+)\\)", value_string))[[1]]
      values <- as.numeric(gsub("np\\.float32\\(([0-9.eE+-]+)\\)", "\\1", numbers))
    } else if (grepl("\\[.*\\]", value_string)) {
      # Handle regular list format - could be comma or space separated
      clean_string <- gsub("\\[|\\]", "", value_string)
      
      # Check if comma-separated or space-separated
      if (grepl(",", clean_string)) {
        # Comma-separated: [1.0, 2.0, 3.0, ...]
        values <- as.numeric(unlist(strsplit(clean_string, ",")))
      } else {
        # Space-separated: [1.0 2.0 3.0 ...] (ERA5 format)
        values <- as.numeric(unlist(strsplit(clean_string, "\\s+")))
      }
      values <- values[!is.na(values)]
    } else {
      # Handle simple comma-separated format
      values <- as.numeric(unlist(strsplit(gsub("[^0-9,.-]", "", value_string), ",")))
      values <- values[!is.na(values)]
    }
    
    return(values)
    
  }, error = function(e) {
    return(numeric(0))
  })
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate maximum 24-hour intensification for CRCM5 data structure
calculate_max_intensification_rate <- function(storm_data, metric = "min_slp", time_window_hours = 24) {
  intensification_rates <- list()
  
  # Convert time window to number of timesteps (assuming 6-hourly data)
  timesteps_per_6hr <- 1
  lag_timesteps <- time_window_hours / 6
  
  if (lag_timesteps != round(lag_timesteps)) {
    stop("time_window_hours must be a multiple of 6 (e.g., 6, 12, 18, 24, 30, etc.)")
  }
  lag_timesteps <- as.integer(lag_timesteps)
  
  for (model in names(storm_data)) {
    cat(paste("Processing", model, "- Time window:", time_window_hours, "hours...\n"))
    intensification_rates[[model]] <- list()
    
    for (year_str in names(storm_data[[model]])) {
      year <- as.numeric(year_str)
      year_data <- storm_data[[model]][[year_str]]
      
      if (!is.null(year_data) && nrow(year_data) > 0) {
        year_rates <- list()
        
        for (i in 1:nrow(year_data)) {
          storm_id <- year_data[i, "ID"]
          
          # Extract coordinates
          cx_coords <- extract_coordinates_robust(year_data[[i, "cx_p"]])
          cy_coords <- extract_coordinates_robust(year_data[[i, "cy_p"]])
          
          if (metric == "min_slp") {
            # Extract and parse the pressure string
            values <- extract_numpy_values(year_data[[i, "c_p (Pa)"]])
            invert <- TRUE
          } else if (metric == "max_wind") {
            # Wind data is already numeric, not a string
            values <- year_data[[i, "max_wind_in_box"]]
            # Ensure it's a proper numeric vector
            if (!is.numeric(values)) {
              values <- as.numeric(values)
            }
            invert <- FALSE
          }
          
          # Need at least (lag_timesteps + 1) values to calculate the rate
          if (length(values) >= (lag_timesteps + 1)) {
            # Calculate changes over specified time window
            changes <- diff(values, lag = lag_timesteps)
            
            if (invert) {
              # For pressure: max intensification = most negative change (biggest drop)
              # Negative changes = pressure dropping = intensification
              intensification_changes <- changes[changes < 0]  # Only pressure drops
              
              if (length(intensification_changes) > 0) {
                max_intensification <- min(intensification_changes)  # most negative
                max_intensification_value <- -max_intensification  # convert to positive drop
                max_intens_idx <- which.min(changes) + lag_timesteps
              } else {
                # Storm never intensified
                max_intensification_value <- 0
                max_intens_idx <- NA
              }
              
              # For weakening (pressure rising)
              weakening_changes <- changes[changes > 0]  # Only pressure rises
              if (length(weakening_changes) > 0) {
                max_weakening <- max(weakening_changes)
              } else {
                max_weakening <- 0
              }
              
            } else {
              # For wind: max intensification = most positive change (biggest increase)
              # Positive changes = wind increasing = intensification
              intensification_changes <- changes[changes > 0]  # Only wind increases
              
              if (length(intensification_changes) > 0) {
                max_intensification <- max(intensification_changes)
                max_intensification_value <- max_intensification
                max_intens_idx <- which.max(changes) + lag_timesteps
              } else {
                # Storm never intensified
                max_intensification_value <- 0
                max_intens_idx <- NA
              }
              
              # For weakening (wind decreasing)
              weakening_changes <- changes[changes < 0]  # Only wind decreases
              if (length(weakening_changes) > 0) {
                max_weakening <- -min(weakening_changes)  # convert to positive
              } else {
                max_weakening <- 0
              }
            }
            
            # Find location at max intensification
            if (!is.na(max_intens_idx)) {
              value_at_peak_intens <- values[max_intens_idx]
              
              # Get location at max intensification
              if (max_intens_idx <= length(cx_coords) && max_intens_idx <= length(cy_coords)) {
                cx_at_max_intens <- cx_coords[max_intens_idx]
                cy_at_max_intens <- cy_coords[max_intens_idx]
              } else {
                cx_at_max_intens <- NA
                cy_at_max_intens <- NA
              }
            } else {
              value_at_peak_intens <- NA
              cx_at_max_intens <- NA
              cy_at_max_intens <- NA
            }
            
            # Get peak value for storm
            if (metric == "min_slp") {
              peak_value <- min(values, na.rm = TRUE)
            } else {
              peak_value <- max(values, na.rm = TRUE)
            }
            
            year_rates[[paste0("storm_", i)]] <- list(
              storm_id = storm_id,
              time_window_hours = time_window_hours,
              max_intensification = max_intensification_value,
              max_weakening = max_weakening,
              value_at_max_intens = value_at_peak_intens,
              peak_value = peak_value,
              n_timesteps = length(values),
              cx_at_max_intens = cx_at_max_intens,
              cy_at_max_intens = cy_at_max_intens
            )
          }
        }
        
        if (length(year_rates) > 0) {
          intensification_rates[[model]][[year_str]] <- year_rates
        }
      }
    }
  }
  
  return(intensification_rates)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate quantiles for CFSR and ERA5 intensification rates
#   (one-tailed test, H1 = increasing rates)
plot_era5_intensification_quantile_trends <- function(intensification_data, metric = "min_slp", quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.97, 0.99), file_path = NULL) {
  
  all_years <- sort(as.numeric(names(intensification_data)))
  
  yearly_quants <- matrix(NA, nrow = length(all_years), ncol = length(quantiles))
  rownames(yearly_quants) <- all_years
  colnames(yearly_quants) <- paste0("Q", quantiles)
  
  for (i in seq_along(all_years)) {
    year_str <- as.character(all_years[i])
    values <- c()
    
    if (year_str %in% names(intensification_data)) {
      storms <- intensification_data[[year_str]]
      for (storm in storms) {
        if (!is.null(storm$max_intensification) && is.finite(storm$max_intensification)) {
          values <- c(values, storm$max_intensification)
        }
      }
    }
    
    values <- values[is.finite(values)]
    if (length(values) > 0) {
      yearly_quants[i, ] <- quantile(values, probs = quantiles)
    }
  }
  
  if (!is.null(file_path)) {
    pdf(file_path, family = "Helvetica-Narrow", width = 10, height = 6)
  }
  
  if (metric == "min_slp") {
    ylabel <- "Max 24hr Pressure Drop (Pa)"
  } else if (metric == "max_wind") {
    ylabel <- "Max 24hr Wind Increase (m/s)"
  } else {
    ylabel <- "Max 24hr Intensification"
  }
  
  plot(all_years, yearly_quants[, 1], type = "n",
       ylim = range(yearly_quants, na.rm = TRUE),
       main = "ERA5 Intensification Quantile Trends",
       xlab = "Year", ylab = ylabel)
  
  colors <- viridis(length(quantiles))
  p_values <- numeric(length(quantiles))
  
  for (j in seq_along(quantiles)) {
    lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
    
    valid_idx <- !is.na(yearly_quants[, j])
    if (sum(valid_idx) > 2) {
      trend_data <- data.frame(
        year = all_years[valid_idx],
        quant = yearly_quants[valid_idx, j]
      )
      lm_fit <- lm(quant ~ year, data = trend_data)
      abline(lm_fit, col = colors[j], lty = 2, lwd = 1.5)
      
      # ONE-SIDED p-value: H1 = intensification rates increasing (positive slope)
      slope       <- summary(lm_fit)$coefficients[2, 1]
      p_two_sided <- summary(lm_fit)$coefficients[2, 4]
      p_values[j] <- if (slope > 0) p_two_sided / 2 else 1 - p_two_sided / 2
    } else {
      p_values[j] <- NA
    }
  }
  
  legend_text <- paste0(quantiles * 100, "th %ile (p=", 
                        ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
  legend("topright", inset = c(0, -0.175), xpd = TRUE,
         legend = legend_text, ncol = 2,
         col = colors, lwd = 2, bty = "n", cex = 0.7,
         title = "Quantiles (solid=data, dash=trend)")
  
  if (!is.null(file_path)) dev.off()
  
  return(list(quantiles = yearly_quants, p_values = p_values))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate CRCM5-CMIP6 ensemble quantiles for intensification rates (pooling all models together)
#     (one-tailed test, H1 = increasing rates)
plot_ensemble_intensification_quantile_trends <- function(intensification_data, quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.97, 0.99), file_path = NULL) {
  
  all_years <- sort(unique(unlist(lapply(intensification_data, function(m) as.numeric(names(m))))))
  
  yearly_quants <- matrix(NA, nrow = length(all_years), ncol = length(quantiles))
  rownames(yearly_quants) <- all_years
  colnames(yearly_quants) <- paste0("Q", quantiles)
  
  for (i in seq_along(all_years)) {
    year_str <- as.character(all_years[i])
    values <- c()
    
    for (model in names(intensification_data)) {
      if (year_str %in% names(intensification_data[[model]])) {
        storms <- intensification_data[[model]][[year_str]]
        for (storm in storms) {
          if (!is.null(storm$max_intensification) && is.finite(storm$max_intensification)) {
            values <- c(values, storm$max_intensification)
          }
        }
      }
    }
    
    values <- values[is.finite(values)]
    if (length(values) > 0) {
      yearly_quants[i, ] <- quantile(values, probs = quantiles)
    }
  }
  
  if (!is.null(file_path)) {
    pdf(file_path, family = "Helvetica-Narrow", width = 10, height = 6)
  }
  
  plot(all_years, yearly_quants[, 1], type = "n",
       ylim = range(yearly_quants, na.rm = TRUE),
       main = "", xlab = "", ylab = "")
  
  colors <- viridis(length(quantiles))
  p_values <- numeric(length(quantiles))
  
  for (j in seq_along(quantiles)) {
    lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
    
    valid_idx <- !is.na(yearly_quants[, j])
    if (sum(valid_idx) > 2) {
      trend_data <- data.frame(
        year = all_years[valid_idx],
        quant = yearly_quants[valid_idx, j]
      )
      lm_fit <- lm(quant ~ year, data = trend_data)
      abline(lm_fit, col = colors[j], lty = 2, lwd = 1.5)
      
      # ONE-SIDED p-value: H1 = intensification rates increasing (positive slope)
      slope       <- summary(lm_fit)$coefficients[2, 1]
      p_two_sided <- summary(lm_fit)$coefficients[2, 4]
      p_values[j] <- if (slope > 0) p_two_sided / 2 else 1 - p_two_sided / 2
    } else {
      p_values[j] <- NA
    }
  }
  
  legend_text <- paste0(quantiles * 100, "th %ile (p=", 
                        ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
  legend("topright", inset = c(0, -0.175), xpd = TRUE,
         legend = legend_text, ncol = 2,
         col = colors, lwd = 2, bty = "n", cex = 0.7,
         title = "Quantiles (solid=data, dash=trend)")
  
  if (!is.null(file_path)) dev.off()
  
  return(list(quantiles = yearly_quants, p_values = p_values))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Plot intensification rate quantiles over time for each CRCM5-CMIP6 pilot model
#   (one-tailed test, H1 = increasing rates)
plot_intensification_quantile_trends <- function(intensification_data, metric = "min_slp", quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9), file_path = "file_name") {
  all_years <- sort(unique(unlist(lapply(intensification_data, function(m) as.numeric(names(m))))))
  
  # Initialize lists to collect results for all models
  all_model_quants <- list()
  all_model_pvals <- list()
  
  pdf(file_path, family = "Helvetica-Narrow", width = 6, height = 12)
  par(mfrow = c(4, 1), mar = c(4, 4, 2, 1), oma = c(2, 2, 2, 1))
  
  for (model in names(intensification_data)) {
    yearly_quants <- matrix(NA, nrow = length(all_years), ncol = length(quantiles))
    rownames(yearly_quants) <- all_years
    colnames(yearly_quants) <- paste0("Q", quantiles)
    
    for (i in seq_along(all_years)) {
      year_str <- as.character(all_years[i])
      if (year_str %in% names(intensification_data[[model]])) {
        storms <- intensification_data[[model]][[year_str]]
        values <- c()
        
        for (storm in storms) {
          if (!is.null(storm$max_intensification) && is.finite(storm$max_intensification)) {
            values <- c(values, storm$max_intensification)
          }
        }
        
        values <- values[is.finite(values)]
        if (length(values) > 0) {
          yearly_quants[i, ] <- quantile(values, probs = quantiles)
        }
      }
    }
    
    # Save this model's quantiles
    all_model_quants[[model]] <- yearly_quants
    
    plot(all_years, yearly_quants[, 1], type = "n",
         ylim = range(yearly_quants, na.rm = TRUE),
         xlab = "", ylab = "")
    
    y_label <- if (metric == "min_slp") "Max 24hr Pressure Drop (Pa)" else "Max 24hr Wind Increase (m/s)"
    mtext(paste("Intensification Quantile Trends:", model), 3, 1, adj = 0, family = "Helvetica", font = 2, cex = 1)
    mtext(y_label, 2, 2.5, family = "Helvetica", font = 2, cex = 1)
    mtext("Year", 1, 2.5, family = "Helvetica", font = 2, cex = 1)
    
    colors <- viridis(length(quantiles))
    p_values <- numeric(length(quantiles))
    
    for (j in seq_along(quantiles)) {
      lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
      
      valid_idx <- !is.na(yearly_quants[, j])
      if (sum(valid_idx) > 2) {
        trend_data <- data.frame(
          year = all_years[valid_idx],
          quant = yearly_quants[valid_idx, j]
        )
        lm_fit <- lm(quant ~ year, data = trend_data)
        abline(lm_fit, col = colors[j], lty = 2, lwd = 1.5)
        
        # ONE-SIDED p-value: H1 = intensification rates increasing (positive slope)
        slope       <- summary(lm_fit)$coefficients[2, 1]
        p_two_sided <- summary(lm_fit)$coefficients[2, 4]
        p_values[j] <- if (slope > 0) p_two_sided / 2 else 1 - p_two_sided / 2
      } else {
        p_values[j] <- NA
      }
    }
    
    # Save this model's p-values
    all_model_pvals[[model]] <- p_values
    
    legend_text <- paste0(quantiles * 100, "th %ile (p=",
                          ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
    legend("topright", inset = c(-0.05, -0.225), xpd = "True",
           legend = legend_text, ncol = 2,
           col = colors, lwd = 2, bty = "n", cex = 0.7,
           title = "Quantiles (solid=data, dash=trend)")
  }
  
  dev.off()
  par(mfrow = c(1, 1))
  
  return(list(
    yearly_quants = all_model_quants,
    p_values      = all_model_pvals
  ))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Counts the number of nor'easter occurrences in each grid cell of a specified
#   grid for CFSR and ERA5
#   - No model_names loop; ERA5 is a single dataset
#   - Input is a flat year-keyed list, not nested under model names
#   - Returns a single list(x, y, z) directly
#   - Coordinate parsing handles both comma-separated and space-separated numpy formats
count_storms_in_grid_era5 <- function(noreaster_data, n_lat, n_lon, era5_lat, era5_lon) {
  track_count <- matrix(0, nrow = n_lat, ncol = n_lon)
  total_storms <- 0
  
  for (year_name in names(noreaster_data)) {
    year_data <- noreaster_data[[year_name]]
    if (is.null(year_data) || nrow(year_data) == 0) next
    n_storms <- nrow(year_data)
    
    for (i in 1:n_storms) {
      total_storms <- total_storms + 1
      
      cx_str <- as.character(year_data$cx_p[i])
      cy_str <- as.character(year_data$cy_p[i])
      
      # Strip numpy/bracket formatting
      cx_clean <- gsub("\\[|\\]|np\\.int64\\(|\\)|'f'", "", cx_str)
      cy_clean <- gsub("\\[|\\]|np\\.int64\\(|\\)|'f'", "", cy_str)
      
      # Handle both comma-separated and space-separated formats
      parse_coords <- function(s) {
        if (grepl(",", s)) {
          as.numeric(unlist(strsplit(s, ",")))
        } else {
          as.numeric(unlist(strsplit(trimws(s), "\\s+")))
        }
      }
      
      cx_values <- parse_coords(cx_clean)
      cy_values <- parse_coords(cy_clean)
      
      valid <- !is.na(cx_values) & !is.na(cy_values) &
        cx_values > 0 & cy_values > 0 &
        cx_values <= n_lon & cy_values <= n_lat
      
      if (sum(valid) > 0) {
        cx_valid <- cx_values[valid]
        cy_valid <- cy_values[valid]
        
        cells_visited <- matrix(FALSE, nrow = n_lat, ncol = n_lon)
        
        for (j in seq_along(cx_valid)) {
          lon_idx <- round(cx_valid[j])
          lat_idx <- round(cy_valid[j])
          
          if (lon_idx > 0 && lon_idx <= n_lon &&
              lat_idx > 0 && lat_idx <= n_lat) {
            cells_visited[lat_idx, lon_idx] <- TRUE
          }
        }
        
        track_count <- track_count + cells_visited
      }
    }
  }
  
  cat(sprintf("  Processed %d storms\n", total_storms))
  
  return(list(
    x = era5_lon,
    y = era5_lat,
    z = track_count
  ))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Converts matrix of grid indices into regular lat/lon
mat_to_regular_grid_rect <- function(data_list, lon_out, lat_out) {
  lon_vec <- as.vector(data_list$x)
  lat_vec <- as.vector(data_list$y)
  z_vec   <- as.vector(data_list$z)
  
  i_lon <- sapply(lon_vec, function(l) which.min(abs(lon_out - l)))
  i_lat <- sapply(lat_vec, function(l) which.min(abs(lat_out - l)))
  
  z_out <- matrix(0, nrow = length(lon_out), ncol = length(lat_out))
  
  for (k in seq_along(z_vec)) {
    z_out[i_lon[k], i_lat[k]] <- z_out[i_lon[k], i_lat[k]] + z_vec[k]
  }
  
  list(x = lon_out, y = lat_out, z = z_out)  # explicitly return the output grid vectors
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Count tracks passing through grid cells for nested CRCM5 data structure
count_storms_in_grid <- function(noreaster_data, model_names, n_lat, n_lon, crcm5_lat, crcm5_lon) {
  result <- list()
  
  for (model_name in model_names) {
    cat(paste("Processing", model_name, "...\n"))
    
    track_count <- matrix(0, nrow = n_lat, ncol = n_lon)
    model_data <- noreaster_data[[model_name]]
    total_storms <- 0
    
    for (year_name in names(model_data)) {
      year_data <- model_data[[year_name]]
      n_storms <- nrow(year_data)
      
      for (i in 1:n_storms) {
        total_storms <- total_storms + 1
        
        cx_str <- year_data$cx_p[i]
        cy_str <- year_data$cy_p[i]
        
        cx_clean <- gsub("\\[|\\]|np\\.int64\\(|\\)|'f'", "", cx_str)
        cy_clean <- gsub("\\[|\\]|np\\.int64\\(|\\)|'f'", "", cy_str)
        
        cx_values <- as.numeric(unlist(strsplit(cx_clean, ",")))
        cy_values <- as.numeric(unlist(strsplit(cy_clean, ",")))
        
        valid <- !is.na(cx_values) & !is.na(cy_values) & 
          cx_values > 0 & cy_values > 0 &
          cx_values <= n_lon & cy_values <= n_lat
        
        if (sum(valid) > 0) {
          cx_valid <- cx_values[valid]
          cy_valid <- cy_values[valid]
          
          cells_visited <- matrix(FALSE, nrow = n_lat, ncol = n_lon)
          
          for (j in seq_along(cx_valid)) {
            lon_idx <- round(cx_valid[j])
            lat_idx <- round(cy_valid[j])
            
            if (lon_idx > 0 && lon_idx <= n_lon && 
                lat_idx > 0 && lat_idx <= n_lat) {
              cells_visited[lat_idx, lon_idx] <- TRUE
            }
          }
          
          track_count <- track_count + cells_visited
        }
      }
    }
    
    # Store 2D coordinate matrices and untransposed count matrix
    result[[model_name]] <- list(
      x = crcm5_lon,    # 2D matrix [628 x 655]
      y = crcm5_lat,    # 2D matrix [628 x 655]
      z = track_count   # [628 x 655], no transpose needed
    )
    cat(sprintf("  Processed %d storms\n", total_storms))
  }
  
  return(result)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Converts matrix of grid indices to normal lat/lon grid, specifically for CRCM5
#   data (handles grid cells in the rectangular grid that fall outside the
#   rotated native CRCM5 domain)
mat_to_regular_grid <- function(data_list, lon_out = seq(-110, -7, by = 0.5),
                                lat_out = seq(15, 82, by = 0.5)) {
  lon_vec <- as.vector(data_list$x)
  lat_vec <- as.vector(data_list$y)
  z_vec   <- as.vector(data_list$z)
  
  i_lon <- sapply(lon_vec, function(l) which.min(abs(lon_out - l)))
  i_lat <- sapply(lat_vec, function(l) which.min(abs(lat_out - l)))
  
  z_out <- matrix(0, nrow = length(lon_out), ncol = length(lat_out))  # 0 not NA
  
  for (k in seq_along(z_vec)) {
    z_out[i_lon[k], i_lat[k]] <- z_out[i_lon[k], i_lat[k]] + z_vec[k]
  }
  
  # Set cells outside the CRCM5 domain back to NA
  visited <- matrix(FALSE, nrow = length(lon_out), ncol = length(lat_out))
  visited[cbind(i_lon, i_lat)] <- TRUE
  z_out[!visited] <- NA
  
  list(x = lon_out, y = lat_out, z = z_out)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to extract minimum SLP locations from ERA5 & CFSR data
extract_era5_min_slp_locations <- function(era5_storm_data_1980_2014, lon_grid, lat_grid) {
  cat("=== EXTRACTING ERA5 MINIMUM SLP LOCATIONS ===\n")
  
  all_min_locations <- list()
  
  for (year in names(era5_storm_data_1980_2014)) {
    year_data <- era5_storm_data_1980_2014[[year]]
    
    if (nrow(year_data) == 0) {
      cat("No data for year", year, "\n")
      next
    }
    
    cat("Processing year", year, ":", nrow(year_data), "storms\n")
    
    year_min_locations <- list()
    
    for (i in 1:nrow(year_data)) {
      storm <- year_data[i, ]
      storm_id <- storm$ID
      
      # Parse coordinate arrays - check for both 6-hourly and hourly versions
      cx_p_values <- parse_array_string(storm$cx_p)
      cy_p_values <- parse_array_string(storm$cy_p)
      
      # Try to get hourly coordinate arrays if available
      cx_p_h_values <- NULL
      cy_p_h_values <- NULL
      if ("cx_p_h" %in% names(storm)) {
        cx_p_h_values <- parse_array_string(storm$cx_p_h)
      }
      if ("cy_p_h" %in% names(storm)) {
        cy_p_h_values <- parse_array_string(storm$cy_p_h)
      }
      
      # Parse pressure array
      min_slp_values <- parse_array_string(storm$`c_p (Pa)`)
      
      # Check if we have valid pressure data
      if (is.null(min_slp_values) || length(min_slp_values) == 0) {
        next
      }
      
      # Check if coordinate lengths match pressure data length
      use_hourly_coords <- FALSE
      
      # First check if 6-hourly coordinates match pressure data length
      if (length(cx_p_values) != length(min_slp_values) || length(cy_p_values) != length(min_slp_values)) {
        # Check if hourly coordinates are available and match
        if (!is.null(cx_p_h_values) && !is.null(cy_p_h_values)) {
          if (length(cx_p_h_values) == length(min_slp_values) && length(cy_p_h_values) == length(min_slp_values)) {
            cx_p_values <- cx_p_h_values
            cy_p_values <- cy_p_h_values
            use_hourly_coords <- TRUE
          } else {
            next  # Skip this storm
          }
        } else {
          next
        }
      }
      
      # Ensure we have valid coordinate data
      if (length(cx_p_values) == 0 || length(cy_p_values) == 0) {
        next
      }
      
      # Find the index of minimum SLP
      min_slp_idx <- which.min(min_slp_values)
      
      if (length(min_slp_idx) == 0) {
        next
      }
      
      # Extract minimum SLP location
      min_slp_location <- list(
        storm_id = storm_id,
        year = as.numeric(year),
        min_slp = min_slp_values[min_slp_idx],
        rlat_idx = cy_p_values[min_slp_idx],
        rlon_idx = cx_p_values[min_slp_idx],
        model = "ERA5_raw",
        model_member = "ERA5_raw",
        member = "raw",
        used_hourly_coords = use_hourly_coords,
        n_time_steps = length(min_slp_values)
      )
      
      year_min_locations[[length(year_min_locations) + 1]] <- min_slp_location
    }
    
    all_min_locations[[year]] <- year_min_locations
  }
  
  # Convert to data.table
  min_locations_list <- unlist(all_min_locations, recursive = FALSE)
  
  if (length(min_locations_list) == 0) {
    stop("No valid ERA5 storm data found")
  }
  
  min_locations_dt <- rbindlist(min_locations_list, fill = TRUE)
  
  cat("Extracted", nrow(min_locations_dt), "ERA5 minimum SLP locations\n")
  if (any(min_locations_dt$used_hourly_coords, na.rm = TRUE)) {
    n_hourly <- sum(min_locations_dt$used_hourly_coords, na.rm = TRUE)
    cat("Used hourly coordinates for", n_hourly, "storms\n")
  }
  
  # Convert grid indices to lat/lon coordinates
  cat("Converting ERA5 grid indices to lat/lon coordinates...\n")
  
  min_locations_dt[, c("longitude", "latitude") := {
    coords <- mapply(lat_lon, rlon_idx, rlat_idx, 
                     MoreArgs = list(lon_grid = lon_grid, lat_grid = lat_grid),
                     SIMPLIFY = FALSE)
    coords_matrix <- do.call(rbind, coords)
    list(coords_matrix[,1], coords_matrix[,2])
  }]
  
  # Remove locations with invalid coordinates
  valid_coords <- min_locations_dt[!is.na(longitude) & !is.na(latitude)]
  
  cat("Valid ERA5 coordinate conversions:", nrow(valid_coords), "out of", nrow(min_locations_dt), "\n")
  
  # Summary
  summary_era5 <- valid_coords[, .(
    n_storms = .N,
    mean_min_slp = mean(min_slp, na.rm = TRUE),
    mean_longitude = mean(longitude, na.rm = TRUE),
    mean_latitude = mean(latitude, na.rm = TRUE),
    lon_range = paste(round(range(longitude, na.rm = TRUE), 2), collapse = " to "),
    lat_range = paste(round(range(latitude, na.rm = TRUE), 2), collapse = " to "),
    n_used_hourly = sum(used_hourly_coords, na.rm = TRUE)
  )]
  
  cat("\nERA5 raw minimum SLP locations summary:\n")
  print(summary_era5)
  
  return(valid_coords)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to safely parse array strings from CSV files
parse_array_string <- function(array_str) {
  if (is.na(array_str) || array_str == "" || is.null(array_str)) {
    return(numeric(0))
  }
  
  # Remove brackets and split by whitespace or commas
  clean_str <- gsub("\\[|\\]", "", array_str)
  clean_str <- gsub(",", " ", clean_str)
  
  # Split and convert to numeric
  values <- as.numeric(unlist(strsplit(clean_str, "\\s+")))
  values <- values[!is.na(values)]
  
  return(values)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to extract maximum wind locations from CFSR & ERA5 data
extract_era5_max_wind_locations <- function(era5_storm_data_1980_2014, wind_variable, lon_grid, lat_grid) {
  cat("=== EXTRACTING ERA5 MAXIMUM WIND LOCATIONS ===\n")
  
  all_max_locations <- list()
  
  for (year in names(era5_storm_data_1980_2014)) {
    year_data <- era5_storm_data_1980_2014[[year]]
    
    if (nrow(year_data) == 0) {
      cat("No data for year", year, "\n")
      next
    }
    
    cat("Processing year", year, ":", nrow(year_data), "storms\n")
    
    year_max_locations <- list()
    
    for (i in 1:nrow(year_data)) {
      storm <- year_data[i, ]
      storm <- as.data.frame(storm)
      storm_id <- storm$ID
      
      # Parse coordinate arrays - check for both 6-hourly and hourly versions
      cx_p_values <- parse_array_string(storm$cx_p)
      cy_p_values <- parse_array_string(storm$cy_p)
      
      # Try to get hourly coordinate arrays if available
      cx_p_h_values <- NULL
      cy_p_h_values <- NULL
      if ("cx_p_h" %in% names(storm)) {
        cx_p_h_values <- parse_array_string(storm$cx_p_h)
      }
      else {
        cx_p_values <- parse_array_string(storm$cx_p)
      }
      if ("cy_p_h" %in% names(storm)) {
        cy_p_h_values <- parse_array_string(storm$cy_p_h)
      }
      else {
        cy_p_values <- parse_array_string(storm$cy_p)
      }
      
      # Parse wind speed array - try different possible column names
      max_w_values <- parse_array_string(storm[[wind_variable]])
      
      wind_val <- storm[[wind_variable]]
      if (is.list(wind_val)) {
        max_w_values <- unlist(wind_val)
      } else {
        max_w_values <- parse_array_string(wind_val)
      }
      
      # Check if we have valid wind data
      if (is.null(max_w_values) || length(max_w_values) == 0) {
        cat("Warning: No valid wind speed data for storm", storm_id, "in year", year, "\n")
        next
      }
      
      # Implement the checkpoint you mentioned: check if coordinate lengths match wind length
      use_hourly_coords <- FALSE
      
      # First check if 6-hourly coordinates match wind data length
      if (length(cx_p_values) != length(max_w_values) || length(cy_p_values) != length(max_w_values)) {
        cat("6-hourly coordinates length mismatch for storm", storm_id, 
            ": cx_p=", length(cx_p_values), ", cy_p=", length(cy_p_values), 
            ", max_w=", length(max_w_values), "\n")
        
        # Check if hourly coordinates are available and match
        if (!is.null(cx_p_h_values) && !is.null(cy_p_h_values)) {
          if (length(cx_p_h_values) == length(max_w_values) && length(cy_p_h_values) == length(max_w_values)) {
            cat("Using hourly coordinates (cx_p_h, cy_p_h) instead\n")
            cx_p_values <- cx_p_h_values
            cy_p_values <- cy_p_h_values
            use_hourly_coords <- TRUE
          } else {
            cat("Hourly coordinates also don't match: cx_p_h=", length(cx_p_h_values), 
                ", cy_p_h=", length(cy_p_h_values), "\n")
            next  # Skip this storm
          }
        } else {
          cat("No hourly coordinates available, skipping storm\n")
          next
        }
      }
      
      # Ensure we have valid coordinate data
      if (length(cx_p_values) == 0 || length(cy_p_values) == 0) {
        next
      }
      
      # Find the index of maximum wind speed
      max_wind_idx <- which.max(max_w_values)
      
      if (length(max_wind_idx) == 0) {
        next
      }
      
      # Extract maximum wind location
      max_wind_location <- list(
        storm_id = storm_id,
        year = as.numeric(year),
        max_wind = max_w_values[max_wind_idx],
        rlat_idx = cy_p_values[max_wind_idx],
        rlon_idx = cx_p_values[max_wind_idx],
        model = "ERA5_raw",
        model_member = "ERA5_raw",
        member = "raw",
        used_hourly_coords = use_hourly_coords,
        n_time_steps = length(max_w_values)
      )
      
      year_max_locations[[length(year_max_locations) + 1]] <- max_wind_location
    }
    
    all_max_locations[[year]] <- year_max_locations
  }
  
  # Convert to data.table
  max_locations_list <- unlist(all_max_locations, recursive = FALSE)
  
  if (length(max_locations_list) == 0) {
    stop("No valid ERA5 storm data found")
  }
  
  # Convert to data.table
  max_locations_dt <- rbindlist(max_locations_list, fill = TRUE)
  
  cat("Extracted", nrow(max_locations_dt), "ERA5 maximum wind locations\n")
  if (any(max_locations_dt$used_hourly_coords, na.rm = TRUE)) {
    n_hourly <- sum(max_locations_dt$used_hourly_coords, na.rm = TRUE)
    cat("Used hourly coordinates for", n_hourly, "storms\n")
  }
  
  # Convert grid indices to lat/lon coordinates
  cat("Converting ERA5 grid indices to lat/lon coordinates...\n")
  
  max_locations_dt[, c("longitude", "latitude") := {
    coords <- mapply(lat_lon, rlon_idx, rlat_idx, 
                     MoreArgs = list(lon_grid = lon_grid, lat_grid = lat_grid),
                     SIMPLIFY = FALSE)
    coords_matrix <- do.call(rbind, coords)
    list(coords_matrix[,1], coords_matrix[,2])
  }]
  
  # Remove locations with invalid coordinates
  valid_coords <- max_locations_dt[!is.na(longitude) & !is.na(latitude)]
  
  cat("Valid ERA5 coordinate conversions:", nrow(valid_coords), "out of", nrow(max_locations_dt), "\n")
  
  # Summary
  summary_era5 <- valid_coords[, .(
    n_storms = .N,
    mean_max_wind = mean(max_wind, na.rm = TRUE),
    mean_longitude = mean(longitude, na.rm = TRUE),
    mean_latitude = mean(latitude, na.rm = TRUE),
    lon_range = paste(round(range(longitude, na.rm = TRUE), 2), collapse = " to "),
    lat_range = paste(round(range(latitude, na.rm = TRUE), 2), collapse = " to "),
    n_used_hourly = sum(used_hourly_coords, na.rm = TRUE)
  )]
  
  cat("\nERA5 raw maximum wind locations summary:\n")
  print(summary_era5)
  
  return(valid_coords)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to extract minimum SLP locations from CRCM5 storm data
extract_min_slp_locations <- function(noreaster_data, lon_grid, lat_grid) {
  cat("=== EXTRACTING MINIMUM SLP LOCATIONS ===\n")
  
  # Helper function to parse numpy arrays from CSV
  parse_numpy_array <- function(array_str) {
    if (is.na(array_str) || array_str == "" || array_str == "'f'") {
      return(numeric(0))
    }
    
    # Remove numpy array notation and brackets
    clean_str <- gsub("\\[|\\]", "", array_str)
    clean_str <- gsub("np\\.(int64|float32)\\(([^)]+)\\)", "\\2", clean_str)
    clean_str <- gsub("'f'", "", clean_str)
    clean_str <- trimws(clean_str)
    
    if (clean_str == "" || clean_str == "f") {
      return(numeric(0))
    }
    
    # Split by comma and convert to numeric
    values <- trimws(strsplit(clean_str, ",")[[1]])
    values <- values[values != "" & values != "f" & !is.na(values)]
    
    if (length(values) == 0) {
      return(numeric(0))
    }
    
    return(as.numeric(values))
  }
  
  # Initialize storage for results
  all_results <- list()
  result_counter <- 1
  
  # Process each model
  for (model_name in names(noreaster_data)) {
    cat(paste("Processing model:", model_name, "\n"))
    
    model_data <- noreaster_data[[model_name]]
    model_storms <- 0
    
    # Process each year
    for (year in names(model_data)) {
      year_data <- model_data[[year]]
      
      if (nrow(year_data) == 0) next
      
      # Process each storm in this year
      for (storm_idx in 1:nrow(year_data)) {
        storm <- year_data[storm_idx, ]
        
        tryCatch({
          # Parse coordinate and pressure arrays
          cx_p <- parse_numpy_array(storm$`cx_p`)
          cy_p <- parse_numpy_array(storm$`cy_p`)
          c_p <- parse_numpy_array(storm$`c_p (Pa)`)
          
          # Skip if we don't have valid data
          if (length(cx_p) == 0 || length(cy_p) == 0 || length(c_p) == 0) {
            next
          }
          
          # Ensure all arrays are the same length
          min_length <- min(length(cx_p), length(cy_p), length(c_p))
          if (min_length == 0) next
          
          cx_p <- cx_p[1:min_length]
          cy_p <- cy_p[1:min_length]
          c_p <- c_p[1:min_length]
          
          # Find minimum pressure
          min_idx <- which.min(c_p)
          
          if (length(min_idx) > 0 && !is.na(min_idx)) {
            # Store the result
            all_results[[result_counter]] <- list(
              model = model_name,
              model_member = toupper(model_name),
              member = "1",
              storm_id_clean = paste0(model_name, "_", year, "_", storm_idx),
              min_slp = c_p[min_idx],
              rlat_idx = cy_p[min_idx],
              rlon_idx = cx_p[min_idx],
              data_year = as.numeric(year),
              season = paste0(as.numeric(year)-1, "-", year)
            )
            
            result_counter <- result_counter + 1
            model_storms <- model_storms + 1
          }
          
        }, error = function(e) {
          # Silently skip problematic storms
        })
      }
    }
    
    cat(paste("  ", model_name, ":", model_storms, "storms processed\n"))
  }
  
  
  # Convert results to data.table
  if (length(all_results) == 0) {
    stop("No valid storms found in the data")
  }
  
  min_slp_locations <- rbindlist(all_results, fill = TRUE)
  
  cat("Extracted minimum SLP locations for", nrow(min_slp_locations), "storms\n")
  
  # Convert grid indices to lat/lon coordinates
  cat("Converting grid indices to lat/lon coordinates...\n")
  
  min_slp_locations[, c("longitude", "latitude") := {
    coords <- mapply(lat_lon, rlon_idx, rlat_idx, 
                     MoreArgs = list(lon_grid = lon_grid, lat_grid = lat_grid),
                     SIMPLIFY = FALSE)
    coords_matrix <- do.call(rbind, coords)
    list(coords_matrix[,1], coords_matrix[,2])
  }]
  
  # Remove locations with invalid coordinates
  valid_coords <- min_slp_locations[!is.na(longitude) & !is.na(latitude)]
  
  cat("Valid coordinate conversions:", nrow(valid_coords), "out of", nrow(min_slp_locations), "\n")
  
  # Summary by model
  summary_by_model <- valid_coords[, .(
    n_storms = .N,
    mean_min_slp = mean(min_slp, na.rm = TRUE),
    mean_longitude = mean(longitude, na.rm = TRUE),
    mean_latitude = mean(latitude, na.rm = TRUE),
    lon_range = paste(round(range(longitude, na.rm = TRUE), 2), collapse = " to "),
    lat_range = paste(round(range(latitude, na.rm = TRUE), 2), collapse = " to ")
  ), by = .(model)]
  
  cat("\nMinimum SLP locations summary by model:\n")
  print(summary_by_model)
  
  return(valid_coords)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate kernel density estimate for maximum intensity, handles CRCM5
#     nested data structure
calculate_max_intensity_kde <- function(max_intensity_locations, kde_resolution = 100, 
                                   era5_lon_min, era5_lat_min, 
                                   era5_lon_max, era5_lat_max,
                                   lon = NULL, lat = NULL, average_kdes = NULL,
                                   multi_member_models = c("MPI", "ERA5")) {
  
  cat("=== CALCULATING KDE FOR MAXIMUM WIND LOCATIONS ===\n")
  
  # Define ERA5 domain bounds
  era5_x_min <- era5_lon_min
  era5_x_max <- era5_lon_max
  era5_y_min <- era5_lat_min
  era5_y_max <- era5_lat_max
  
  # Get lon/lat coordinates
  lon_col <- "longitude"
  lat_col <- "latitude"
  
  # Get unique models
  models <- unique(max_intensity_locations$model)
  cat("Models found:", paste(models, collapse = ", "), "\n")
  
  kde_results <- list()
  
  for (model_name in models) {
    cat("\nProcessing model:", model_name, "\n")
    
    # Extract all data for this model
    model_data <- max_intensity_locations[model == model_name]
    
    # Check if this is a multi-member model
    if (toupper(model_name) %in% multi_member_models) {
      cat("  Multi-member model detected - calculating ensemble-mean KDE\n")
      
      # Get unique members for this model
      members <- unique(model_data$model_member)
      cat("  Members found:", paste(members, collapse = ", "), "\n")
      
      member_kdes <- list()
      
      for (member_name in members) {
        cat("    Processing member:", member_name, "\n")
        
        # Extract data for this specific member
        member_data <- model_data[model_member == member_name]
        
        # Get coordinates
        lon_coords <- member_data[[lon_col]]
        lat_coords <- member_data[[lat_col]]
        
        # Remove any remaining NA coordinates
        valid_coords <- !is.na(lon_coords) & !is.na(lat_coords)
        lon_coords <- lon_coords[valid_coords]
        lat_coords <- lat_coords[valid_coords]
        
        cat("      Valid coordinates:", length(lon_coords), "\n")
        
        # Check if we have enough points for KDE
        if (length(lon_coords) > 5) {
          tryCatch({
            # Calculate KDE for this member
            kde_result <- kde2d(lon_coords, lat_coords,
                                n = kde_resolution,
                                lims = c(era5_x_min, era5_x_max, era5_y_min, era5_y_max))
            
            member_kdes[[member_name]] <- kde_result
            cat("      ✓ Success: KDE calculated for member", member_name, 
                "with", length(lon_coords), "points\n")
            
          }, error = function(e) {
            cat("      ✗ Error in KDE calculation for member", member_name, ":", e$message, "\n")
          })
        } else {
          cat("      ✗ Too few points for KDE:", length(lon_coords), "\n")
        }
      }
      
      # Average the member KDEs if we have any
      if (length(member_kdes) > 0) {
        if (!is.null(average_kdes)) {
          kde_results[[tolower(model_name)]] <- average_kdes(member_kdes)
          cat("  ✅ Ensemble-averaged KDE for", model_name, "using", length(member_kdes), "members\n")
        } else {
          cat("  ✗ average_kdes function not provided - cannot create ensemble mean\n")
        }
        
        # Show member breakdown
        cat("    Member breakdown:\n")
        for (member in names(member_kdes)) {
          cat("      ", member, ": included in ensemble\n")
        }
      } else {
        cat("  ✗ No valid member KDEs for", model_name, "\n")
      }
      
    } else {
      # Single-member model - use original approach
      cat("  Single-member model - using standard KDE approach\n")
      
      # Get coordinates
      lon_coords <- model_data[[lon_col]]
      lat_coords <- model_data[[lat_col]]
      
      # Remove any remaining NA coordinates
      valid_coords <- !is.na(lon_coords) & !is.na(lat_coords)
      lon_coords <- lon_coords[valid_coords]
      lat_coords <- lat_coords[valid_coords]
      
      cat("  Valid coordinates:", length(lon_coords), "\n")
      
      # Check if we have enough points for KDE
      if (length(lon_coords) > 5) {
        tryCatch({
          # Calculate KDE
          kde_result <- kde2d(lon_coords, lat_coords,
                              n = kde_resolution,
                              lims = c(era5_x_min, era5_x_max, era5_y_min, era5_y_max))
          
          kde_results[[tolower(model_name)]] <- kde_result
          cat("  ✓ Success: KDE calculated for", model_name, 
              "with", length(lon_coords), "points\n")
          
        }, error = function(e) {
          cat("  ✗ Error in KDE calculation for", model_name, ":", e$message, "\n")
        })
      } else {
        cat("  ✗ Too few points for KDE:", length(lon_coords), "\n")
      }
    }
  }
  
  cat("\n=== KDE CALCULATION SUMMARY ===\n")
  cat("Successfully calculated KDE for", length(kde_results), "models:\n")
  for (model in names(kde_results)) {
    cat("  -", toupper(model), "\n")
  }
  
  return(kde_results)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Extract max. wind locations from nested CRCM5 data structure
extract_max_wind_locations <- function(strong_noreasters, lon_grid, lat_grid) {
  cat("=== EXTRACTING MAXIMUM WIND LOCATIONS (NOR'EASTERS) ===\n")
  
  # Flatten nested list → single data.table
  all_rows <- list()
  
  for (model_name in names(strong_noreasters)) {
    model_data <- strong_noreasters[[model_name]]
    
    for (year_name in names(model_data)) {
      year_df <- model_data[[year_name]]
      
      if (is.null(year_df) || nrow(year_df) == 0) next
      
      dt <- as.data.table(year_df)
      dt[, model     := model_name]
      dt[, data_year := as.integer(year_name)]
      
      all_rows[[length(all_rows) + 1]] <- dt
    }
  }
  
  if (length(all_rows) == 0) stop("No data found in strong_noreasters.")
  
  data <- rbindlist(all_rows, fill = TRUE)
  cat("Combined data dimensions:", nrow(data), "rows,", ncol(data), "cols\n")
  cat("Models found:", paste(unique(data$model), collapse = ", "), "\n")
  
  # Parse Python numpy string columns
  cat("Parsing Python numpy string columns...\n")
  for (col in c("cx_p", "cy_p")) {
    data[, (col) := lapply(get(col), function(x) {
      vals <- extract_coordinates_robust(x)
      if (length(vals) > 0) fix_corrupted_values(vals) else vals
    })]
    cat("  Parsed:", col, "\n")
  }
  
  # Extract lifetime maximum wind and its location
  #     max_wind_in_box is a vector across each storm's lifetime - we find the
  #     timestep of the peak wind, then extract the pressure center coordinates
  #     (cx_p, cy_p) at that same timestep for comparability with ERA5
  cat("Extracting lifetime maximum wind speed and location per storm...\n")
  data[, peak_timestep  := sapply(max_wind_in_box, which.max)]
  data[, max_wind_val   := as.numeric(sapply(max_wind_in_box, max, na.rm = TRUE))]
  data[, max_wind_x_idx := as.numeric(mapply(function(x, i) x[[i]], cx_p, peak_timestep))]
  data[, max_wind_y_idx := as.numeric(mapply(function(y, i) y[[i]], cy_p, peak_timestep))]
  
  # Construct unique storm ID
  data[, storm_id_clean := paste(model, data_year, ID, sep = "_")]
  
  # Verify required columns
  required_cols <- c("storm_id_clean", "max_wind_val",
                     "max_wind_x_idx", "max_wind_y_idx",
                     "season_start_year", "season_end_year")
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         "\nAvailable: ", paste(names(data), collapse = ", "))
  }
  
  # Diagnostics
  cat("\nRow counts by model:\n")
  print(data[, .N, by = model])
  
  missing_by_model <- data[, .(
    missing_storm_id  = sum(is.na(storm_id_clean)),
    missing_max_wind  = sum(is.na(max_wind_val)),
    missing_x_idx     = sum(is.na(max_wind_x_idx)),
    missing_y_idx     = sum(is.na(max_wind_y_idx)),
    total_rows        = .N
  ), by = model]
  cat("\nMissing values by model:\n")
  print(missing_by_model)
  
  # Drop incomplete rows
  data_clean <- data[
    !is.na(storm_id_clean) &
      !is.na(max_wind_val) &
      !is.na(max_wind_x_idx) &
      !is.na(max_wind_y_idx)
  ]
  cat("\nAfter removing missing values:", nrow(data_clean), "rows\n")
  
  cat("\nUnique storms per model:\n")
  print(data_clean[, .(unique_storms = uniqueN(storm_id_clean)), by = model])
  
  # Select columns to keep
  keep_cols <- intersect(
    c("model", "data_year", "storm_id_clean", "ID",
      "max_wind_val", "max_wind_x_idx", "max_wind_y_idx",
      "season_start_year", "season_end_year",
      "rlat_idx", "rlon_idx"),
    names(data_clean)
  )
  max_wind_locations <- data_clean[, keep_cols, with = FALSE]
  
  cat("\nRows after column selection:", nrow(max_wind_locations), "\n")
  cat("By model:\n")
  print(max_wind_locations[, .N, by = model])
  
  # Convert grid indices → lat/lon
  cat("\nConverting grid indices to lat/lon coordinates...\n")
  
  max_wind_locations[, c("longitude", "latitude") := {
    coords <- mapply(
      lat_lon,
      max_wind_x_idx, max_wind_y_idx,
      MoreArgs = list(lon_grid = lon_grid, lat_grid = lat_grid),
      SIMPLIFY = FALSE
    )
    coords_matrix <- do.call(rbind, coords)
    list(coords_matrix[, 1], coords_matrix[, 2])
  }]
  
  # Drop invalid coordinate conversions
  valid_coords <- max_wind_locations[!is.na(longitude) & !is.na(latitude)]
  cat("Valid coordinate conversions:", nrow(valid_coords),
      "out of", nrow(max_wind_locations), "\n")
  
  if (nrow(valid_coords) < nrow(max_wind_locations)) {
    cat("Coordinate conversion losses by model:\n")
    print(max_wind_locations[, .(
      total = .N,
      valid = sum(!is.na(longitude) & !is.na(latitude)),
      lost  = sum(is.na(longitude)  | is.na(latitude))
    ), by = model])
    
    failed <- max_wind_locations[is.na(longitude) | is.na(latitude)]
    cat("Examples of failed conversions:\n")
    print(failed[1:min(5, .N),
                 .(model, storm_id_clean, max_wind_x_idx, max_wind_y_idx)])
  }
  
  # Summary
  summary_by_model <- valid_coords[, .(
    n_storms       = .N,
    mean_max_wind  = mean(max_wind_val,  na.rm = TRUE),
    mean_longitude = mean(longitude,     na.rm = TRUE),
    mean_latitude  = mean(latitude,      na.rm = TRUE),
    lon_range      = paste(round(range(longitude, na.rm = TRUE), 2), collapse = " to "),
    lat_range      = paste(round(range(latitude,  na.rm = TRUE), 2), collapse = " to ")
  ), by = model]
  
  cat("\nSummary by model:\n")
  print(summary_by_model)
  
  return(valid_coords)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Extract maximum intensification rate locations from CFSR and ERA5 data
extract_era5_max_intensification_locations <- function(era5_intens_data, lon_grid, lat_grid) {
  cat("=== EXTRACTING ERA5 MAX INTENSIFICATION RATE LOCATIONS ===\n")
  
  all_results <- list()
  result_counter <- 1
  total_storms <- 0
  
  for (year in names(era5_intens_data)) {
    year_data <- era5_intens_data[[year]]
    
    for (storm_name in names(year_data)) {
      storm <- year_data[[storm_name]]
      
      tryCatch({
        cx         <- storm$cx_at_max_intens
        cy         <- storm$cy_at_max_intens
        max_intens <- storm$max_intensification
        
        if (!is.null(cx) && !is.null(cy) && !is.null(max_intens) &&
            !is.na(cx)  && !is.na(cy)  && !is.na(max_intens)) {
          
          all_results[[result_counter]] <- list(
            storm_id            = as.character(storm$storm_id$ID),
            year                = as.numeric(year),
            max_intensification = max_intens,
            time_window_hours   = storm$time_window_hours,
            rlat_idx            = cy,
            rlon_idx            = cx,
            model               = "era5",
            model_member        = "ERA5",
            member              = "1"
          )
          
          result_counter <- result_counter + 1
          total_storms   <- total_storms + 1
        }
        
      }, error = function(e) {})
    }
  }
  
  if (length(all_results) == 0) {
    stop("No valid ERA5 intensification data found")
  }
  
  intens_locations_dt <- rbindlist(all_results, fill = TRUE)
  cat("Extracted", nrow(intens_locations_dt), "ERA5 max intensification locations\n")
  
  # Convert grid indices to lat/lon coordinates
  cat("Converting ERA5 grid indices to lat/lon coordinates...\n")
  
  intens_locations_dt[, c("longitude", "latitude") := {
    coords <- mapply(lat_lon, rlon_idx, rlat_idx,
                     MoreArgs = list(lon_grid = lon_grid, lat_grid = lat_grid),
                     SIMPLIFY = FALSE)
    coords_matrix <- do.call(rbind, coords)
    list(coords_matrix[, 1], coords_matrix[, 2])
  }]
  
  valid_coords <- intens_locations_dt[!is.na(longitude) & !is.na(latitude)]
  cat("Valid coordinate conversions:", nrow(valid_coords), "out of", nrow(intens_locations_dt), "\n")
  
  summary_era5 <- valid_coords[, .(
    n_storms        = .N,
    mean_max_intens = mean(max_intensification, na.rm = TRUE),
    mean_longitude  = mean(longitude, na.rm = TRUE),
    mean_latitude   = mean(latitude, na.rm = TRUE),
    lon_range       = paste(round(range(longitude, na.rm = TRUE), 2), collapse = " to "),
    lat_range       = paste(round(range(latitude, na.rm = TRUE), 2), collapse = " to ")
  )]
  
  cat("\nERA5 max intensification locations summary:\n")
  print(summary_era5)
  
  return(valid_coords)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Extract locations of maximum intensification from CRCM5 intensification rate data.
# Sibling of extract_min_slp_locations but operates on the list structure produced
# by calculate_max_intensification_rate (nested: model -> year -> storm list).
# Output fields beyond coordinates: max_intensification, peak_value.
extract_max_intensification_locations <- function(intens_data, lon_grid, lat_grid) {
  cat("=== EXTRACTING MAXIMUM INTENSIFICATION LOCATIONS ===\n")
  
  all_results <- list()
  result_counter <- 1
  
  for (model_name in names(intens_data)) {
    cat(paste("Processing model:", model_name, "\n"))
    
    model_data <- intens_data[[model_name]]
    model_storms <- 0
    
    for (year in names(model_data)) {
      year_data <- model_data[[year]]
      
      if (is.null(year_data) || length(year_data) == 0) next
      
      for (storm_name in names(year_data)) {
        storm <- year_data[[storm_name]]
        
        tryCatch({
          max_intens <- storm$max_intensification
          peak_val   <- storm$peak_value
          cx         <- storm$cx_at_max_intens
          cy         <- storm$cy_at_max_intens
          storm_id   <- as.character(storm$storm_id)
          
          # Skip storms that never intensified or have no location
          if (is.null(max_intens) || !is.finite(max_intens)) next
          if (is.null(cx) || is.null(cy) || is.na(cx) || is.na(cy)) next
          
          all_results[[result_counter]] <- list(
            model              = model_name,
            model_member       = toupper(model_name),
            storm_id_clean     = paste0(model_name, "_", year, "_", storm_name),
            storm_id           = storm_id,
            max_intensification = max_intens,
            peak_value         = peak_val,
            rlat_idx           = cy,
            rlon_idx           = cx,
            data_year          = as.numeric(year),
            season             = paste0(as.numeric(year)-1, "-", year)
          )
          result_counter <- result_counter + 1
          model_storms   <- model_storms + 1
          
        }, error = function(e) {
          # Silently skip problematic storms
        })
      }
    }
    
    cat(paste(" ", model_name, ":", model_storms, "storms processed\n"))
  }
  
  if (length(all_results) == 0) {
    stop("No valid storms found in the data")
  }
  
  intens_locations <- rbindlist(all_results, fill = TRUE)
  
  cat("Extracted maximum intensification locations for", nrow(intens_locations), "storms\n")
  
  # Convert grid indices to lat/lon coordinates
  cat("Converting grid indices to lat/lon coordinates...\n")
  
  intens_locations[, c("longitude", "latitude") := {
    coords <- mapply(lat_lon, rlon_idx, rlat_idx,
                     MoreArgs = list(lon_grid = lon_grid, lat_grid = lat_grid),
                     SIMPLIFY = FALSE)
    coords_matrix <- do.call(rbind, coords)
    list(coords_matrix[,1], coords_matrix[,2])
  }]
  
  valid_coords <- intens_locations[!is.na(longitude) & !is.na(latitude)]
  
  cat("Valid coordinate conversions:", nrow(valid_coords), "out of", nrow(intens_locations), "\n")
  
  # Summary by model
  summary_by_model <- valid_coords[, .(
    n_storms            = .N,
    mean_max_intens     = mean(max_intensification, na.rm = TRUE),
    mean_peak_value     = mean(peak_value, na.rm = TRUE),
    mean_longitude      = mean(longitude, na.rm = TRUE),
    mean_latitude       = mean(latitude, na.rm = TRUE),
    lon_range           = paste(round(range(longitude, na.rm = TRUE), 2), collapse = " to "),
    lat_range           = paste(round(range(latitude, na.rm = TRUE), 2), collapse = " to ")
  ), by = .(model)]
  
  cat("\nMaximum intensification locations summary by model:\n")
  print(summary_by_model)
  
  return(valid_coords)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function for extracting specified intensity values from MPI member nor'easters:
extract_mpi_intensity_values <- function(storm_data_1980_2014, metric) {
  results <- list()
  
  members <- c("sim1", "sim2", "sim3", "sim4", "sim5")
  
  for (model in members) {
    if (model %in% names(storm_data_1980_2014)) {
      cat(paste("Processing", model, "...\n"))
      results[[model]] <- list()
      
      for (year_str in names(storm_data_1980_2014[[model]])) {
        year <- as.numeric(year_str)
        #if (year >= 1951 && year <= 2014) {
        year_data <- storm_data_1980_2014[[model]][[year_str]]
        
        if (metric == "min_slp") {
          if (is.data.frame(year_data) && nrow(year_data) > 0) {
            year_storms <- list()
            for (i in 1:nrow(year_data)) {
              min_slp <- year_data[i, "min_p (Pa)"]
              year_storms[[paste0("storm_", i)]] <- list(min_slp = year_data[i, "min_p (Pa)"])
            }
            results[[model]][[year_str]] <- year_storms
          }
        }
        if (metric == "max_wind") {
          if (is.data.frame(year_data) && nrow(year_data) > 0) {
            year_storms <- list()
            for (i in 1:nrow(year_data)) {
              max_wind <- max(year_data[[i, "max_wind_in_box"]])
              year_storms[[paste0("storm_", i)]] <- list(max_wind = max_wind)
            }
            results[[model]][[year_str]] <- year_storms
          }
        }
      }
    }
  }
  
  return(results)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to calculate pooled MPI intensity quantiles for each member
calculate_mpi_intensity_quantiles <- function(member_data, metric, probs = my.probs) {
  all_values <- numeric()
  
  # Loop through years
  for (year in names(member_data)) {
    storms <- member_data[[year]]
    if (is.null(storms) || length(storms) == 0) next
    
    # Loop through storms
    for (storm in storms) {
      if (metric == "min_slp") {
        val <- storm[["min_slp"]]
        if (is.null(val)) val <- storm[["min_p (Pa)"]]
      } else if (metric == "max_wind") {
        val <- storm[["max_wind"]]
      }
      
      if (!is.null(val) && !is.na(val)) {
        all_values <- c(all_values, as.numeric(val))
      }
    }
  }
  
  # Calculate quantiles from combined values
  if (length(all_values) > 0) {
    return(quantile(all_values, probs = probs, na.rm = TRUE))
  } else {
    return(NA)
  }
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to find storm by partial ID match (used in TC filter testing)
find_storm_by_id <- function(data, partial_id) {
  for (year in names(data)) {
    year_data <- data[[year]]
    if (is.data.frame(year_data)) {
      matches <- grep(partial_id, year_data$ID, value = FALSE)
      if (length(matches) > 0) {
        return(list(
          year = year,
          storm_idx = matches,
          storm_data = year_data[matches, ]
        ))
      }
    }
  }
  return(NULL)
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to plot a single CFSR storm track (used in TC filter testing)
plot_cfsr_storm_track <- function(storm_data, storm_name, lon_grid, lat_grid) {
  
  # Extract track coordinates - use [[1]] to unwrap list column
  cx_p_coords <- extract_coordinates_robust(storm_data$cx_p[[1]])
  cy_p_coords <- extract_coordinates_robust(storm_data$cy_p[[1]])
  
  # Fix corrupted values
  cx_p_coords <- fix_corrupted_values(cx_p_coords)
  cy_p_coords <- fix_corrupted_values(cy_p_coords)
  
  # Convert to lat/lon
  track_lons <- numeric(length(cx_p_coords))
  track_lats <- numeric(length(cy_p_coords))
  
  for (i in seq_along(cx_p_coords)) {
    coords <- lat_lon(cx_p_coords[i], cy_p_coords[i], lon_grid, lat_grid)
    track_lons[i] <- coords[1]
    track_lats[i] <- coords[2]
  }
  
  # Extract wind speeds from list column
  winds_all <- as.numeric(unlist(storm_data$max_wind_in_box))
  winds_all <- winds_all[!is.na(winds_all)]
  
  # Guard against empty wind data
  if (length(winds_all) == 0) {
    winds_all <- rep(NA, length(track_lons))
  }
  
  # Trim to match track length if needed (should already be 1:1)
  n_points <- min(length(track_lons), length(winds_all))
  track_lons <- track_lons[1:n_points]
  track_lats <- track_lats[1:n_points]
  winds_plot <- winds_all[1:n_points]
  
  # Set fixed map bounds for Atlantic basin
  xlim <- c(-100, -10)
  ylim <- c(10, 70)
  
  # Create base map with corrected aspect ratio
  map("worldHires", xlim = xlim, ylim = ylim, 
      col = "gray90", fill = TRUE, border = "gray60",
      asp = 1.3)
  map.axes()
  
  # Add grid
  abline(h = seq(0, 90, 10), col = "gray80", lty = 2, lwd = 0.5)
  abline(v = seq(-180, 0, 10), col = "gray80", lty = 2, lwd = 0.5)
  
  # Add 30°N reference line (TC filter threshold)
  abline(h = 30, col = "darkgreen", lty = 2, lwd = 2)
  text(-95, 31, "30°N", col = "darkgreen", cex = 1.25, font = 2)
  
  # Plot track lines connecting points
  lines(track_lons, track_lats, col = "gray40", lwd = 1.75)
  
  # Plot each time step as a colored point based on intensity
  for (i in seq_along(track_lons)) {
    category <- get_ss_category(winds_plot[i])
    color <- get_ss_color(category)
    
    points(track_lons[i], track_lats[i], pch = 21, 
           col = "black", bg = color, cex = 1.5, lwd = 1)
  }
  
  # Mark genesis (first point) with special marker
  points(track_lons[1], track_lats[1], pch = 23, 
         col = "black", bg = "green", cex = 2.5, lwd = 2)
  
  # Add title with storm info
  title(main = storm_name, font.main = 2, cex.main = 2)
  
  # Add info text
  info_text <- paste0(
    "Genesis: ", round(track_lats[1], 1), "°N\n",
    "Min P: ", round(storm_data$`min_p (Pa)`/100, 1), " hPa\n",
    "Max Wind: ", round(max(unlist(storm_data$max_wind_in_box), na.rm = TRUE), 1), " m/s"
  )
  
  legend("topleft", legend = info_text, bty = "n", cex = 1.5)
  
  # Add Saffir-Simpson legend
  legend("bottomright", 
         legend = c("TD", "TS", "Cat 1", "Cat 2", "Cat 3", "Cat 4", "Cat 5", "Genesis"), 
         pch = c(rep(21, 7), 23),
         pt.bg = c(get_ss_color("TD"), get_ss_color("TS"), 
                   get_ss_color("Cat1"), get_ss_color("Cat2"),
                   get_ss_color("Cat3"), get_ss_color("Cat4"),
                   get_ss_color("Cat5"), "green"),
         col = "black",
         pt.cex = 1.75,
         pt.lwd = 1.25,
         ncol = 2,
         bty = "n", 
         cex = 1.5,
         title = "Intensity")
  
  # Add box
  box()
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Function to plot a single ERA5 storm track (used in TC filter testing)
plot_era5_storm_track <- function(storm_data, storm_name, lon_grid, lat_grid) {
  
  # Extract track coordinates
  cx_p_coords <- extract_coordinates_robust(storm_data$cx_p)
  cy_p_coords <- extract_coordinates_robust(storm_data$cy_p)
  
  # Fix corrupted values
  cx_p_coords <- fix_corrupted_values(cx_p_coords)
  cy_p_coords <- fix_corrupted_values(cy_p_coords)
  
  # Convert to lat/lon
  track_lons <- numeric(length(cx_p_coords))
  track_lats <- numeric(length(cy_p_coords))
  
  for (i in seq_along(cx_p_coords)) {
    coords <- lat_lon(cx_p_coords[i], cy_p_coords[i], lon_grid, lat_grid)
    track_lons[i] <- coords[1]
    track_lats[i] <- coords[2]
  }
  
  # Extract wind speeds - subsample to match 6-hourly track positions
  wind_string <- gsub("\\[|\\]", "", storm_data$`max_w_h_750 (m/s)`)
  winds_all <- as.numeric(unlist(strsplit(wind_string, "\\s+")))
  
  # Guard against empty wind data
  if (length(winds_all) == 0) {
    winds_all <- rep(NA, length(track_lons))
  }
  
  # Subsample winds to match track (every 6th value for 6-hourly data)
  wind_indices <- seq(1, length(winds_all), by = 6)
  winds_6hourly <- winds_all[wind_indices]
  
  # Trim to match track length if needed
  n_points <- min(length(track_lons), length(winds_6hourly))
  track_lons <- track_lons[1:n_points]
  track_lats <- track_lats[1:n_points]
  winds_6hourly <- winds_6hourly[1:n_points]
  
  # Set fixed map bounds for Atlantic basin
  xlim <- c(-100, -10)
  ylim <- c(10, 70)
  
  # Create base map with corrected aspect ratio
  map("worldHires", xlim = xlim, ylim = ylim, 
      col = "gray90", fill = TRUE, border = "gray60",
      asp = 1.3)
  map.axes()
  
  # Add grid
  abline(h = seq(0, 90, 10), col = "gray80", lty = 2, lwd = 0.5)
  abline(v = seq(-180, 0, 10), col = "gray80", lty = 2, lwd = 0.5)
  
  # Add 30°N reference line (TC filter threshold)
  abline(h = 30, col = "darkgreen", lty = 2, lwd = 2)
  text(-95, 31, "30°N", col = "darkgreen", cex = 1.25, font = 2)
  
  # Plot track lines connecting points
  lines(track_lons, track_lats, col = "gray40", lwd = 1.75)
  
  # Plot each time step as a colored point based on intensity
  for (i in seq_along(track_lons)) {
    category <- get_ss_category(winds_6hourly[i])
    color <- get_ss_color(category)
    
    points(track_lons[i], track_lats[i], pch = 21, 
           col = "black", bg = color, cex = 1.5, lwd = 1)
  }
  
  # Mark genesis (first point) with special marker
  points(track_lons[1], track_lats[1], pch = 23, 
         col = "black", bg = "green", cex = 2.5, lwd = 2)
  
  # Add title with storm info
  title(main = storm_name, font.main = 2, cex.main = 2)
  
  # Add info text
  info_text <- paste0(
    "Genesis: ", round(track_lats[1], 1), "°N\n",
    "Min P: ", round(storm_data$`min_p_h (Pa)`/100, 1), " hPa\n",
    "Max Wind: ", round(storm_data$`max_w_all_h_750 (m/s)`, 1), " m/s"
  )
  
  legend("topleft", legend = info_text, bty = "n", cex = 1.5)
  
  # Add Saffir-Simpson legend
  legend("bottomright", 
         legend = c("TD", "TS", "Cat 1", "Cat 2", "Cat 3", "Cat 4", "Cat 5", "Genesis"), 
         pch = c(rep(21, 7), 23),
         pt.bg = c(get_ss_color("TD"), get_ss_color("TS"), 
                   get_ss_color("Cat1"), get_ss_color("Cat2"),
                   get_ss_color("Cat3"), get_ss_color("Cat4"),
                   get_ss_color("Cat5"), "green"),
         col = "black",
         pt.cex = 1.75,
         pt.lwd = 1.25,
         ncol = 2,
         bty = "n", 
         cex = 1.5,
         title = "Intensity")
  
  # Add box
  box()
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Define Saffir-Simpson scale (based on max sustained wind speed in m/s)
#     TD: < 17 m/s, TS: 17-32 m/s, Cat1: 33-42 m/s, Cat2: 43-49 m/s, 
#     Cat3: 50-58 m/s, Cat4: 59-70 m/s, Cat5: > 70 m/s
#     (used in TC filter testing)
get_ss_category <- function(wind_ms) {
  if (is.na(wind_ms)) return("TD")
  if (wind_ms < 17) return("TD")
  if (wind_ms < 33) return("TS")
  if (wind_ms < 43) return("Cat1")
  if (wind_ms < 50) return("Cat2")
  if (wind_ms < 58) return("Cat3")
  if (wind_ms < 70) return("Cat4")
  return("Cat5")
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Specifies colors for each Saffir-Simpson intensity level (used in TC filter testing)
get_ss_color <- function(category) {
  colors <- c(
    "TD" = "#5EBAFF",      # Light blue (Tropical Depression)
    "TS" = "#00FAF4",      # Cyan (Tropical Storm)
    "Cat1" = "#FFFFCC",    # Pale yellow (Category 1)
    "Cat2" = "#FFE775",    # Yellow (Category 2)
    "Cat3" = "#FFC140",    # Orange (Category 3)
    "Cat4" = "#FF8F20",    # Dark orange (Category 4)
    "Cat5" = "#FF6060"     # Red (Category 5)
  )
  return(colors[category])
}






