# ==============================================================================

# Revision 1: Additional analyses of historical nor'easter intensity characteristics 
#   in CFSR, ERA5, CRCM5, and CRCM5-CMIP6

# ==============================================================================

# Associated publication: Weaver, M.M. et al. (2026). Historical nor'easter
#                         intensity trends in reanalysis and regional climate
#                         simulations: A cross-dataset perspective.

# This file contains additional analyses performed after the first round of 
#   revisions for the associated publication. This includes Mann-Kendall testing
#   of nor'easter intensity and intensification rate trends, extension of ERA5
#   and CFSR intensity and intensification rate trends to present, Figures
#   S7-S19 (Mann-Kendall trends) and Figures S29-30 (CFSR/ERA5 trends to present).
#     - Fig. S8-12 & S14-19: Mann-Kendall intensity and intensification rate trends
#       (different datasets, seasons, time periods, pilot models, etc.)
#     - Fig. 7 & S13: Mann-Kendall intensity and intensification rate p-values
#       w/ and w/o TCs included
#     - Fig. S29-30: intensity and intensification rate trends for ERA5 and CFSR
#       extended to present (1950-2025 and 1979-2025)

# ==============================================================================

# Load necessary .R script(s): contains all previous analyses except TC filtering
source("/Path/to/analysis/file/hist_noreaster_intensity_analysis.R")

# Install necessary packages (new)
install.packages(c("Kendall", "trend"))

# ==============================================================================

# New functions: Mann-Kendall intensity & intensification rate trends

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate (& plot) specific quantiles over time for ERA5 & CFSR with ONE-TAILED
# Mann-Kendall tests (H1 = increasing intensity, either decreasing min. SLP or
# increasing max. wind)
plot_era5_intensity_quantile_trends_mk <- function(intensity_data, metric,
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
        
        # guard against length-zero before is.finite()
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
  
  mtext("Quantile Trends: ERA5 (Mann-Kendall)", 3, 1, adj = 0, 
        family = "Helvetica", font = 2, cex = 1)
  mtext(y_label, 2, 2.5, family = "Helvetica", font = 2, cex = 1)
  mtext("Year", 1, 2.5, family = "Helvetica", font = 2, cex = 1)
  
  colors <- viridis(length(quantiles))
  p_values <- numeric(length(quantiles))
  tau_values <- numeric(length(quantiles))
  
  for (j in seq_along(quantiles)) {
    lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
    
    valid_idx <- !is.na(yearly_quants[, j])
    if (sum(valid_idx) > 2) {
      trend_data <- data.frame(
        year = all_years[valid_idx],
        quant = yearly_quants[valid_idx, j]
      )
      
      # Mann-Kendall test (rank-based, nonparametric)
      mk_fit <- Kendall::MannKendall(trend_data$quant)
      tau <- as.numeric(mk_fit$tau)
      p_two_tailed <- as.numeric(mk_fit$sl[1])
      
      # Sen's slope for the visual trend line
      sen <- trend::sens.slope(trend_data$quant)
      slope <- as.numeric(sen$estimates)
      intercept <- median(trend_data$quant) - slope * median(trend_data$year)
      abline(a = intercept, b = slope, col = colors[j], lty = 2, lwd = 1.5)
      
      tau_values[j] <- tau
      
      # ONE-SIDED p-value based on the sign of tau
      if (metric == "max_wind") {
        p_values[j] <- if (tau > 0) p_two_tailed / 2 else 1 - (p_two_tailed / 2)
      } else if (metric == "min_slp") {
        p_values[j] <- if (tau < 0) p_two_tailed / 2 else 1 - (p_two_tailed / 2)
      }
    } else {
      p_values[j] <- NA
      tau_values[j] <- NA
    }
  }
  
  legend_text <- paste0(quantiles * 100, "th %ile (p=",
                        ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
  legend("topright", inset = c(-0.05, -0.225), xpd = "True",
         legend = legend_text, ncol = 2,
         col = colors, lwd = 2, bty = "n", cex = 0.7,
         title = "Quantiles (solid=data, dash=Sen's slope)")
  
  dev.off()
  par(mfrow = c(1, 1))
  
  return(list(
    quantiles = yearly_quants,
    p_values = p_values,
    tau = tau_values
  ))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate & plot CRCM5-CMIP6 ensemble quantiles (pooling all models together)
# One-tailed Mann-Kendall test (H1 = increasing intensity)
plot_ensemble_quantile_trends_mk <- function(intensity_data, metric,
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
    
    for (model in names(intensity_data)) {
      if (year_str %in% names(intensity_data[[model]])) {
        storms <- intensity_data[[model]][[year_str]]
        
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
    
    values <- values[is.finite(values)]
    if (length(values) > 0) {
      yearly_quants[i, ] <- quantile(values, probs = quantiles)
    }
  }
  
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
  tau_values <- numeric(length(quantiles))
  
  for (j in seq_along(quantiles)) {
    lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
    
    valid_idx <- !is.na(yearly_quants[, j])
    if (sum(valid_idx) > 2) {
      trend_data <- data.frame(
        year = all_years[valid_idx],
        quant = yearly_quants[valid_idx, j]
      )
      
      mk_fit <- Kendall::MannKendall(trend_data$quant)
      tau <- as.numeric(mk_fit$tau)
      p_two_sided <- as.numeric(mk_fit$sl[1])
      
      sen <- trend::sens.slope(trend_data$quant)
      slope <- as.numeric(sen$estimates)
      intercept <- median(trend_data$quant) - slope * median(trend_data$year)
      abline(a = intercept, b = slope, col = colors[j], lty = 2, lwd = 1.5)
      
      tau_values[j] <- tau
      
      # H1: min_slp decreases (tau < 0) or max_wind increases (tau > 0)
      if (metric == "min_slp") {
        p_values[j] <- if (tau < 0) p_two_sided / 2 else 1 - p_two_sided / 2
      } else {
        p_values[j] <- if (tau > 0) p_two_sided / 2 else 1 - p_two_sided / 2
      }
    } else {
      p_values[j] <- NA
      tau_values[j] <- NA
    }
  }
  
  legend_text <- paste0(quantiles * 100, "th %ile (p=",
                        ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
  
  legend("topright", inset = c(0, -0.175), xpd = "True",
         legend = legend_text, ncol = 2,
         col = colors, lwd = 2, bty = "n", cex = 0.7,
         title = "Quantiles (solid=data, dash=Sen's slope)")
  
  if (!is.null(file_path)) {
    dev.off()
  }
  
  return(list(quantiles = yearly_quants, p_values = p_values, tau = tau_values))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate & plot specific CRCM5 quantiles over time BY MODEL (not pooled ensemble)
# One-tailed Mann-Kendall test (H1 = increasing intensity)
plot_intensity_quantile_trends_mk <- function(intensity_data, metric, 
                                              quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9), 
                                              file_path = "file_name") {
  all_years <- sort(unique(unlist(lapply(intensity_data, function(m) as.numeric(names(m))))))
  
  all_model_quants <- list()
  all_model_pvals <- list()
  all_model_tau <- list()
  
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
    
    all_model_quants[[model]] <- yearly_quants
    
    plot(all_years, yearly_quants[, 1], type = "n",
         ylim = range(yearly_quants, na.rm = TRUE),
         xlab = "", ylab = "")
    
    mtext(paste("Quantile Trends (Mann-Kendall):", model), 3, 1, adj = 0, 
          family = "Helvetica", font = 2, cex = 1)
    mtext("Minimum Pressure (Pa)", 2, 2.5, family = "Helvetica", font = 2, cex = 1)
    mtext("Year", 1, 2.5, family = "Helvetica", font = 2, cex = 1)
    
    colors <- viridis(length(quantiles))
    p_values <- numeric(length(quantiles))
    tau_values <- numeric(length(quantiles))
    
    for (j in seq_along(quantiles)) {
      lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
      
      valid_idx <- !is.na(yearly_quants[, j])
      if (sum(valid_idx) > 2) {
        trend_data <- data.frame(
          year = all_years[valid_idx],
          quant = yearly_quants[valid_idx, j]
        )
        
        mk_fit <- Kendall::MannKendall(trend_data$quant)
        tau <- as.numeric(mk_fit$tau)
        p_two_sided <- as.numeric(mk_fit$sl[1])
        
        sen <- trend::sens.slope(trend_data$quant)
        slope <- as.numeric(sen$estimates)
        intercept <- median(trend_data$quant) - slope * median(trend_data$year)
        abline(a = intercept, b = slope, col = colors[j], lty = 2, lwd = 1.5)
        
        tau_values[j] <- tau
        
        if (metric == "max_wind") {
          p_values[j] <- if (tau > 0) p_two_sided / 2 else 1 - p_two_sided / 2
        } else {
          p_values[j] <- if (tau < 0) p_two_sided / 2 else 1 - p_two_sided / 2
        }
      } else {
        p_values[j] <- NA
        tau_values[j] <- NA
      }
    }
    
    all_model_pvals[[model]] <- p_values
    all_model_tau[[model]] <- tau_values
    
    legend_text <- paste0(quantiles * 100, "th %ile (p=",
                          ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
    legend("topright", inset = c(-0.05, -0.225), xpd = "True",
           legend = legend_text, ncol = 2,
           col = colors, lwd = 2, bty = "n", cex = 0.7,
           title = "Quantiles (solid=data, dash=Sen's slope)")
  }
  
  dev.off()
  par(mfrow = c(1, 1))
  
  return(list(
    yearly_quants = all_model_quants,
    p_values      = all_model_pvals,
    tau           = all_model_tau
  ))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate quantiles for CFSR and ERA5 intensification rates
#   (one-tailed Mann-Kendall test, H1 = increasing rates)
plot_era5_intensification_quantile_trends_mk <- function(intensification_data, 
                                                         metric = "min_slp", 
                                                         quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.97, 0.99), 
                                                         file_path = NULL) {
  
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
       main = "ERA5 Intensification Quantile Trends (Mann-Kendall)",
       xlab = "Year", ylab = ylabel)
  
  colors <- viridis(length(quantiles))
  p_values <- numeric(length(quantiles))
  tau_values <- numeric(length(quantiles))
  
  for (j in seq_along(quantiles)) {
    lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
    
    valid_idx <- !is.na(yearly_quants[, j])
    if (sum(valid_idx) > 2) {
      trend_data <- data.frame(
        year = all_years[valid_idx],
        quant = yearly_quants[valid_idx, j]
      )
      
      # Mann-Kendall test (rank-based, nonparametric)
      mk_fit <- Kendall::MannKendall(trend_data$quant)
      tau <- as.numeric(mk_fit$tau)
      p_two_sided <- as.numeric(mk_fit$sl[1])
      
      # Sen's slope for the visual trend line
      sen <- trend::sens.slope(trend_data$quant)
      slope <- as.numeric(sen$estimates)
      intercept <- median(trend_data$quant) - slope * median(trend_data$year)
      abline(a = intercept, b = slope, col = colors[j], lty = 2, lwd = 1.5)
      
      tau_values[j] <- tau
      
      # ONE-SIDED p-value: H1 = intensification rates increasing (tau > 0)
      p_values[j] <- if (tau > 0) p_two_sided / 2 else 1 - p_two_sided / 2
    } else {
      p_values[j] <- NA
      tau_values[j] <- NA
    }
  }
  
  legend_text <- paste0(quantiles * 100, "th %ile (p=", 
                        ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
  legend("topright", inset = c(0, -0.175), xpd = TRUE,
         legend = legend_text, ncol = 2,
         col = colors, lwd = 2, bty = "n", cex = 0.7,
         title = "Quantiles (solid=data, dash=Sen's slope)")
  
  if (!is.null(file_path)) dev.off()
  
  return(list(quantiles = yearly_quants, p_values = p_values, tau = tau_values))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Calculate CRCM5-CMIP6 ensemble quantiles for intensification rates (pooling all models together)
#     (one-tailed Mann-Kendall test, H1 = increasing rates)
plot_ensemble_intensification_quantile_trends_mk <- function(intensification_data, 
                                                             quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.97, 0.99), 
                                                             file_path = NULL) {
  
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
  tau_values <- numeric(length(quantiles))
  
  for (j in seq_along(quantiles)) {
    lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
    
    valid_idx <- !is.na(yearly_quants[, j])
    if (sum(valid_idx) > 2) {
      trend_data <- data.frame(
        year = all_years[valid_idx],
        quant = yearly_quants[valid_idx, j]
      )
      
      mk_fit <- Kendall::MannKendall(trend_data$quant)
      tau <- as.numeric(mk_fit$tau)
      p_two_sided <- as.numeric(mk_fit$sl[1])
      
      sen <- trend::sens.slope(trend_data$quant)
      slope <- as.numeric(sen$estimates)
      intercept <- median(trend_data$quant) - slope * median(trend_data$year)
      abline(a = intercept, b = slope, col = colors[j], lty = 2, lwd = 1.5)
      
      tau_values[j] <- tau
      
      # ONE-SIDED p-value: H1 = intensification rates increasing (tau > 0)
      p_values[j] <- if (tau > 0) p_two_sided / 2 else 1 - p_two_sided / 2
    } else {
      p_values[j] <- NA
      tau_values[j] <- NA
    }
  }
  
  legend_text <- paste0(quantiles * 100, "th %ile (p=", 
                        ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
  legend("topright", inset = c(0, -0.175), xpd = TRUE,
         legend = legend_text, ncol = 2,
         col = colors, lwd = 2, bty = "n", cex = 0.7,
         title = "Quantiles (solid=data, dash=Sen's slope)")
  
  if (!is.null(file_path)) dev.off()
  
  return(list(quantiles = yearly_quants, p_values = p_values, tau = tau_values))
}

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- -- ---

# Plot intensification rate quantiles over time for each CRCM5-CMIP6 pilot model
#   (one-tailed Mann-Kendall test, H1 = increasing rates)
plot_intensification_quantile_trends_mk <- function(intensification_data, 
                                                    metric = "min_slp", 
                                                    quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9), 
                                                    file_path = "file_name") {
  all_years <- sort(unique(unlist(lapply(intensification_data, function(m) as.numeric(names(m))))))
  
  # Initialize lists to collect results for all models
  all_model_quants <- list()
  all_model_pvals <- list()
  all_model_tau <- list()
  
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
    mtext(paste("Intensification Quantile Trends:", model), 3, 1, adj = 0, 
          family = "Helvetica", font = 2, cex = 1)
    mtext(y_label, 2, 2.5, family = "Helvetica", font = 2, cex = 1)
    mtext("Year", 1, 2.5, family = "Helvetica", font = 2, cex = 1)
    
    colors <- viridis(length(quantiles))
    p_values <- numeric(length(quantiles))
    tau_values <- numeric(length(quantiles))
    
    for (j in seq_along(quantiles)) {
      lines(all_years, yearly_quants[, j], col = colors[j], lwd = 2)
      
      valid_idx <- !is.na(yearly_quants[, j])
      if (sum(valid_idx) > 2) {
        trend_data <- data.frame(
          year = all_years[valid_idx],
          quant = yearly_quants[valid_idx, j]
        )
        
        mk_fit <- Kendall::MannKendall(trend_data$quant)
        tau <- as.numeric(mk_fit$tau)
        p_two_sided <- as.numeric(mk_fit$sl[1])
        
        sen <- trend::sens.slope(trend_data$quant)
        slope <- as.numeric(sen$estimates)
        intercept <- median(trend_data$quant) - slope * median(trend_data$year)
        abline(a = intercept, b = slope, col = colors[j], lty = 2, lwd = 1.5)
        
        tau_values[j] <- tau
        
        # ONE-SIDED p-value: H1 = intensification rates increasing (tau > 0)
        p_values[j] <- if (tau > 0) p_two_sided / 2 else 1 - p_two_sided / 2
      } else {
        p_values[j] <- NA
        tau_values[j] <- NA
      }
    }
    
    # Save this model's p-values and tau
    all_model_pvals[[model]] <- p_values
    all_model_tau[[model]] <- tau_values
    
    legend_text <- paste0(quantiles * 100, "th %ile (p=",
                          ifelse(is.na(p_values), "NA", sprintf("%.3f", p_values)), ")")
    legend("topright", inset = c(-0.05, -0.225), xpd = "True",
           legend = legend_text, ncol = 2,
           col = colors, lwd = 2, bty = "n", cex = 0.7,
           title = "Quantiles (solid=data, dash=Sen's slope)")
  }
  
  dev.off()
  par(mfrow = c(1, 1))
  
  return(list(
    yearly_quants = all_model_quants,
    p_values      = all_model_pvals,
    tau           = all_model_tau
  ))
}

# ==============================================================================

# Calculate intensity quantiles over time & trends for each dataset, using
#   Mann-Kendall test

my.probs = seq(from = 0, to = 1, by = 0.01)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CFSR, 1979-2014

# Full dataset
cfsr_slp_quant_trends_mk_1979_2014 <- plot_era5_intensity_quantile_trends_mk(cfsr_min_slp_1979_2014,
                                                                       "min_slp", 
                                                                       quantiles = my.probs)
cfsr_wind_quant_trends_mk_1979_2014 <- plot_era5_intensity_quantile_trends_mk(cfsr_max_wind_1979_2014, 
                                                                        "max_wind", 
                                                                        quantiles = my.probs)
# SON
cfsr_slp_quant_trends_mk_1979_2014_SON <- plot_era5_intensity_quantile_trends_mk(cfsr_min_slp_1979_2014_SON, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
cfsr_wind_quant_trends_mk_1979_2014_SON <- plot_era5_intensity_quantile_trends_mk(cfsr_max_wind_1979_2014_SON, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# DJF
cfsr_slp_quant_trends_mk_1979_2014_DJF <- plot_era5_intensity_quantile_trends_mk(cfsr_min_slp_1979_2014_DJF, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
cfsr_wind_quant_trends_mk_1979_2014_DJF <- plot_era5_intensity_quantile_trends_mk(cfsr_max_wind_1979_2014_DJF, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# MA
cfsr_slp_quant_trends_mk_1979_2014_MA <- plot_era5_intensity_quantile_trends_mk(cfsr_min_slp_1979_2014_MA, 
                                                                          "min_slp", 
                                                                          quantiles = my.probs)
cfsr_wind_quant_trends_mk_1979_2014_MA <- plot_era5_intensity_quantile_trends_mk(cfsr_max_wind_1979_2014_MA, 
                                                                           "max_wind", 
                                                                           quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# ERA5

# 1950-2014
era5_slp_quant_trends_mk_1950_2014 <- plot_era5_intensity_quantile_trends_mk(era5_min_slp_1950_2014,
                                                                       "min_slp", 
                                                                       quantiles = my.probs)
era5_wind_quant_trends_mk_1950_2014 <- plot_era5_intensity_quantile_trends_mk(era5_max_wind_1950_2014, 
                                                                        "max_wind", 
                                                                        quantiles = my.probs)
# SON
era5_slp_quant_trends_mk_1950_2014_SON <- plot_era5_intensity_quantile_trends_mk(era5_min_slp_1950_2014_SON, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
era5_wind_quant_trends_mk_1950_2014_SON <- plot_era5_intensity_quantile_trends_mk(era5_max_wind_1950_2014_SON, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# DJF
era5_slp_quant_trends_mk_1950_2014_DJF <- plot_era5_intensity_quantile_trends_mk(era5_min_slp_1950_2014_DJF, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
era5_wind_quant_trends_mk_1950_2014_DJF <- plot_era5_intensity_quantile_trends_mk(era5_max_wind_1950_2014_DJF, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# MA
era5_slp_quant_trends_mk_1950_2014_MA <- plot_era5_intensity_quantile_trends_mk(era5_min_slp_1950_2014_MA, 
                                                                          "min_slp", 
                                                                          quantiles = my.probs)
era5_wind_quant_trends_mk_1950_2014_MA <- plot_era5_intensity_quantile_trends_mk(era5_max_wind_1950_2014_MA, 
                                                                           "max_wind", 
                                                                           quantiles = my.probs)

# 1979-2014
era5_slp_quant_trends_mk_1979_2014 <- plot_era5_intensity_quantile_trends_mk(era5_min_slp_1979_2014,
                                                                       "min_slp", 
                                                                       quantiles = my.probs)
era5_wind_quant_trends_mk_1979_2014 <- plot_era5_intensity_quantile_trends_mk(era5_max_wind_1979_2014, 
                                                                        "max_wind", 
                                                                        quantiles = my.probs)
# SON
era5_slp_quant_trends_mk_1979_2014_SON <- plot_era5_intensity_quantile_trends_mk(era5_min_slp_1979_2014_SON, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
era5_wind_quant_trends_mk_1979_2014_SON <- plot_era5_intensity_quantile_trends_mk(era5_max_wind_1979_2014_SON, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# DJF
era5_slp_quant_trends_mk_1979_2014_DJF <- plot_era5_intensity_quantile_trends_mk(era5_min_slp_1979_2014_DJF, 
                                                                           "min_slp", 
                                                                           quantiles = my.probs)
era5_wind_quant_trends_mk_1979_2014_DJF <- plot_era5_intensity_quantile_trends_mk(era5_max_wind_1979_2014_DJF, 
                                                                            "max_wind", 
                                                                            quantiles = my.probs)
# MA
era5_slp_quant_trends_mk_1979_2014_MA <- plot_era5_intensity_quantile_trends_mk(era5_min_slp_1979_2014_MA, 
                                                                          "min_slp", 
                                                                          quantiles = my.probs)
era5_wind_quant_trends_mk_1979_2014_MA <- plot_era5_intensity_quantile_trends_mk(era5_max_wind_1979_2014_MA, 
                                                                           "max_wind", 
                                                                           quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-ERA5, 1979-2014

# Full dataset
crcm5_era5_slp_quant_trends_mk_1979_2014 <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_min_slp_1979_2014,
                                                                             "min_slp", 
                                                                             quantiles = my.probs)
crcm5_era5_wind_quant_trends_mk_1979_2014 <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_max_wind_1979_2014, 
                                                                              "max_wind", 
                                                                              quantiles = my.probs)
# SON
crcm5_era5_slp_quant_trends_mk_1979_2014_SON <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_min_slp_1979_2014_SON, 
                                                                                 "min_slp", 
                                                                                 quantiles = my.probs)
crcm5_era5_wind_quant_trends_mk_1979_2014_SON <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_max_wind_1979_2014_SON, 
                                                                                  "max_wind", 
                                                                                  quantiles = my.probs)
# DJF
crcm5_era5_slp_quant_trends_mk_1979_2014_DJF <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_min_slp_1979_2014_DJF, 
                                                                                 "min_slp", 
                                                                                 quantiles = my.probs)
crcm5_era5_wind_quant_trends_mk_1979_2014_DJF <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_max_wind_1979_2014_DJF, 
                                                                                  "max_wind", 
                                                                                  quantiles = my.probs)
# MA
crcm5_era5_slp_quant_trends_mk_1979_2014_MA <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_min_slp_1979_2014_MA, 
                                                                                "min_slp", 
                                                                                quantiles = my.probs)
crcm5_era5_wind_quant_trends_mk_1979_2014_MA <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_max_wind_1979_2014_MA, 
                                                                                 "max_wind", 
                                                                                 quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-CMIP6 ensemble (pooled pilot models)

# 1950-2014
crcm5_cmip6_ensemble_slp_quant_trends_mk_1950_2014 <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_min_slp_1950_2014, 
                                                                                 metric = "min_slp", 
                                                                                 quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_mk_1950_2014 <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_max_wind_1950_2014, 
                                                                                  metric = "max_wind", 
                                                                                  quantiles = my.probs)
# SON
crcm5_cmip6_ensemble_slp_quant_trends_mk_1950_2014_SON <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_min_slp_1950_2014_SON, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_mk_1950_2014_SON <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_max_wind_1950_2014_SON, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# DJF
crcm5_cmip6_ensemble_slp_quant_trends_mk_1950_2014_DJF <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_min_slp_1950_2014_DJF, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_mk_1950_2014_DJF <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_max_wind_1950_2014_DJF, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# MA
crcm5_cmip6_ensemble_slp_quant_trends_mk_1950_2014_MA <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_min_slp_1950_2014_MA, 
                                                                                    metric = "min_slp", 
                                                                                    quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_mk_1950_2014_MA <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_max_wind_1950_2014_MA, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)

# 1979-2014
crcm5_cmip6_ensemble_slp_quant_trends_mk_1979_2014 <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_min_slp_1979_2014, 
                                                                                 metric = "min_slp", 
                                                                                 quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_mk_1979_2014 <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_max_wind_1979_2014, 
                                                                                  metric = "max_wind", 
                                                                                  quantiles = my.probs)
# SON
crcm5_cmip6_ensemble_slp_quant_trends_mk_1979_2014_SON <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_min_slp_1979_2014_SON, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_mk_1979_2014_SON <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_max_wind_1979_2014_SON, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# DJF
crcm5_cmip6_ensemble_slp_quant_trends_mk_1979_2014_DJF <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_min_slp_1979_2014_DJF, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_mk_1979_2014_DJF <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_max_wind_1979_2014_DJF, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)
# MA
crcm5_cmip6_ensemble_slp_quant_trends_mk_1979_2014_MA <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_min_slp_1979_2014_MA, 
                                                                                    metric = "min_slp", 
                                                                                    quantiles = my.probs)
crcm5_cmip6_ensemble_wind_quant_trends_mk_1979_2014_MA <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_max_wind_1979_2014_MA, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-CMIP6 individual models (saved in data frames)

# 1950-2014
crcm5_cmip6_models_slp_quant_trends_mk_1950_2014 <- plot_intensity_quantile_trends_mk(crcm5_cmip6_min_slp_1950_2014, 
                                                                                metric = "min_slp", 
                                                                                quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_mk_1950_2014 <- plot_intensity_quantile_trends_mk(crcm5_cmip6_max_wind_1950_2014, 
                                                                                 metric = "max_wind", 
                                                                                 quantiles = my.probs)
# SON
crcm5_cmip6_models_slp_quant_trends_mk_1950_2014_SON <- plot_intensity_quantile_trends_mk(crcm5_cmip6_min_slp_1950_2014_SON, 
                                                                                    metric = "min_slp", 
                                                                                    quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_mk_1950_2014_SON <- plot_intensity_quantile_trends_mk(crcm5_cmip6_max_wind_1950_2014_SON, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)
# DJF
crcm5_cmip6_models_slp_quant_trends_mk_1950_2014_DJF <- plot_intensity_quantile_trends_mk(crcm5_cmip6_min_slp_1950_2014_DJF, 
                                                                                    metric = "min_slp", 
                                                                                    quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_mk_1950_2014_DJF <- plot_intensity_quantile_trends_mk(crcm5_cmip6_max_wind_1950_2014_DJF, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)
# MA
crcm5_cmip6_models_slp_quant_trends_mk_1950_2014_MA <- plot_intensity_quantile_trends_mk(crcm5_cmip6_min_slp_1950_2014_MA, 
                                                                                   metric = "min_slp", 
                                                                                   quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_mk_1950_2014_MA <- plot_intensity_quantile_trends_mk(crcm5_cmip6_max_wind_1950_2014_MA, 
                                                                                    metric = "max_wind", 
                                                                                    quantiles = my.probs)

# 1979-2014
crcm5_cmip6_models_slp_quant_trends_mk_1979_2014 <- plot_intensity_quantile_trends_mk(crcm5_cmip6_min_slp_1979_2014, 
                                                                                metric = "min_slp", 
                                                                                quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_mk_1979_2014 <- plot_intensity_quantile_trends_mk(crcm5_cmip6_max_wind_1979_2014, 
                                                                                 metric = "max_wind", 
                                                                                 quantiles = my.probs)
# SON
crcm5_cmip6_models_slp_quant_trends_mk_1979_2014_SON <- plot_intensity_quantile_trends_mk(crcm5_cmip6_min_slp_1979_2014_SON, 
                                                                                    metric = "min_slp", 
                                                                                    quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_mk_1979_2014_SON <- plot_intensity_quantile_trends_mk(crcm5_cmip6_max_wind_1979_2014_SON, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)
# DJF
crcm5_cmip6_models_slp_quant_trends_mk_1979_2014_DJF <- plot_intensity_quantile_trends_mk(crcm5_cmip6_min_slp_1979_2014_DJF, 
                                                                                    metric = "min_slp", 
                                                                                    quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_mk_1979_2014_DJF <- plot_intensity_quantile_trends_mk(crcm5_cmip6_max_wind_1979_2014_DJF, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)
# MA
crcm5_cmip6_models_slp_quant_trends_mk_1979_2014_MA <- plot_intensity_quantile_trends_mk(crcm5_cmip6_min_slp_1979_2014_MA, 
                                                                                   metric = "min_slp", 
                                                                                   quantiles = my.probs)
crcm5_cmip6_models_wind_quant_trends_mk_1979_2014_MA <- plot_intensity_quantile_trends_mk(crcm5_cmip6_max_wind_1979_2014_MA, 
                                                                                    metric = "max_wind", 
                                                                                    quantiles = my.probs)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Calculate intensity quantiles & Mann-Kendall trends over time for each MPI member, 1979-2014

mpi_slp_quant_trends_mk_1979_2014 <- plot_intensity_quantile_trends_mk(mpi_min_slp_1979_2014, 
                                                                 "min_slp", 
                                                                 quantiles = my.probs)
mpi_wind_quant_trends_mk_1979_2014 <- plot_intensity_quantile_trends_mk(mpi_max_wind_1979_2014, 
                                                                  "max_wind", 
                                                                  quantiles = my.probs)

# ==============================================================================

# Calculate intensification rate quantile trends using Mann-Kendall

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CFSR, 1979-2014

# Full dataset
cfsr_slp_intens_quant_trends_mk_1979_2014 <- plot_era5_intensification_quantile_trends_mk(cfsr_intens_rates_slp_1979_2014, 
                                                                                    "min_slp", 
                                                                                    quantiles = my.probs)
cfsr_wind_intens_quant_trends_mk_1979_2014 <- plot_era5_intensification_quantile_trends_mk(cfsr_intens_rates_wind_1979_2014, 
                                                                                     "max_wind", 
                                                                                     quantiles = my.probs)
# SON
cfsr_slp_intens_quant_trends_mk_1979_2014_SON <- plot_era5_intensification_quantile_trends_mk(cfsr_intens_rates_slp_1979_2014_SON, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
cfsr_wind_intens_quant_trends_mk_1979_2014_SON <- plot_era5_intensification_quantile_trends_mk(cfsr_intens_rates_wind_1979_2014_SON, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# DJF
cfsr_slp_intens_quant_trends_mk_1979_2014_DJF <- plot_era5_intensification_quantile_trends_mk(cfsr_intens_rates_slp_1979_2014_DJF, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
cfsr_wind_intens_quant_trends_mk_1979_2014_DJF <- plot_era5_intensification_quantile_trends_mk(cfsr_intens_rates_wind_1979_2014_DJF, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# MA
cfsr_slp_intens_quant_trends_mk_1979_2014_MA <- plot_era5_intensification_quantile_trends_mk(cfsr_intens_rates_slp_1979_2014_MA, 
                                                                                       "min_slp", 
                                                                                       quantiles = my.probs)
cfsr_wind_intens_quant_trends_mk_1979_2014_MA <- plot_era5_intensification_quantile_trends_mk(cfsr_intens_rates_wind_1979_2014_MA, 
                                                                                        "max_wind", 
                                                                                        quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# ERA5

# 1950-2014
era5_slp_intens_quant_trends_mk_1950_2014 <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_slp_1950_2014, 
                                                                                    "min_slp", 
                                                                                    quantiles = my.probs)
era5_wind_intens_quant_trends_mk_1950_2014 <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_wind_1950_2014, 
                                                                                     "max_wind", 
                                                                                     quantiles = my.probs)
# SON
era5_slp_intens_quant_trends_mk_1950_2014_SON <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_slp_1950_2014_SON, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
era5_wind_intens_quant_trends_mk_1950_2014_SON <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_wind_1950_2014_SON, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# DJF
era5_slp_intens_quant_trends_mk_1950_2014_DJF <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_slp_1950_2014_DJF, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
era5_wind_intens_quant_trends_mk_1950_2014_DJF <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_wind_1950_2014_DJF, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# MA
era5_slp_intens_quant_trends_mk_1950_2014_MA <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_slp_1950_2014_MA, 
                                                                                       "min_slp", 
                                                                                       quantiles = my.probs)
era5_wind_intens_quant_trends_mk_1950_2014_MA <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_wind_1950_2014_MA, 
                                                                                        "max_wind", 
                                                                                        quantiles = my.probs)

# 1979-2014
era5_slp_intens_quant_trends_mk_1979_2014 <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_slp_1979_2014, 
                                                                                    "min_slp", 
                                                                                    quantiles = my.probs)
era5_wind_intens_quant_trends_mk_1979_2014 <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_wind_1979_2014, 
                                                                                     "max_wind", 
                                                                                     quantiles = my.probs)
# SON
era5_slp_intens_quant_trends_mk_1979_2014_SON <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_slp_1979_2014_SON, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
era5_wind_intens_quant_trends_mk_1979_2014_SON <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_wind_1979_2014_SON, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# DJF
era5_slp_intens_quant_trends_mk_1979_2014_DJF <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_slp_1979_2014_DJF, 
                                                                                        "min_slp", 
                                                                                        quantiles = my.probs)
era5_wind_intens_quant_trends_mk_1979_2014_DJF <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_wind_1979_2014_DJF, 
                                                                                         "max_wind", 
                                                                                         quantiles = my.probs)
# MA
era5_slp_intens_quant_trends_mk_1979_2014_MA <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_slp_1979_2014_MA, 
                                                                                       "min_slp", 
                                                                                       quantiles = my.probs)
era5_wind_intens_quant_trends_mk_1979_2014_MA <- plot_era5_intensification_quantile_trends_mk(era5_intens_rates_wind_1979_2014_MA, 
                                                                                        "max_wind", 
                                                                                        quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-ERA5, 1979-2014

# Full dataset
crcm5_era5_slp_intens_quant_trends_mk_1979_2014 <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_intens_rates_slp_1979_2014$era5, 
                                                                                          "min_slp", 
                                                                                          quantiles = my.probs)
crcm5_era5_wind_intens_quant_trends_mk_1979_2014 <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_intens_rates_wind_1979_2014$era5, 
                                                                                           "max_wind", 
                                                                                           quantiles = my.probs)
# SON
crcm5_era5_slp_intens_quant_trends_mk_1979_2014_SON <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_intens_rates_slp_1979_2014_SON, 
                                                                                              "min_slp", 
                                                                                              quantiles = my.probs)
crcm5_era5_wind_intens_quant_trends_mk_1979_2014_SON <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_intens_rates_wind_1979_2014_SON, 
                                                                                               "max_wind", 
                                                                                               quantiles = my.probs)
# DJF
crcm5_era5_slp_intens_quant_trends_mk_1979_2014_DJF <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_intens_rates_slp_1979_2014_DJF, 
                                                                                              "min_slp", 
                                                                                              quantiles = my.probs)
crcm5_era5_wind_intens_quant_trends_mk_1979_2014_DJF <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_intens_rates_wind_1979_2014_DJF, 
                                                                                               "max_wind", 
                                                                                               quantiles = my.probs)
# MA
crcm5_era5_slp_intens_quant_trends_mk_1979_2014_MA <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_intens_rates_slp_1979_2014_MA, 
                                                                                             "min_slp", 
                                                                                             quantiles = my.probs)
crcm5_era5_wind_intens_quant_trends_mk_1979_2014_MA <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_intens_rates_wind_1979_2014_MA, 
                                                                                              "max_wind", 
                                                                                              quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-CMIP6 ensemble (pooled pilot models)

# 1950-2014
crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1950_2014 <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1950_2014,  
                                                                                                        quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1950_2014 <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1950_2014, 
                                                                                                         quantiles = my.probs)
# SON
crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1950_2014_SON <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1950_2014_SON,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1950_2014_SON <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1950_2014_SON, 
                                                                                                             quantiles = my.probs)
# DJF
crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1950_2014_DJF <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1950_2014_DJF,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1950_2014_DJF <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1950_2014_DJF, 
                                                                                                             quantiles = my.probs)
# MA
crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1950_2014_MA <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1950_2014_MA,  
                                                                                                           quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1950_2014_MA <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1950_2014_MA, 
                                                                                                            quantiles = my.probs)

# 1979-2014
crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1979_2014 <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1979_2014, 
                                                                                                        quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1979_2014 <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1979_2014, 
                                                                                                         quantiles = my.probs)
# SON
crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1979_2014_SON <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1979_2014_SON,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1979_2014_SON <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1979_2014_SON, 
                                                                                                             quantiles = my.probs)
# DJF
crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1979_2014_DJF <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1979_2014_DJF,  
                                                                                                            quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1979_2014_DJF <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1979_2014_DJF, 
                                                                                                             quantiles = my.probs)
# MA
crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1979_2014_MA <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1979_2014_MA,  
                                                                                                           quantiles = my.probs)
crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1979_2014_MA <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1979_2014_MA, 
                                                                                                            quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# CRCM5-CMIP6 intensification rate quantiles by model

# 1950-2014
crcm5_cmip6_models_slp_intens_quant_trends_mk_1950_2014 <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1950_2014,  
                                                                                             quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_mk_1950_2014 <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1950_2014, 
                                                                                              quantiles = my.probs)
# SON
crcm5_cmip6_models_slp_intens_quant_trends_mk_1950_2014_SON <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1950_2014_SON,  
                                                                                                 quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_mk_1950_2014_SON <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1950_2014_SON, 
                                                                                                  quantiles = my.probs)
# DJF
crcm5_cmip6_models_slp_intens_quant_trends_mk_1950_2014_DJF <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1950_2014_DJF,  
                                                                                                 quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_mk_1950_2014_DJF <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1950_2014_DJF, 
                                                                                                  quantiles = my.probs)
# MA
crcm5_cmip6_models_slp_intens_quant_trends_mk_1950_2014_MA <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1950_2014_MA,  
                                                                                                quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_mk_1950_2014_MA <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1950_2014_MA, 
                                                                                                 quantiles = my.probs)

# 1979-2014
crcm5_cmip6_models_slp_intens_quant_trends_mk_1979_2014 <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1979_2014, 
                                                                                             quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_mk_1979_2014 <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1979_2014, 
                                                                                              quantiles = my.probs)
# SON
crcm5_cmip6_models_slp_intens_quant_trends_mk_1979_2014_SON <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1979_2014_SON,  
                                                                                                 quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_mk_1979_2014_SON <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1979_2014_SON, 
                                                                                                  quantiles = my.probs)
# DJF
crcm5_cmip6_models_slp_intens_quant_trends_mk_1979_2014_DJF <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1979_2014_DJF,  
                                                                                                 quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_mk_1979_2014_DJF <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1979_2014_DJF, 
                                                                                                  quantiles = my.probs)
# MA
crcm5_cmip6_models_slp_intens_quant_trends_mk_1979_2014_MA <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_slp_1979_2014_MA,  
                                                                                                quantiles = my.probs)
crcm5_cmip6_models_wind_intens_quant_trends_mk_1979_2014_MA <- plot_intensification_quantile_trends_mk(crcm5_cmip6_intens_rates_wind_1979_2014_MA, 
                                                                                                 quantiles = my.probs)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Calculate intensification rate quantiles & Mann-Kendall trends over timefor each MPI 
#   member, 1979-2014

mpi_slp_intens_quant_trends_mk_1979_2014 <- plot_intensification_quantile_trends_mk(mpi_intens_rates_slp_1979_2014, 
                                                                              "min_slp", 
                                                                              quantiles = my.probs)
mpi_wind_intens_quant_trends_mk_1979_2014 <- plot_intensification_quantile_trends_mk(mpi_intens_rates_wind_1979_2014, 
                                                                               "max_wind", 
                                                                               quantiles = my.probs)

# ==============================================================================

# Mann-Kendall statistical significance figures (S8-12 & S14-19, *S7 & S13 are
#   TC filter figures and plotted in the TC filter section of this document)

# ==============================================================================

# M-K FIGURE S8

# Percentile vs. p-values intensity fig.
#   - min. SLP in top row, max. wind in bottom row
#   - 1950-2014 on left (ERA5 & CRCM5-CMIP6), 1979-2014 on right (all datasets)
#   - one line for each dataset on each plot

pdf("/Path/to/output/figures/fig_s8_intensity_mk_p_values.pdf", 
    family = "Helvetica-Narrow", width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), oma = c(3, 6, 6, 1))

# 1950-2014 min. SLP
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_mk_1950_2014$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_mk_1950_2014$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("topright", inset = c(-1.15, -0.325), xpd = "True", 
       legend = c("CFSR", "ERA5", "CRCM5-ERA5", "CRCM5-CMIP6"), 
       fill = c("#E69F00", "#00A9A5", "#CC79A7", "#882255"), cex = 1.75, ncol = 4, bty = "n")
mtext("Min. SLP", 2, 6.5, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.25, family = 'Helvetica', font = 2, cex = 2)
mtext("1950-2014", 3, 1.4, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.025, labels = "a", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1979-2014 min. SLP
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_mk_1979_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#00A9A5")
lines(rev(my.probs), crcm5_era5_slp_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#CC79A7")
lines(rev(my.probs), cfsr_slp_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("1979-2014", 3, 1.4, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.025, labels = "b", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1950-2014 max. wind
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_mk_1950_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_mk_1950_2014$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.25, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind", 2, 6.5, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.25, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "c", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1979-2014 max. wind
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_mk_1979_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.25, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "d", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# M-K FIGURE S9

# 1979-2014 seasonal intensity quantiles vs. p-values

pdf("/Path/to/output/figures/fig_s9_seasonal_intensity_mk_p_values_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(3, 3, 1, 1), oma = c(3, 7, 7, 1))

# min. SLP SON
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_mk_1979_2014_SON$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#00A9A5")
lines(rev(my.probs), crcm5_era5_slp_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#CC79A7")
lines(rev(my.probs), cfsr_slp_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Min. SLP", 2, 6.75, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Autumn (SON)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.03, labels = "a", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# min. SLP DJF
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_mk_1979_2014_DJF$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
lines(rev(my.probs), crcm5_era5_slp_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#CC79A7")
lines(rev(my.probs), cfsr_slp_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.285), xpd = "True", legend = c("CFSR", "ERA5", "CRCM5-ERA5", "CRCM5-CMIP6"), 
       fill = c("#E69F00", "#00A9A5", "#CC79A7", "#882255"), cex = 2.5, ncol = 4, bty = "n")
mtext("Winter (DJF)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.03, labels = "c", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# min. SLP MA
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_mk_1979_2014_MA$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#00A9A5")
lines(rev(my.probs), crcm5_era5_slp_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#CC79A7")
lines(rev(my.probs), cfsr_slp_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Spring (MA)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.03, labels = "e", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind SON
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_mk_1979_2014_SON$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind", 2, 6.75, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.03, labels = "b", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind DJF
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_mk_1979_2014_DJF$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.03, labels = "d", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind MA
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_mk_1979_2014_MA$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.03, labels = "f", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# M-K FIGURE S10

# Percentile vs. p-values intensification rate fig.
#   - min. SLP in top row, max. wind in bottom row
#   - 1950-2014 on left (ERA5 & CRCM5-CMIP6), 1979-2014 on right (all datasets)
#   - one line for each dataset on each plot

pdf("/Path/to/output/figures/fig_s10_intensification_mk_p_values.pdf", 
    family = "Helvetica-Narrow", width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), oma = c(3, 6, 6, 1))

# 1950-2014 min. SLP
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1950_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_mk_1950_2014$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("topright", inset = c(-1.15, -0.325), xpd = "True", 
       legend = c("CFSR", "ERA5", "CRCM5-ERA5", "CRCM5-CMIP6"), 
       fill = c("#E69F00", "#00A9A5", "#CC79A7","#882255"), cex = 1.75, ncol = 4, bty = "n")
mtext("Min. SLP Rate", 2, 6.5, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.25, family = 'Helvetica', font = 2, cex = 2)
mtext("1950-2014", 3, 1.4, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.035, labels = "a", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1979-2014 min. SLP
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1979_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_slp_intens_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_slp_intens_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("1979-2014", 3, 1.4, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.035, labels = "b", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1950-2014 max. wind
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1950_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_mk_1950_2014$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.25, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind Rate", 2, 6.5, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.25, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.035, labels = "c", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1979-2014 max. wind
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1979_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_intens_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_intens_quant_trends_mk_1979_2014$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.25, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.035, labels = "d", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# M-K FIGURE S11

# 1979-2014 seasonal intensification rate percentiles vs. p-values

pdf("/Path/to/output/figures/fig_s11_seasonal_intensification_mk_p_values_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(3, 3, 1, 1), oma = c(3, 7, 7, 1))

# min. SLP SON
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1979_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_slp_intens_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_slp_intens_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Min. SLP Rate", 2, 6.75, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Autumn (SON)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.025, labels = "a", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# min. SLP DJF
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1979_2014_DJF$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_slp_intens_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_slp_intens_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.285), xpd = "True", 
       legend = c("CFSR", "ERA5", "CRCM5-ERA5", "CRCM5-CMIP6"), 
       fill = c("#E69F00", "#00A9A5", "#CC79A7", "#882255"), cex = 2.5, ncol = 4, bty = "n")
mtext("Winter (DJF)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.025, labels = "c", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# min. SLP MA
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1979_2014_MA$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_slp_intens_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_slp_intens_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Spring (MA)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.025, labels = "e", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind SON
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1979_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_intens_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_intens_quant_trends_mk_1979_2014_SON$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind Rate", 2, 6.75, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "b", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind DJF
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1979_2014_DJF$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_intens_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_intens_quant_trends_mk_1979_2014_DJF$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "d", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind MA
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1979_2014_MA$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_intens_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_intens_quant_trends_mk_1979_2014_MA$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "f", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# M-K FIGURE S12

# MPI member survival functions, intensity p-values, & intensification rate p-values, 1979-2014

pdf("/Path/to/output/figures/fig_s12_crcm5_mpi_members_mk_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 12, height = 16)
par(mfrow = c(3, 2), mar = c(6, 4, 4, 1), oma = c(3, 4, 3, 2))

# Min. SLP survival function
plot(mpi_sim1_min_slp_quants_1979_2014, my.probs*100, col = "#1B3A6B", type = "l", 
     main = "", lwd = 3, lty = 1, xlim = rev(c(90000, 98000)), ylim = c(0, 100), 
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(mpi_sim2_min_slp_quants_1979_2014, my.probs*100, col = "#2E6DB4", lwd = 3, lty = 1)
lines(mpi_sim3_min_slp_quants_1979_2014, my.probs*100, col = "#4D9FD6", lwd = 3, lty = 1)
lines(mpi_sim4_min_slp_quants_1979_2014, my.probs*100, col = "#89C4E1", lwd = 3, lty = 1)
lines(mpi_sim5_min_slp_quants_1979_2014, my.probs*100, col = "#C2DFF0", lwd = 3, lty = 1)
leg <- legend("topright", inset = c(-0.9, -0.2), xpd = "True", 
              legend = c("                ", "                ", "                ", "                ", "                "), 
              fill = c("#1B3A6B", "#2E6DB4", "#4D9FD6", "#89C4E1", "#C2DFF0"), cex = 2.5, ncol = 5, bty = "n")
text(leg$text$x, leg$text$y,
     labels = c("r1i1p1f1", "r2i1p1f1", "r3i1p1f1", "r4i1p1f1", "r5i1p1f1"),
     family = "Helvetica", cex = 2.5, adj = 0, xpd = NA)
text(x = 90100, y = 97, labels = "a", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
mtext("Exceedance Probability", 2, 4, family = 'Helvetica', font = 2, cex = 2)
mtext("Minimum SLP (hPa)", 1, 3.75, family = 'Helvetica', font = 2, cex = 2)
axis(side = 1, at = seq(89500, 101000, 500), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 90000, to = 100000, by = 2000), labels = c("900", "920", "940", "960", "980", "1000"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2.5, 
     at = seq(from = 0, to = 100, by = 20), labels = c("0", "20", "40", "60", "80", "100"))
box()

# Max. wind survival function
plot(rev(mpi_sim1_max_wind_quants_1979_2014), my.probs*100, col = "#1B3A6B", type = "l", 
     main = "", lwd = 3, lty = 1, xlim = c(15, 50), ylim = c(0, 100), 
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(mpi_sim2_max_wind_quants_1979_2014), my.probs*100, col = "#2E6DB4", lwd = 3, lty = 1)
lines(rev(mpi_sim3_max_wind_quants_1979_2014), my.probs*100, col = "#4D9FD6", lwd = 3, lty = 1)
lines(rev(mpi_sim4_max_wind_quants_1979_2014), my.probs*100, col = "#89C4E1", lwd = 3, lty = 1)
lines(rev(mpi_sim5_max_wind_quants_1979_2014), my.probs*100, col = "#C2DFF0", lwd = 3, lty = 1)
text(x = 49.5, y = 97, labels = "b", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
mtext("Maximum Wind Speed (m/s)", 1, 3.75, family = 'Helvetica', font = 2, cex = 2)
axis(side = 1, at = seq(0, 60, 2), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 50, by = 10), labels = c("0", "10", "20", "30", "40", "50"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2.5, 
     at = seq(from = 0, to = 100, by = 20), labels = c("0", "20", "40", "60", "80", "100"))
box()

# 1979-2014 min. SLP quantile trends
plot(rev(my.probs), mpi_slp_quant_trends_mk_1979_2014$p_values$sim1, type = "l", 
     col = ("#1B3A6B"), lwd = 3, xlim = c(0, 1), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), mpi_slp_quant_trends_mk_1979_2014$p_values$sim2, lwd = 3, col = "#2E6DB4")
lines(rev(my.probs), mpi_slp_quant_trends_mk_1979_2014$p_values$sim3, lwd = 3, col = "#4D9FD6")
lines(rev(my.probs), mpi_slp_quant_trends_mk_1979_2014$p_values$sim4, lwd = 3, col = "#89C4E1")
lines(rev(my.probs), mpi_slp_quant_trends_mk_1979_2014$p_values$sim5, lwd = 3, col = "#C2DFF0")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
text(x = 0.98, y = 0.97, labels = "c", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
mtext("Percentile", 1, 3.75, family = 'Helvetica', font = 2, cex = 2)
mtext("p-value", 2, 4, family = 'Helvetica', font = 2, cex = 2)
mtext("Min. SLP Intensity", 3, 0.5, family = 'Helvetica', font = 2, cex = 2)
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1979-2014 max. wind quantile trends
plot(my.probs, mpi_wind_quant_trends_mk_1979_2014$p_values$sim1, type = "l", 
     col = ("#1B3A6B"), lwd = 3, xlim = c(0, 1), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, mpi_wind_quant_trends_mk_1979_2014$p_values$sim2, lwd = 3, col = "#2E6DB4")
lines(my.probs, mpi_wind_quant_trends_mk_1979_2014$p_values$sim3, lwd = 3, col = "#4D9FD6")
lines(my.probs, mpi_wind_quant_trends_mk_1979_2014$p_values$sim4, lwd = 3, col = "#89C4E1")
lines(my.probs, mpi_wind_quant_trends_mk_1979_2014$p_values$sim5, lwd = 3, col = "#C2DFF0")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
text(x = 0.98, y = 0.97, labels = "d", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
mtext("Percentile", 1, 3.75, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Wind Speed Intensity", 3, 0.5, family = 'Helvetica', font = 2, cex = 2)
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1979-2014 min. SLP intens. rate quantile trends
plot(my.probs, mpi_slp_intens_quant_trends_mk_1979_2014$p_values$sim1, type = "l", 
     col = ("#1B3A6B"), lwd = 3, xlim = c(0, 1), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, mpi_slp_intens_quant_trends_mk_1979_2014$p_values$sim2, lwd = 3, col = "#2E6DB4")
lines(my.probs, mpi_slp_intens_quant_trends_mk_1979_2014$p_values$sim3, lwd = 3, col = "#4D9FD6")
lines(my.probs, mpi_slp_intens_quant_trends_mk_1979_2014$p_values$sim4, lwd = 3, col = "#89C4E1")
lines(my.probs, mpi_slp_intens_quant_trends_mk_1979_2014$p_values$sim5, lwd = 3, col = "#C2DFF0")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
text(x = 0.98, y = 0.97, labels = "e", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
mtext("Percentile", 1, 3.75, family = 'Helvetica', font = 2, cex = 2)
mtext("p-value", 2, 4, family = 'Helvetica', font = 2, cex = 2)
mtext("Min. SLP Intensification Rate", 3, 0.5, family = 'Helvetica', font = 2, cex = 2)
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1979-2014 max. wind intens. rate quantile trends
plot(my.probs, mpi_wind_intens_quant_trends_mk_1979_2014$p_values$sim1, type = "l", 
     col = ("#1B3A6B"), lwd = 3, xlim = c(0, 1), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, mpi_wind_intens_quant_trends_mk_1979_2014$p_values$sim2, lwd = 3, col = "#2E6DB4")
lines(my.probs, mpi_wind_intens_quant_trends_mk_1979_2014$p_values$sim3, lwd = 3, col = "#4D9FD6")
lines(my.probs, mpi_wind_intens_quant_trends_mk_1979_2014$p_values$sim4, lwd = 3, col = "#89C4E1")
lines(my.probs, mpi_wind_intens_quant_trends_mk_1979_2014$p_values$sim5, lwd = 3, col = "#C2DFF0")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
text(x = 0.98, y = 0.97, labels = "f", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
mtext("Percentile", 1, 3.75, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Wind Intensification Rate", 3, 0.5, family = 'Helvetica', font = 2, cex = 2)
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# M-K FIGURE S14

# Plot intensity quantiles vs. p-values by model: 1950-2014

pdf("/Path/to/output/figures/fig_s14_crcm5_cmip6_pilots_intensity_mk_p_values_1950_2014.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 5, 2))

# CanESM5 min SLP
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_mk_1950_2014$p_values$canesm5, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Min. SLP", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("CRCM5-CanESM5", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "a", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CNRM min SLP
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_mk_1950_2014$p_values$cnrm, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-CNRM", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "b", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# MPI min SLP
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_mk_1950_2014$p_values$mpi, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-MPI", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "c", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# NorESM2 min SLP
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_mk_1950_2014$p_values$noresm2, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-NorESM2", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "d", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CanESM5 max wind
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_mk_1950_2014$p_values$canesm5, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "e", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CNRM max wind
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_mk_1950_2014$p_values$cnrm, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "f", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# MPI max wind
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_mk_1950_2014$p_values$mpi, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "g", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# NorESM2 max wind
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_mk_1950_2014$p_values$noresm2, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# M-K FIGURE S15

# Plot intensity quantiles vs. p-values by model: 1979-2014

pdf("/Path/to/output/figures/fig_s15_crcm5_cmip6_pilots_intensity_mk_p_values_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 5, 2))

# CanESM5 min SLP
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_mk_1979_2014$p_values$canesm5, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Min. SLP", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("CRCM5-CanESM5", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "a", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CNRM min SLP
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_mk_1979_2014$p_values$cnrm, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-CNRM", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "b", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# MPI min SLP
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_mk_1979_2014$p_values$mpi, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-MPI", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "c", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# NorESM2 min SLP
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_mk_1979_2014$p_values$noresm2, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-NorESM2", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "d", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CanESM5 max wind
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_mk_1979_2014$p_values$canesm5, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "e", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CNRM max wind
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_mk_1979_2014$p_values$cnrm, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "f", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# MPI max wind
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_mk_1979_2014$p_values$mpi, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "g", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# NorESM2 max wind
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_mk_1979_2014$p_values$noresm2, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# M-K FIGURE S16

# Plot intensification rate quantiles vs. p-values by model: 1950-2014

pdf("/Path/to/output/figures/fig_s16_crcm5_cmip6_pilots_intensification_mk_p_values_1950_2014.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 5, 2))

# CanESM5 min SLP
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_mk_1950_2014$p_values$canesm5, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Min. SLP Rate", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("CRCM5-CanESM5", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "a", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CNRM min SLP
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_mk_1950_2014$p_values$cnrm, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-CNRM", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "b", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# MPI min SLP
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_mk_1950_2014$p_values$mpi, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-MPI", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "c", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# NorESM2 min SLP
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_mk_1950_2014$p_values$noresm2, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-NorESM2", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "d", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CanESM5 max wind
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_mk_1950_2014$p_values$canesm5, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind Rate", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "e", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CNRM max wind
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_mk_1950_2014$p_values$cnrm, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "f", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# MPI max wind
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_mk_1950_2014$p_values$mpi, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "g", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# NorESM2 max wind
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_mk_1950_2014$p_values$noresm2, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# M-K FIGURE S17

# Plot intensification rate quantiles vs. p-values by model: 1979-2014

pdf("/Path/to/output/figures/fig_s17_crcm5_cmip6_pilots_intensification_mk_p_values_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 5, 2))

# CanESM5 min SLP
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_mk_1979_2014$p_values$canesm5, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Min. SLP Rate", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("CRCM5-CanESM5", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "a", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CNRM min SLP
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_mk_1979_2014$p_values$cnrm, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-CNRM", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "b", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# MPI min SLP
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_mk_1979_2014$p_values$mpi, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-MPI", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "c", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# NorESM2 min SLP
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_mk_1979_2014$p_values$noresm2, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("CRCM5-NorESM2", 3, 1.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "d", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CanESM5 max wind
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_mk_1979_2014$p_values$canesm5, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind Rate", 2, 7.5, family = 'Helvetica', font = 2, cex = 2.5)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "e", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# CNRM max wind
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_mk_1979_2014$p_values$cnrm, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "f", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# MPI max wind
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_mk_1979_2014$p_values$mpi, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "g", cex = 3.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# NorESM2 max wind
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_mk_1979_2014$p_values$noresm2, 
     type = "l", lwd = 3, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# M-K FIGURE S18

# 1950-2014 seasonal intensity quantiles vs. p-values

pdf("/Path/to/output/figures/fig_s18_seasonal_intensity_mk_p_values_1950_2014.pdf", 
    family = "Helvetica-Narrow", width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(3, 3, 1, 1), oma = c(3, 7, 7, 1))

# min. SLP SON
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_mk_1950_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_mk_1950_2014_SON$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Min. SLP", 2, 6.75, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Autumn (SON)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.03, labels = "a", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# min. SLP DJF
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_mk_1950_2014_DJF$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_mk_1950_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.285), xpd = "True", legend = c("ERA5", "CRCM5-CMIP6"), 
       fill = c("#00A9A5", "#882255"), cex = 2.5, ncol = 2, bty = "n")
mtext("Winter (DJF)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.03, labels = "c", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# min. SLP MA
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_mk_1950_2014_MA$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_mk_1950_2014_MA$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Spring (MA)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.03, labels = "e", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind SON
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_mk_1950_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_mk_1950_2014_SON$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind", 2, 6.75, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.03, labels = "b", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind DJF
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_mk_1950_2014_DJF$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_mk_1950_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.03, labels = "d", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind MA
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_mk_1950_2014_MA$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_mk_1950_2014_MA$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.03, labels = "f", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# M-K FIGURE S19

# 1950-2014 seasonal intensification rate quantiles vs. p-values

pdf("/Path/to/output/figures/fig_s19_seasonal_intensification_mk_p_values_1950_2014.pdf", 
    family = "Helvetica-Narrow", width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(3, 3, 1, 1), oma = c(3, 7, 7, 1))

# min. SLP SON
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1950_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_mk_1950_2014_SON$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Min. SLP Rate", 2, 6.75, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Autumn (SON)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.03, labels = "a", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# min. SLP DJF
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1950_2014_DJF$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_mk_1950_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("top", inset = c(0, -0.285), xpd = "True", legend = c("ERA5", "CRCM5-CMIP6"), 
       fill = c("#00A9A5", "#882255"), cex = 2.5, ncol = 2, bty = "n")
mtext("Winter (DJF)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.03, labels = "c", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# min. SLP MA
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_mk_1950_2014_MA$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_mk_1950_2014_MA$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Spring (MA)", 3, 1.5, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.03, labels = "e", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind SON
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1950_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_mk_1950_2014_SON$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind Rate", 2, 6.75, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.03, labels = "b", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind DJF
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1950_2014_DJF$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_mk_1950_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.04, labels = "d", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# max. wind MA
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_mk_1950_2014_MA$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_mk_1950_2014_MA$p_values, lwd = 2, col = "#00A9A5")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.03, labels = "f", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.5, 0), cex.axis = 2.5, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

dev.off()

# ==============================================================================

# Mann-Kendall trend comparison of nor'easter datasets with and without tropical
#   cyclone (TC) filtering

# ==============================================================================

# Load TC filtering workflow
source("/Path/to/analysis/file/TC_filter_hist_noreaster_intensity.R")

# ==============================================================================

# Mann-Kendall intensity quantile testing

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Calculate min. SLP quantile trends (M-K) for all datasets, full time period available for
#   each dataset

# CFSR
cfsr_slp_quant_trends_mk_wTC <- plot_era5_intensity_quantile_trends_mk(cfsr_min_slp_wTC, 
                                                                 metric = "min_slp", 
                                                                 quantiles = my.probs)
cfsr_slp_quant_trends_mk_noTC <- plot_era5_intensity_quantile_trends_mk(cfsr_min_slp_noTC, 
                                                                  metric = "min_slp", 
                                                                  quantiles = my.probs)

# ERA5
era5_slp_quant_trends_mk_wTC <- plot_era5_intensity_quantile_trends_mk(era5_min_slp_wTC, 
                                                                 metric = "min_slp", 
                                                                 quantiles = my.probs)
era5_slp_quant_trends_mk_noTC <- plot_era5_intensity_quantile_trends_mk(era5_min_slp_noTC, 
                                                                  metric = "min_slp", 
                                                                  quantiles = my.probs)

# CRCM5-ERA5
crcm5_era5_slp_quant_trends_mk_wTC <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_min_slp_wTC, 
                                                                       metric = "min_slp", 
                                                                       quantiles = my.probs)
crcm5_era5_slp_quant_trends_mk_noTC <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_min_slp_noTC, 
                                                                        metric = "min_slp", 
                                                                        quantiles = my.probs)

# CRCM5
crcm5_cmip6_slp_quant_trends_mk_wTC <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_min_slp_wTC, 
                                                                  metric = "min_slp", 
                                                                  quantiles = my.probs)
crcm5_cmip6_slp_quant_trends_mk_noTC <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_min_slp_noTC, 
                                                                   metric = "min_slp", 
                                                                   quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Calculate max. wind quantile trends_mk (M-K) for all datasets, full time period available for
#   each dataset

# CFSR
cfsr_wind_quant_trends_mk_wTC <- plot_era5_intensity_quantile_trends_mk(cfsr_max_wind_wTC, 
                                                                  metric = "max_wind", 
                                                                  quantiles = my.probs)
cfsr_wind_quant_trends_mk_noTC <- plot_era5_intensity_quantile_trends_mk(cfsr_max_wind_noTC, 
                                                                   metric = "max_wind", 
                                                                   quantiles = my.probs)

# ERA5
era5_wind_quant_trends_mk_wTC <- plot_era5_intensity_quantile_trends_mk(era5_max_wind_wTC, 
                                                                  metric = "max_wind", 
                                                                  quantiles = my.probs)
era5_wind_quant_trends_mk_noTC <- plot_era5_intensity_quantile_trends_mk(era5_max_wind_noTC, 
                                                                   metric = "max_wind", 
                                                                   quantiles = my.probs)

# CRCM5-ERA5
crcm5_era5_wind_quant_trends_mk_wTC <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_max_wind_wTC, 
                                                                        metric = "max_wind", 
                                                                        quantiles = my.probs)
crcm5_era5_wind_quant_trends_mk_noTC <- plot_era5_intensity_quantile_trends_mk(crcm5_era5_max_wind_noTC, 
                                                                         metric = "max_wind", 
                                                                         quantiles = my.probs)

# CRCM5
crcm5_cmip6_wind_quant_trends_mk_wTC <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_max_wind_wTC, 
                                                                   metric = "max_wind", 
                                                                   quantiles = my.probs)
crcm5_cmip6_wind_quant_trends_mk_noTC <- plot_ensemble_quantile_trends_mk(crcm5_cmip6_max_wind_noTC, 
                                                                    metric = "max_wind", 
                                                                    quantiles = my.probs)

# ==============================================================================

# Mann-Kendall intensification rate quantile testing

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Calculate min. SLP rate quantile trends (M-K) for all datasets, full time period available
#   for each dataset

# CFSR
cfsr_slp_intens_quant_trends_mk_wTC <- plot_era5_intensification_quantile_trends_mk(cfsr_slp_intens_rates_wTC, 
                                                                              metric = "min_slp", 
                                                                              quantiles = my.probs)
cfsr_slp_intens_quant_trends_mk_noTC <- plot_era5_intensification_quantile_trends_mk(cfsr_slp_intens_rates_noTC, 
                                                                               metric = "min_slp", 
                                                                               quantiles = my.probs)

# ERA5
era5_slp_intens_quant_trends_mk_wTC <- plot_era5_intensification_quantile_trends_mk(era5_slp_intens_rates_wTC, 
                                                                              metric = "min_slp", 
                                                                              quantiles = my.probs)
era5_slp_intens_quant_trends_mk_noTC <- plot_era5_intensification_quantile_trends_mk(era5_slp_intens_rates_noTC, 
                                                                               metric = "min_slp", 
                                                                               quantiles = my.probs)

# CRCM5-ERA5
crcm5_era5_slp_intens_quant_trends_mk_wTC <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_slp_intens_rates_wTC, 
                                                                                    metric = "min_slp", 
                                                                                    quantiles = my.probs)
crcm5_era5_slp_intens_quant_trends_mk_noTC <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_slp_intens_rates_noTC, 
                                                                                     metric = "min_slp", 
                                                                                     quantiles = my.probs)

# CRCM5
crcm5_cmip6_slp_intens_quant_trends_mk_wTC <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_slp_intens_rates_wTC, 
                                                                                         quantiles = my.probs)
crcm5_cmip6_slp_intens_quant_trends_mk_noTC <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_slp_intens_rates_noTC, 
                                                                                          quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Calculate max. wind rate quantile trends_mk (M-K) for all datasets, full time period available
#   for each dataset

# CFSR
cfsr_wind_intens_quant_trends_mk_wTC <- plot_era5_intensification_quantile_trends_mk(cfsr_wind_intens_rates_wTC, 
                                                                               metric = "max_wind", 
                                                                               quantiles = my.probs)
cfsr_wind_intens_quant_trends_mk_noTC <- plot_era5_intensification_quantile_trends_mk(cfsr_wind_intens_rates_noTC, 
                                                                                metric = "max_wind", 
                                                                                quantiles = my.probs)

# ERA5
era5_wind_intens_quant_trends_mk_wTC <- plot_era5_intensification_quantile_trends_mk(era5_wind_intens_rates_wTC, 
                                                                               metric = "max_wind", 
                                                                               quantiles = my.probs)
era5_wind_intens_quant_trends_mk_noTC <- plot_era5_intensification_quantile_trends_mk(era5_wind_intens_rates_noTC, 
                                                                                metric = "max_wind", 
                                                                                quantiles = my.probs)

# CRCM5-ERA5
crcm5_era5_wind_intens_quant_trends_mk_wTC <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_wind_intens_rates_wTC, 
                                                                                     metric = "max_wind", 
                                                                                     quantiles = my.probs)
crcm5_era5_wind_intens_quant_trends_mk_noTC <- plot_era5_intensification_quantile_trends_mk(crcm5_era5_wind_intens_rates_noTC, 
                                                                                      metric = "max_wind", 
                                                                                      quantiles = my.probs)

# CRCM5
crcm5_cmip6_wind_intens_quant_trends_mk_wTC <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_wind_intens_rates_wTC, 
                                                                                          quantiles = my.probs)
crcm5_cmip6_wind_intens_quant_trends_mk_noTC <- plot_ensemble_intensification_quantile_trends_mk(crcm5_cmip6_wind_intens_rates_noTC, 
                                                                                           quantiles = my.probs)

# ==============================================================================

# Mann-Kendall trend versions of Main text Figure 1 & Supplemental Figure S10

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# M-K FIGURE S7

# Percentile vs. p-values intensity fig. with & w/o TCs
#   - min. SLP in top row, max. wind in bottom row
#   - each column is a different dataset

pdf("/Path/to/output/figures/fig_s7_TC_intensity_mk_p_values.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 7, 2))

# CFSR min. SLP
plot(rev(my.probs), cfsr_slp_quant_trends_mk_wTC$p_values, type = "l", col = ("#E69F00"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), cfsr_slp_quant_trends_mk_noTC$p_values, lwd = 3, col = "#E69F00")
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
plot(rev(my.probs), era5_slp_quant_trends_mk_wTC$p_values, type = "l", col = ("#00A9A5"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_mk_noTC$p_values, lwd = 3, col = "#00A9A5")
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
plot(rev(my.probs), crcm5_era5_slp_quant_trends_mk_wTC$p_values, type = "l", col = ("#CC79A7"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), crcm5_era5_slp_quant_trends_mk_noTC$p_values, lwd = 3, col = "#CC79A7")
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
plot(rev(my.probs), crcm5_cmip6_slp_quant_trends_mk_wTC$p_values, type = "l", col = ("#882255"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), crcm5_cmip6_slp_quant_trends_mk_noTC$p_values, lwd = 3, col = "#882255")
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
plot(my.probs, cfsr_wind_quant_trends_mk_wTC$p_values, type = "l", col = ("#E69F00"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, cfsr_wind_quant_trends_mk_noTC$p_values, lwd = 3, col = "#E69F00")
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
plot(my.probs, era5_wind_quant_trends_mk_wTC$p_values, type = "l", col = ("#00A9A5"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_mk_noTC$p_values, lwd = 3, col = "#00A9A5")
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
plot(my.probs, crcm5_era5_wind_quant_trends_mk_wTC$p_values, type = "l", col = ("#CC79A7"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_era5_wind_quant_trends_mk_noTC$p_values, lwd = 3, col = "#CC79A7")
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
plot(my.probs, crcm5_cmip6_wind_quant_trends_mk_wTC$p_values, type = "l", col = ("#882255"), 
     lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_cmip6_wind_quant_trends_mk_noTC$p_values, lwd = 3, col = "#882255")
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

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# M-K FIGURE S13

# CRCM5-CMIP6. plot quantile vs. p-values intensification fig.
#   - min. SLP in top row, max. wind in bottom row
#   - each column is a different dataset

pdf("/Path/to/output/figures/fig_s13_TC_intensification_mk_p_values.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 7, 2))

# CFSR min. SLP
plot(my.probs, cfsr_slp_intens_quant_trends_mk_wTC$p_values, type = "l", 
     col = ("#E69F00"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, cfsr_slp_intens_quant_trends_mk_noTC$p_values, lwd = 3, col = "#E69F00")
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
plot(my.probs, era5_slp_intens_quant_trends_mk_wTC$p_values, type = "l", 
     col = ("#00A9A5"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_mk_noTC$p_values, lwd = 3, col = "#00A9A5")
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
plot(my.probs, crcm5_era5_slp_intens_quant_trends_mk_wTC$p_values, type = "l", 
     col = ("#CC79A7"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_era5_slp_intens_quant_trends_mk_noTC$p_values, lwd = 3, col = "#CC79A7")
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
plot(my.probs, crcm5_cmip6_slp_intens_quant_trends_mk_wTC$p_values, type = "l", 
     col = ("#882255"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_cmip6_slp_intens_quant_trends_mk_noTC$p_values, lwd = 3, col = "#882255")
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
plot(my.probs, cfsr_wind_intens_quant_trends_mk_wTC$p_values, type = "l", 
     col = ("#E69F00"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, cfsr_wind_intens_quant_trends_mk_noTC$p_values, lwd = 3, col = "#E69F00")
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
plot(my.probs, era5_wind_intens_quant_trends_mk_wTC$p_values, type = "l", 
     col = ("#00A9A5"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_mk_noTC$p_values, lwd = 3, col = "#00A9A5")
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
plot(my.probs, crcm5_era5_wind_intens_quant_trends_mk_wTC$p_values, type = "l", 
     col = ("#CC79A7"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_era5_wind_intens_quant_trends_mk_noTC$p_values, lwd = 3, col = "#CC79A7")
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
plot(my.probs, crcm5_cmip6_wind_intens_quant_trends_mk_wTC$p_values, type = "l", 
     col = ("#882255"), lwd = 3, lty = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, crcm5_cmip6_wind_intens_quant_trends_mk_noTC$p_values, lwd = 3, col = "#882255")
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

# CFSR & ERA5 intensity & intensification rate trend analysis: 1950-present & 
#   1979-present

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# Filter CFSR & ERA5 storm datasets to appropriate year ranges

# CFSR data starts in 1979, so no analysis of 1950-2025 trends

# Filter CFSR storms for 1979-2025
cfsr_strong_noreasters_1979_2025 <- filter_years_era5(cfsr_strong_noreasters,
                                                      start_year = 1980,  # CFSR nor'easter seasons are
                                                      end_year = 2025)    # keyed to the year they end


# Filter ERA5 storms for 1950-2025
era5_strong_noreasters_1950_2025 <- filter_years_era5(era5_strong_noreasters,
                                                      start_year = 1950,  # ERA5 nor'easter seasons are
                                                      end_year = 2024)    # keyed to the year they begin
# Filter ERA5 storms for 1979-2025
era5_strong_noreasters_1979_2025 <- filter_years_era5(era5_strong_noreasters,
                                                      start_year = 1979,  # ERA5 nor'easter seasons are
                                                      end_year = 2024)    # keyed to the year they begin

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# Extract intensity values for each dataset/time interval/intensity metric

# CFSR
cfsr_min_slp_1979_2025 <- process_era5_for_intensity(cfsr_strong_noreasters_1979_2025, 
                                                     metric = "min_slp")
cfsr_max_wind_1979_2025 <- process_era5_for_intensity(cfsr_strong_noreasters_1979_2025, 
                                                      metric = "max_wind")

# ERA5
era5_min_slp_1950_2025 <- process_era5_for_intensity(era5_strong_noreasters_1950_2025, 
                                                     metric = "min_slp")
era5_max_wind_1950_2025 <- process_era5_for_intensity(era5_strong_noreasters_1950_2025, 
                                                      metric = "max_wind")
era5_min_slp_1979_2025 <- process_era5_for_intensity(era5_strong_noreasters_1979_2025, 
                                                     metric = "min_slp")
era5_max_wind_1979_2025 <- process_era5_for_intensity(era5_strong_noreasters_1979_2025, 
                                                      metric = "max_wind")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# Calculate intensity trends

# CFSR
cfsr_slp_quant_trends_1979_2025 <- plot_era5_intensity_quantile_trends(cfsr_min_slp_1979_2025, # plots time series of each quantile
                                                                       "min_slp", 
                                                                       quantiles = my.probs)
cfsr_wind_quant_trends_1979_2025 <- plot_era5_intensity_quantile_trends(cfsr_max_wind_1979_2025, 
                                                                        "max_wind", 
                                                                        quantiles = my.probs)

# ERA5
era5_slp_quant_trends_1950_2025 <- plot_era5_intensity_quantile_trends(era5_min_slp_1950_2025,
                                                                       "min_slp", 
                                                                       quantiles = my.probs)
era5_wind_quant_trends_1950_2025 <- plot_era5_intensity_quantile_trends(era5_max_wind_1950_2025, 
                                                                        "max_wind", 
                                                                        quantiles = my.probs)
era5_slp_quant_trends_1979_2025 <- plot_era5_intensity_quantile_trends(era5_min_slp_1979_2025,
                                                                       "min_slp", 
                                                                       quantiles = my.probs)
era5_wind_quant_trends_1979_2025 <- plot_era5_intensity_quantile_trends(era5_max_wind_1979_2025, 
                                                                        "max_wind", 
                                                                        quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# Calculate intensification rates

# CFSR
cfsr_intens_rates_slp_1979_2025 <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_1979_2025, 
                                                                           metric ="min_slp", 
                                                                           time_window_hours = 24, 
                                                                           timestep_hours = 6)
cfsr_intens_rates_wind_1979_2025 <- calculate_era5_max_intensification_rate(cfsr_strong_noreasters_1979_2025, 
                                                                            metric ="max_wind", 
                                                                            time_window_hours = 24, 
                                                                            timestep_hours = 6)

# ERA5
era5_intens_rates_slp_1950_2025 <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1950_2025, 
                                                                           metric ="min_slp", 
                                                                           time_window_hours = 24, 
                                                                           timestep_hours = 1)
era5_intens_rates_wind_1950_2025 <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1950_2025, 
                                                                            metric ="max_wind", 
                                                                            time_window_hours = 24, 
                                                                            timestep_hours = 1)
era5_intens_rates_slp_1979_2025 <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1979_2025, 
                                                                           metric ="min_slp", 
                                                                           time_window_hours = 24, 
                                                                           timestep_hours = 1)
era5_intens_rates_wind_1979_2025 <- calculate_era5_max_intensification_rate(era5_strong_noreasters_1979_2025, 
                                                                            metric ="max_wind", 
                                                                            time_window_hours = 24, 
                                                                            timestep_hours = 1)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# Calculate intensification rate trends

# CFSR
cfsr_slp_intens_quant_trends_1979_2025 <- plot_era5_intensification_quantile_trends(cfsr_intens_rates_slp_1979_2025, 
                                                                                    "min_slp", 
                                                                                    quantiles = my.probs)
cfsr_wind_intens_quant_trends_1979_2025 <- plot_era5_intensification_quantile_trends(cfsr_intens_rates_wind_1979_2025, 
                                                                                     "max_wind", 
                                                                                     quantiles = my.probs)

# ERA5
era5_slp_intens_quant_trends_1950_2025 <- plot_era5_intensification_quantile_trends(era5_intens_rates_slp_1950_2025, 
                                                                                    "min_slp", 
                                                                                    quantiles = my.probs)
era5_wind_intens_quant_trends_1950_2025 <- plot_era5_intensification_quantile_trends(era5_intens_rates_wind_1950_2025, 
                                                                                     "max_wind", 
                                                                                     quantiles = my.probs)
era5_slp_intens_quant_trends_1979_2025 <- plot_era5_intensification_quantile_trends(era5_intens_rates_slp_1979_2025, 
                                                                                    "min_slp", 
                                                                                    quantiles = my.probs)
era5_wind_intens_quant_trends_1979_2025 <- plot_era5_intensification_quantile_trends(era5_intens_rates_wind_1979_2025, 
                                                                                     "max_wind", 
                                                                                     quantiles = my.probs)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# FIGURE S29: intensity figure (2x2, time period columns, intensity metric rows)

pdf("/Path/to/output/figures/fig_s29_present_intensity_p_values.pdf", 
    family = "Helvetica-Narrow", width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), oma = c(3, 6, 6, 1))

# 1950-2025 min. SLP
plot(rev(my.probs), era5_slp_quant_trends_1950_2025$p_values, 
     type = "l", col = ("#00A9A5"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("topright", inset = c(-1.15, -0.325), xpd = "True", 
       legend = c("CFSR", "ERA5"), 
       fill = c("#E69F00", "#00A9A5"), cex = 1.75, ncol = 4, bty = "n")
mtext("Min. SLP", 2, 6.5, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.25, family = 'Helvetica', font = 2, cex = 2)
mtext("1950-2025", 3, 1.4, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.025, labels = "a", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1979-2025 min. SLP
plot(rev(my.probs), era5_slp_quant_trends_1979_2025$p_values, type = "l", 
     col = ("#00A9A5"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), cfsr_slp_quant_trends_1979_2025$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("1979-2025", 3, 1.4, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.025, labels = "b", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1950-2025 max. wind
plot(my.probs, era5_wind_quant_trends_1950_2025$p_values, type = "l", 
     col = ("#00A9A5"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.25, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind", 2, 6.5, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.25, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "c", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1979-2025 max. wind
plot(my.probs, era5_wind_quant_trends_1979_2025$p_values, type = "l", 
     col = ("#00A9A5"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, cfsr_wind_quant_trends_1979_2025$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.25, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.025, labels = "d", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

dev.off()

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# FIGURE S30: intensification rate figure (2x2, time period columns, intensity metric rows)

pdf("/Path/to/output/figures/fig_s30_present_intensification_p_values.pdf", 
    family = "Helvetica-Narrow", width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), oma = c(3, 6, 6, 1))

# 1950-2025 min. SLP
plot(my.probs, era5_slp_intens_quant_trends_1950_2025$p_values, type = "l", 
     col = ("#00A9A5"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
legend("topright", inset = c(-1.15, -0.325), xpd = "True", 
       legend = c("CFSR", "ERA5"), 
       fill = c("#E69F00", "#00A9A5"), cex = 1.75, ncol = 4, bty = "n")
mtext("Min. SLP Rate", 2, 6.5, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.25, family = 'Helvetica', font = 2, cex = 2)
mtext("1950-2025", 3, 1.4, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.035, labels = "a", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1979-2025 min. SLP
plot(my.probs, era5_slp_intens_quant_trends_1979_2025$p_values, type = "l", 
     col = ("#00A9A5"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, cfsr_slp_intens_quant_trends_1979_2025$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("1979-2025", 3, 1.4, family = 'Helvetica', font = 2, cex = 2.25)
text(x = 0.02, y = 1.035, labels = "b", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1950-2025 max. wind
plot(my.probs, era5_wind_intens_quant_trends_1950_2025$p_values, type = "l", 
     col = ("#00A9A5"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.25, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Surface Wind Rate", 2, 6.5, family = 'Helvetica', font = 2, cex = 2.25)
mtext("p-value", 2, 3.25, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.035, labels = "c", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

# 1979-2025 max. wind
plot(my.probs, era5_wind_intens_quant_trends_1979_2025$p_values, type = "l", 
     col = ("#00A9A5"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, cfsr_wind_intens_quant_trends_1979_2025$p_values, lwd = 2, col = "#E69F00")
abline(h = c(0.05, 0.1), col = c("gray20", "gray20"), lty = 2, lwd = 1)
mtext("Percentile", 1, 3.25, family = 'Helvetica', font = 2, cex = 2)
text(x = 0.02, y = 1.035, labels = "d", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "20", "40", "60", "80", "100"))
axis(side = 2, at = seq(0, 1, 0.05), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.75, las = 1, 
     at = seq(from = 0, to = 1, by = 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"))

dev.off()




