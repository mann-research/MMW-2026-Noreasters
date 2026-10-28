# ==============================================================================

# Final Figures for Historical Nor'easters: A Cross-Dataset Perspective

# R script to compare CFSR, ERA5, CRCM5-ERA5, & CRCM5-CMIP6 historical intensity
#     characteristics
#     - Fig. 1: intensity survival functions & min. SLP/max. wind comparison for
#       each dataset
#     - Fig. 2 & S22: intensity p-values w/ and w/o TCs included, intensification
#       rate p-values w/ and w/o TCs included (NOT INCLUDED HERE, CODED SEPARATELY
#       IN TC_filter_hist_noreaster_intensity.R)
#     - Fig. 3-6: intensity p-values, seasonal intensity p-values, intensification
#       p-values, seasonal intensification p-values
#     - Fig. 7-9: nor'easter trajectories, locations of peak intensities &
#       intnsification rates
#     - Fig. S1: map of spatial domains for each dataset
#     - Fig. S2: CRCM5-MPI historical members intensity survival functions, 
#       intensity trends, & intensification rate trends
#     - Fig. S3-6: historical tropical cyclone & nor'easter tracks, coded in 
#       tropical_cyclone_filtering.R
#     - Fig. S7-19: Mann-Kendall significance testing results (NOT INCLUDED HERE,
#       CODED SEPARATELY IN hist_noreaster_revision1_additions.R)
#     - Fig. S20-21: frequency & seasonal frequency
#     - Fig. S23-26: intensity & intensification rate p-values by CRCM5-CMIP6
#       pilot model, 1950-2014 & 1979-2014
#     - Fig. S27-28: 1950-2014 seasonal intensity and intensification rate p-values
#     - Fig. S28-30: ERA5 & CFSR intensity/intensification rate trends up to present
#       (NOT INCLUDED HERE, CODED SEPARATELY IN hist_noreaster_revision1_additions.R)

# ==============================================================================

# Load hist_noreaster_intensity_analysis.R, which loads all data and conducts
#     all analyses needed for the following figures
source("/Path/to/analysis/file/hist_noreaster_intensity_analysis.R")

# ==============================================================================

# FIGURE 1

# Intensity survival function w/ each dataset & percentile comparison, 1979-2014

pdf("/Path/to/output/figures/fig_1_intensity_survival_functions.pdf", 
    family = "Helvetica-Narrow", width = 18, height = 6)
par(mfrow = c(1, 3), mar = c(1, 6, 1, 2), oma = c(4, 4, 4, 0.5))

# min. SLP
plot(crcm5_cmip6_min_slp_quants_1979_2014, my.probs*100, col = "#882255", 
     type = "l", main = "", lwd = 3, lty = 1.025, xlim = rev(c(90000, 98000)), 
     ylim = c(0, 100), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(crcm5_era5_min_slp_quants_1979_2014, my.probs*100, col = "#CC79A7", lwd = 3, lty = 1)
lines(era5_min_slp_quants_1979_2014, my.probs*100, col = "#00A9A5", lwd = 3, lty = 1)
lines(cfsr_min_slp_quants_1979_2014, my.probs*100, col = "#E69F00", lwd = 3, lty = 1)
mtext("Exceedance Probability", 2, 3.25, family = 'Helvetica', font = 2, cex = 2)
mtext("Min. SLP (hPa)", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 90000, y = 99, labels = "a", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(89500, 101000, 500), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2, las = 1, 
     at = seq(from = 90000, to = 100000, by = 2000), labels = c("900", "920", "940", "960", "980", "1000"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 0, to = 100, by = 20), labels = c("0", "20", "40", "60", "80", "100"))
box()

# max. wind
plot(rev(crcm5_cmip6_max_wind_quants_1979_2014), my.probs*100, col = "#882255", 
     type = "l", main = "", lwd = 3, lty = 1.025, xlim = c(18, 50), ylim = c(0, 100), 
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(crcm5_era5_max_wind_quants_1979_2014), my.probs*100, col = "#CC79A7", lwd = 3, lty = 1)
lines(rev(era5_max_wind_quants_1979_2014), my.probs*100, col = "#00A9A5", lwd = 3, lty = 1)
lines(rev(cfsr_max_wind_quants_1979_2014), my.probs*100, col = "#E69F00", lwd = 3, lty = 1)
legend("top", inset = c(0, -0.125), xpd = "True", 
       legend = c("CFSR", "ERA5", "CRCM5-ERA5", "CRCM5-CMIP6"), 
       fill = c("#E69F00", "#00A9A5", "#CC79A7", "#882255"), cex = 2.25, ncol = 4, bty = "n")
mtext("Max. Surface Wind Speed (m/s)", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Exceedance Probability", 2, 3.25, family = 'Helvetica', font = 2, cex = 2)
text(x = 50, y = 99, labels = "b", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(0, 60, 2), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2, las = 1, 
     at = seq(from = 0, to = 50, by = 10), labels = c("0", "10", "20", "30", "40", "50"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 0, to = 100, by = 20), labels = c("0", "20", "40", "60", "80", "100"))
box()

# Q-Q plot
plot(rev(crcm5_cmip6_min_slp_quants_1979_2014), 
     crcm5_cmip6_max_wind_quants_1979_2014, col = "#882255", pch = 19, cex = 2,
     main = "", xlim = rev(c(90000, 97500)), ylim = c(18, 50), xaxt = "n", 
     yaxt = "n", xlab = "", ylab = "")
points(rev(crcm5_cmip6_min_slp_quants_1979_2014)[101], 
       crcm5_cmip6_max_wind_quants_1979_2014[101],
       col = "black", bg = "#882255", pch = 21, cex = 2.5, lwd = 2)
points(rev(era5_min_slp_quants_1979_2014), era5_max_wind_quants_1979_2014, 
       col = "#00A9A5", pch = 15, cex = 2)
points(rev(era5_min_slp_quants_1979_2014)[101], era5_max_wind_quants_1979_2014[101],
       col = "black", bg = "#00A9A5", pch = 22, cex = 2.5, lwd = 2)
points(rev(crcm5_era5_min_slp_quants_1979_2014), crcm5_era5_max_wind_quants_1979_2014, 
       col = "#CC79A7", pch = 17, cex = 2)
points(rev(crcm5_era5_min_slp_quants_1979_2014)[101], 
       crcm5_era5_max_wind_quants_1979_2014[101],
       col = "black", bg = "#CC79A7", pch = 24, cex = 2, lwd = 2)
points(rev(cfsr_min_slp_quants_1979_2014), cfsr_max_wind_quants_1979_2014, 
       col = "#E69F00", pch = 18, cex = 2.5)
points(rev(cfsr_min_slp_quants_1979_2014)[101], cfsr_max_wind_quants_1979_2014[101],
       col = "black", bg = "#E69F00", pch = 23, cex = 2.5, lwd = 2)
mtext("Min. SLP Percentile Values", 1, 3.5, family = 'Helvetica', font = 2, cex = 2)
mtext("Max. Wind Percentile Values", 2, 3.25, family = 'Helvetica', font = 2, cex = 2)
text(x = 90000, y = 49.75, labels = "c", cex = 2.5, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, at = seq(89500, 101000, 500), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 1, 0), cex.axis = 2, las = 1, 
     at = seq(from = 90000, to = 100000, by = 2000), 
     labels = c("900", "920", "940", "960", "980", "1000"))
axis(side = 2, at = seq(0, 60, 2), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 10, to = 60, by = 10), labels = c("10", "20", "30", "40", "50", "60"))
box()

dev.off()

# ==============================================================================

# FIGURE 3

# Percentile vs. p-values intensity fig.
#   - min. SLP in top row, max. wind in bottom row
#   - 1950-2014 on left (ERA5 & CRCM5-CMIP6), 1979-2014 on right (all datasets)
#   - one line for each dataset on each plot

pdf("/Path/to/output/figures/fig_3_intensity_p_values.pdf", 
    family = "Helvetica-Narrow", width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), oma = c(3, 6, 6, 1))

# 1950-2014 min. SLP
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_1950_2014$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_1950_2014$p_values, lwd = 2, col = "#00A9A5")
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
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_1979_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_1979_2014$p_values, lwd = 2, col = "#00A9A5")
lines(rev(my.probs), crcm5_era5_slp_quant_trends_1979_2014$p_values, lwd = 2, col = "#CC79A7")
lines(rev(my.probs), cfsr_slp_quant_trends_1979_2014$p_values, lwd = 2, col = "#E69F00")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_1950_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_1950_2014$p_values, lwd = 2, col = "#00A9A5")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_1979_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_1979_2014$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_quant_trends_1979_2014$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_quant_trends_1979_2014$p_values, lwd = 2, col = "#E69F00")
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

# ==============================================================================

# FIGURE 4

# 1979-2014 seasonal intensity quantiles vs. p-values

pdf("/Path/to/output/figures/fig_4_seasonal_intensity_p_values_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(3, 3, 1, 1), oma = c(3, 7, 7, 1))

# min. SLP SON
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_1979_2014_SON$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#00A9A5")
lines(rev(my.probs), crcm5_era5_slp_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#CC79A7")
lines(rev(my.probs), cfsr_slp_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#E69F00")
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
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_1979_2014_DJF$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
lines(rev(my.probs), crcm5_era5_slp_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#CC79A7")
lines(rev(my.probs), cfsr_slp_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#E69F00")
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
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_1979_2014_MA$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#00A9A5")
lines(rev(my.probs), crcm5_era5_slp_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#CC79A7")
lines(rev(my.probs), cfsr_slp_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#E69F00")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_1979_2014_SON$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#E69F00")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_1979_2014_DJF$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#E69F00")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_1979_2014_MA$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#E69F00")
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

# FIGURE 5

# Percentile vs. p-values intensification rate fig.
#   - min. SLP in top row, max. wind in bottom row
#   - 1950-2014 on left (ERA5 & CRCM5-CMIP6), 1979-2014 on right (all datasets)
#   - one line for each dataset on each plot

pdf("/Path/to/output/figures/fig_5_intensification_p_values.pdf", 
    family = "Helvetica-Narrow", width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), oma = c(3, 6, 6, 1))

# 1950-2014 min. SLP
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_1950_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_1950_2014$p_values, lwd = 2, col = "#00A9A5")
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
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_1979_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_1979_2014$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_slp_intens_quant_trends_1979_2014$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_slp_intens_quant_trends_1979_2014$p_values, lwd = 2, col = "#E69F00")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_1950_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_1950_2014$p_values, lwd = 2, col = "#00A9A5")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_1979_2014$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_1979_2014$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_intens_quant_trends_1979_2014$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_intens_quant_trends_1979_2014$p_values, lwd = 2, col = "#E69F00")
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

# ==============================================================================

# FIGURE 6

# 1979-2014 seasonal intensification rate percentiles vs. p-values

pdf("/Path/to/output/figures/fig_6_seasonal_intensification_p_values_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(3, 3, 1, 1), oma = c(3, 7, 7, 1))

# min. SLP SON
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_1979_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_slp_intens_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_slp_intens_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#E69F00")
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
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_1979_2014_DJF$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_slp_intens_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_slp_intens_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#E69F00")
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
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_1979_2014_MA$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_slp_intens_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_slp_intens_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#E69F00")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_1979_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_intens_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_intens_quant_trends_1979_2014_SON$p_values, lwd = 2, col = "#E69F00")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_1979_2014_DJF$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_intens_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_intens_quant_trends_1979_2014_DJF$p_values, lwd = 2, col = "#E69F00")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_1979_2014_MA$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#00A9A5")
lines(my.probs, crcm5_era5_wind_intens_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#CC79A7")
lines(my.probs, cfsr_wind_intens_quant_trends_1979_2014_MA$p_values, lwd = 2, col = "#E69F00")
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

# ==============================================================================

# FIGURE 7

# Track density figure using grid counts for tracks instead of kde2d, separate
#     colorbars for each dataset -- 1979-2014
#     **1 row, 4 columns: columns = CFSR, ERA5, CRCM5-ERA5, CRCM5-CMIP6

# Define color scales for each dataset

# CFSR
cfsr_cols <- colorRampPalette(c(
  "#FFFFFF",   # white
  "#F8EDD7",
  "#F1D8AB",
  "#E8C078",
  "#DDA338",   # original
  "#B1822D",
  "#856222",
  "#584116"    # darkest
))

# ERA5
era5_cols <- colorRampPalette(c(
  "#FFFFFF",   # white
  "#DBEDED",
  "#B3DAD9",
  "#85C2C1",
  "#4BA6A4",   # original
  "#3C8583",
  "#2D6462",
  "#1E4242"    # darkest
))

# CRCM5-ERA5
crcm5_era5_cols <- colorRampPalette(c(
  "#FFFFFF",   # white
  "#F2E5ED",
  "#E5C8D9",
  "#D4A7C2",
  "#C07DA5",   # original
  "#AD7195",
  "#9A6484",
  "#865874"    # darkest
))

# CRCM5-CMIP6
crcm5_cmip6_cols <- colorRampPalette(c(
  "#FFFFFF",   # white
  "#EBDDE4",
  "#D8BCCA",
  "#C093AA",
  "#A66687",
  "#904069",
  "#7C1E4E",   # original
  "#4A122F"    # darkest
))

pdf("/Path/to/output/figures/fig_7_track_locations_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(3, 4, 3, 2), oma = c(1, 1, 1, 1))

# CFSR
image.plot(cfsr_track_counts_1979_2014, col = cfsr_cols(300), 
           breaks = colbar_breaks_cfsr_track_counts_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0), 
           at = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), 
           labels = c("0", "0.1", "0.3", "0.3", "0.4", "0.5"), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("CFSR", 3, 1.25, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "a", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# ERA5
image.plot(era5_track_counts_1979_2014, col = era5_cols(300), 
           breaks = colbar_breaks_era5_track_counts_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0), 
           at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), 
           labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2"), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("ERA5", 3, 1.25, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "b", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# CRCM5-ERA5
image.plot(crcm5_track_counts_1979_2014$era5, col = crcm5_era5_cols(300), 
           breaks = colbar_breaks_crcm5_era5_track_counts_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0), 
           at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), 
           labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2"), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("CRCM5-ERA5", 3, 1.25, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "c", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# CRCM5-CMIP6 ensemble
image.plot(crcm5_cmip6_ensemble_counts_1979_2014, col = crcm5_cmip6_cols(300), 
           breaks = colbar_breaks_crcm5_cmip6_track_counts_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0), 
           at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2), 
           labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0", "1.2"), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("CRCM5-CMIP6", 3, 1.25, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "d", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

dev.off()

# ==============================================================================

# FIGURE 8

# Max. intensity density figure using kde2d -- 1979-2014 DATASET
#     **2 row, 4 columns: columns = CFSR, ERA5, CRCM5-ERA5, CRCM5-CMIP6,
#     rows = min. SLP, max. wind

pdf("/Path/to/output/figures/fig_8_max_intensity_locations_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 9)
par(mfrow = c(2, 4), mar = c(3, 5, 1, 3), oma = c(3, 6, 6, 6))

# -- -- -- -- -- -- -- -- -- -- --

# Min. SLP

# CFSR
image.plot(cfsr_min_slp_kde_1979_2014, col = cfsr_cols(300), 
           breaks = colbar_breaks_cfsr_min_slp_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.0005, 0.001, 0.0015, 0.002), 
           labels = c("0", "", "0.001", "", "0.002"), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("Min. SLP", 2, 4, family = 'Helvetica', font = 2, cex = 2)
mtext("CFSR", 3, 3, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "a", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# ERA5
image.plot(era5_min_slp_kde_1979_2014, col = era5_cols(300), 
           breaks = colbar_breaks_era5_min_slp_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.001, 0.002, 0.003, 0.004, 0.005), 
           labels = c("0", "", "0.002", "", "0.004", ""), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("ERA5", 3, 3, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "b", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# CRCM5-ERA5
image.plot(crcm5_era5_min_slp_kde_1979_2014$era5, col = crcm5_era5_cols(300), 
           breaks = colbar_breaks_crcm5_era5_min_slp_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.0005, 0.001, 0.0015), 
           labels = c("0", "0.0005", "0.001", "0.0015"), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("CRCM5-ERA5", 3, 3, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "c", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# CRCM5-CMIP6 ensemble
image.plot(crcm5_cmip6_ensemble_min_slp_kde_1979_2014, col = crcm5_cmip6_cols(300), 
           breaks = colbar_breaks_crcm5_cmip6_min_slp_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.0005, 0.001, 0.0015, 0.002), 
           labels = c("0", "", "0.001", "", "0.002"), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("CRCM5-CMIP6", 3, 3, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "d", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# -- -- -- -- -- -- -- -- -- -- --

# Max. wind

# CFSR
image.plot(cfsr_max_wind_kde_1979_2014, col = cfsr_cols(300), 
           breaks = colbar_breaks_cfsr_max_wind_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.0005, 0.001, 0.0015), 
           labels = c("0", "0.0005", "0.001", "0.0015"), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("Max. wind", 2, 4, family = 'Helvetica', font = 2, cex = 2)
text(x = -13, y = 77, labels = "e", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# ERA5
image.plot(era5_max_wind_kde_1979_2014, col = era5_cols(300), 
           breaks = colbar_breaks_era5_max_wind_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.001, 0.002, 0.003, 0.004), 
           labels = c("0", "", "0.002", "", "0.004"), box()))
map(database = "world", add = T)
map("state", add = T)
text(x = -13, y = 77, labels = "f", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# CRCM5-ERA5
image.plot(crcm5_era5_max_wind_kde_1979_2014$era5, col = crcm5_era5_cols(300), 
           breaks = colbar_breaks_crcm5_era5_max_wind_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.0005, 0.001, 0.0015), 
           labels = c("0", "0.0005", "0.001", "0.0015"), box()))
map(database = "world", add = T)
map("state", add = T)
text(x = -13, y = 77, labels = "g", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# CRCM5-CMIP6 ensemble
image.plot(crcm5_cmip6_ensemble_max_wind_kde_1979_2014, col = crcm5_cmip6_cols(300), 
           breaks = colbar_breaks_crcm5_cmip6_max_wind_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
            at = c(0, 0.0005, 0.001, 0.0015, 0.002), 
            labels = c("0", "", "0.001", "", "0.002"), box()))
map(database = "world", add = T)
map("state", add = T)
text(x = -13, y = 77, labels = "h", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

dev.off()

# ==============================================================================

# FIGURE 9

# Max. intensification rate density figure using kde2d -- 1979-2014 DATASET
#     **2 row, 4 columns: columns = CFSR, ERA5, CRCM5-ERA5, CRCM5-CMIP6,
#     rows = min. SLP, max. wind

pdf("/Path/to/output/figures/fig_9_max_intensification_locations_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 9)
par(mfrow = c(2, 4), mar = c(3, 5, 1, 3), oma = c(3, 6, 6, 6))

# -- -- -- -- -- -- -- -- -- -- --

# Min. SLP

# CFSR
image.plot(cfsr_slp_intens_kde_1979_2014, col = cfsr_cols(300), 
           breaks = colbar_breaks_cfsr_slp_intens_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.0005, 0.001, 0.0015, 0.002), 
           labels = c("0", "", "0.001", "", "0.002"), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("Min. SLP rate", 2, 4, family = 'Helvetica', font = 2, cex = 2)
mtext("CFSR", 3, 3, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "a", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# ERA5
image.plot(era5_slp_intens_kde_1979_2014, col = era5_cols(300), 
           breaks = colbar_breaks_era5_slp_intens_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.001, 0.002, 0.003, 0.004, 0.005), 
           labels = c("0", "", "0.002", "", "0.004", ""), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("ERA5", 3, 3, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "b", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# CRCM5-ERA5
image.plot(crcm5_era5_slp_intens_kde_1979_2014$era5, col = crcm5_era5_cols(300), 
           breaks = colbar_breaks_crcm5_era5_slp_intens_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.0005, 0.001, 0.0015, 0.002, 0.0025), 
           labels = c("0", "", "0.001", "", "0.002", ""), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("CRCM5-ERA5", 3, 3, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "c", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# CRCM5-CMIP6 ensemble
image.plot(crcm5_cmip6_ensemble_slp_intens_kde_1979_2014, col = crcm5_cmip6_cols(300), 
           breaks = colbar_breaks_crcm5_cmip6_slp_intens_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003), 
           labels = c("0", "", "0.001", "", "0.002", "", "0.003"), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("CRCM5-CMIP6", 3, 3, family = 'Helvetica', font = 2, cex = 1.75)
text(x = -13, y = 77, labels = "d", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# -- -- -- -- -- -- -- -- -- -- --

# Max. wind

# CFSR
image.plot(cfsr_wind_intens_kde_1979_2014, col = cfsr_cols(300), 
           breaks = colbar_breaks_cfsr_wind_intens_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.0005, 0.001, 0.0015, 0.002), 
           labels = c("0", "", "0.001", "", "0.002"), box()))
map(database = "world", add = T)
map("state", add = T)
mtext("Max. wind rate", 2, 4, family = 'Helvetica', font = 2, cex = 2)
text(x = -13, y = 77, labels = "e", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# ERA5
image.plot(era5_wind_intens_kde_1979_2014, col = era5_cols(300), 
           breaks = colbar_breaks_era5_wind_intens_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.001, 0.002, 0.003, 0.004), 
           labels = c("0", "", "0.002", "", "0.004"), box()))
map(database = "world", add = T)
map("state", add = T)
text(x = -13, y = 77, labels = "f", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# CRCM5-ERA5
image.plot(crcm5_era5_wind_intens_kde_1979_2014$era5, col = crcm5_era5_cols(300), 
           breaks = colbar_breaks_crcm5_era5_wind_intens_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.0005, 0.001, 0.0015, 0.002, 0.0025), 
           labels = c("0", "", "0.001", "", "0.002", ""), box()))
map(database = "world", add = T)
map("state", add = T)
text(x = -13, y = 77, labels = "g", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

# CRCM5-CMIP6 ensemble
image.plot(crcm5_cmip6_ensemble_wind_intens_kde_1979_2014, col = crcm5_cmip6_cols(300), 
           breaks = colbar_breaks_crcm5_cmip6_wind_intens_1979_2014, 
           xlim = c(-110, -7), ylim = c(15, 82), yaxt = "n", xaxt = "n", legend.width = 2.5, 
           axis.args = list(tck = -0.4, cex.axis = 2, las = 0, mgp = c(1, 1.2, 0),
           at = c(0, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003), 
           labels = c("0", "", "0.001", "", "0.002", "", "0.003"), box()))
map(database = "world", add = T)
map("state", add = T)
text(x = -13, y = 77, labels = "h", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
axis(side = 1, tck = 0.015, mgp = c(1, 0.8, 0), cex.axis = 2, las = 1, 
     at = seq(from = -100, to = -20, by = 20), labels = c("-100", "-80", "-60", "-40", "-20"))
axis(side = 2, tck = 0.015, mgp=c(1, 0.4, 0), las = 1, cex.axis = 2, 
     at = seq(from = 20, to = 80, by = 20), labels = c("20", "40", "60", "80"))
axis(side = 3, tck = 0.015, labels = FALSE, at = seq(from = -100, to = -40, by = 10))
axis(side = 4, tck = 0.015, labels = FALSE, at = seq(from = 20, to = 80, by = 10))
box()

dev.off()

# ==============================================================================

# FIGURE S1

# Spatial domain map

# Extract CRCM5 domain boundary (perimeter of the rotated grid)
nrow_c <- nrow(crcm5_lat)
ncol_c <- ncol(crcm5_lat)

top    <- cbind(1,            1:ncol_c)
right  <- cbind(2:nrow_c,    ncol_c)
bottom <- cbind(nrow_c,      (ncol_c-1):1)
left   <- cbind((nrow_c-1):2, 1)
perim  <- rbind(top, right, bottom, left)

crcm5_perim_lon <- crcm5_lon[perim]
crcm5_perim_lat <- crcm5_lat[perim]

crcm5_perim_lon_plot <- crcm5_perim_lon
crcm5_perim_lon_plot[crcm5_perim_lon > 0] <- crcm5_perim_lon[crcm5_perim_lon > 0] - 360

keep <- crcm5_perim_lon_plot > -200 & crcm5_perim_lon_plot < 0
crcm5_perim_lon_plot <- crcm5_perim_lon_plot[keep]
crcm5_perim_lat_plot <- crcm5_perim_lat[keep]

pdf("/Path/to/output/figures/fig_s1_spatial_domain_map.pdf",
    family = "Helvetica-Narrow", width = 10, height = 7)

par(xaxs = "i", yaxs = "i", mar = c(1, 1, 1, 1), oma = c(3, 3, 3, 1))
plot.new()
plot.window(xlim = c(-186, -5), ylim = c(5, 85))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
     col = "#F4F9FD", border = NA)
map("world", xlim = c(-186, -5), ylim = c(5, 85),
    col = "gray70", fill = TRUE, bg = "white",
    border = "gray90", lwd = 0.5, add = TRUE)
polygon(crcm5_perim_lon_plot, crcm5_perim_lat_plot,
        border = "#882255", col = adjustcolor("#CC79A7", alpha.f = 0.1),
        lwd = 3)
rect(min(cfsr_grid$lon_seq), min(cfsr_grid$lat_seq),
     max(cfsr_grid$lon_seq), max(cfsr_grid$lat_seq),
     border = "#E69F00", col = adjustcolor("#E69F00", alpha.f = 0.1),
     lwd = 3)
rect(min(era5_grid$lon_seq), min(era5_grid$lat_seq),
     max(era5_grid$lon_seq), max(era5_grid$lat_seq),
     border = "#00A9A5", col = adjustcolor("#00A9A5", alpha.f = 0.1),
     lwd = 3)
polygon(c(-80, -70, -70, -80), c(33, 33, 45, 45),
        border = "black", lwd = 3,
        col = "gray30", density = 10, angle = 45)
axis(side = 1, at = seq(-190, 0, 10), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.5, las = 1,
     at = seq(-180, -10, by = 20),
     labels = paste0(seq(180, 10, by = -20), "°W"))
axis(side = 2, at = seq(0, 90, 5), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp = c(1, 0.4, 0), cex.axis = 1.5, las = 1,
     at = seq(10, 80, by = 10),
     labels = paste0(seq(10, 80, by = 10), "°N"))
legend("top", inset = c(0, -0.1), xpd = NA,
       legend = c("CFSR", "ERA5", "CRCM5", "East Coast"),
       col    = c("#E69F00", "#00A9A5", "#882255", "black"),
       lwd    = 2, bty = "n", cex = 1.25,
       fill   = adjustcolor(c("#E69F00", "#00A9A5", "#882255", "black"), alpha.f = 0.1),
       border = c("#E69F00", "#00A9A5", "#882255", "black"), ncol = 4,
       density = c(NA, NA, NA, 8), angle  = c(NA, NA, NA, 45))
box()
dev.off()

# ==============================================================================

# FIGURE S2

# MPI member survival functions, intensity p-values, & intensification rate p-values, 1979-2014

pdf("/Path/to/output/figures/fig_s2_crcm5_mpi_members_1979_2014.pdf", 
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
plot(rev(my.probs), mpi_slp_quant_trends_1979_2014$p_values$sim1, type = "l", 
     col = ("#1B3A6B"), lwd = 3, xlim = c(0, 1), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), mpi_slp_quant_trends_1979_2014$p_values$sim2, lwd = 3, col = "#2E6DB4")
lines(rev(my.probs), mpi_slp_quant_trends_1979_2014$p_values$sim3, lwd = 3, col = "#4D9FD6")
lines(rev(my.probs), mpi_slp_quant_trends_1979_2014$p_values$sim4, lwd = 3, col = "#89C4E1")
lines(rev(my.probs), mpi_slp_quant_trends_1979_2014$p_values$sim5, lwd = 3, col = "#C2DFF0")
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
plot(my.probs, mpi_wind_quant_trends_1979_2014$p_values$sim1, type = "l", 
     col = ("#1B3A6B"), lwd = 3, xlim = c(0, 1), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, mpi_wind_quant_trends_1979_2014$p_values$sim2, lwd = 3, col = "#2E6DB4")
lines(my.probs, mpi_wind_quant_trends_1979_2014$p_values$sim3, lwd = 3, col = "#4D9FD6")
lines(my.probs, mpi_wind_quant_trends_1979_2014$p_values$sim4, lwd = 3, col = "#89C4E1")
lines(my.probs, mpi_wind_quant_trends_1979_2014$p_values$sim5, lwd = 3, col = "#C2DFF0")
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
plot(my.probs, mpi_slp_intens_quant_trends_1979_2014$p_values$sim1, type = "l", 
     col = ("#1B3A6B"), lwd = 3, xlim = c(0, 1), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, mpi_slp_intens_quant_trends_1979_2014$p_values$sim2, lwd = 3, col = "#2E6DB4")
lines(my.probs, mpi_slp_intens_quant_trends_1979_2014$p_values$sim3, lwd = 3, col = "#4D9FD6")
lines(my.probs, mpi_slp_intens_quant_trends_1979_2014$p_values$sim4, lwd = 3, col = "#89C4E1")
lines(my.probs, mpi_slp_intens_quant_trends_1979_2014$p_values$sim5, lwd = 3, col = "#C2DFF0")
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
plot(my.probs, mpi_wind_intens_quant_trends_1979_2014$p_values$sim1, type = "l", 
     col = ("#1B3A6B"), lwd = 3, xlim = c(0, 1), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, mpi_wind_intens_quant_trends_1979_2014$p_values$sim2, lwd = 3, col = "#2E6DB4")
lines(my.probs, mpi_wind_intens_quant_trends_1979_2014$p_values$sim3, lwd = 3, col = "#4D9FD6")
lines(my.probs, mpi_wind_intens_quant_trends_1979_2014$p_values$sim4, lwd = 3, col = "#89C4E1")
lines(my.probs, mpi_wind_intens_quant_trends_1979_2014$p_values$sim5, lwd = 3, col = "#C2DFF0")
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

# ==============================================================================

# FIGURE S20

# Historical frequency fig.

pdf("/Path/to/output/figures/fig_s20_frequency_time_series.pdf", 
    family = "Helvetica-Narrow", width = 8, height = 12)
par(mfrow = c(3, 1), mar = c(2.5, 3, 1, 1), oma = c(5, 5, 5, 1))

# All storms
plot(crcm5_noreaster_mme_df$year, crcm5_noreaster_mme_df$mean_count, type = "l", lwd = 3, 
     col = "#882255", xlab = "Year", ylab = "Annual Storm Frequency", 
     xlim = c(1940, 2025), ylim = c(10, 50), yaxt = "n", xaxt = "n")
# Add shading for ±1 std dev
polygon(c(crcm5_noreaster_mme_df$year[1:64], rev(crcm5_noreaster_mme_df$year[1:64])),
        c(crcm5_noreaster_mme_df$mean_count[1:64] - crcm5_noreaster_mme_df$std_count[1:64], 
          rev(crcm5_noreaster_mme_df$mean_count[1:64] + crcm5_noreaster_mme_df$std_count[1:64])),
        col = adjustcolor("#88225560", alpha.f = 0.3), border = NA)

lines(crcm5_noreaster_frequencies$era5, col = "#CC79A7", lwd = 3)
lines(era5_noreaster_freq$year, era5_noreaster_freq$count, col = "#00A9A5", lwd = 3)
lines(cfsr_noreaster_freq$year, cfsr_noreaster_freq$count, col = "#E69F00", lwd = 3)

legend("topright", inset = c(0.35, -0.2), xpd = "True", 
       legend = c("CRCM5-CMIP6", "CRCM5-CMIP6 ±1 std dev"), 
       fill = c("#882255", "#88225560"), cex = 1.5, ncol = 1, bty = "n")
legend("topright", inset = c(0.125, -0.2), xpd = "True", 
       legend = c("CRCM5-ERA5", "ERA5"), 
       fill = c("#CC79A7", "#00A9A5"), cex = 1.5, ncol = 1, bty = "n")
legend("topright", inset = c(0, -0.2), xpd = "True", legend = c("CFSR"), 
       fill = c("#E69F00"), cex = 1.5, ncol = 1, bty = "n")

text(x = 2025, y = 48.5, labels = "a", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")

axis(side = 1, at = seq(1935, 2030, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.5, las = 1, 
     at = seq(from = 1940, to = 2020, by = 10), 
     labels = c("1940", "1950", "1960", "1970", "1979", "1990", "2000", "2010", "2020"))
axis(side = 2, at = seq(0, 55, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp=c(1, 0.2, 0), las = 1, cex.axis = 1.5, 
     at = seq(from = 15, to = 55, by = 10), labels = c("15", "25", "35", "45", "55"))
axis(side = 4, tck = 0.007, at = seq(from = 0, to = 40, by = 1), labels = FALSE, col.ticks = "#707070")

box(lwd = 1.5)
#
#
# Weak storms
plot(weak_crcm5_noreaster_mme_df$year, weak_crcm5_noreaster_mme_df$mean_count, type = "l", 
     lwd = 3, col = "#882255", xlab = "Year", ylab = "Annual Storm Frequency", 
     xlim = c(1940, 2025), ylim = c(0, 35), yaxt = "n", xaxt = "n")
# Add shading for ±1 std dev
polygon(c(weak_crcm5_noreaster_mme_df$year[1:64], rev(weak_crcm5_noreaster_mme_df$year[1:64])),
        c(weak_crcm5_noreaster_mme_df$mean_count[1:64] - weak_crcm5_noreaster_mme_df$std_count[1:64], 
          rev(weak_crcm5_noreaster_mme_df$mean_count[1:64] + weak_crcm5_noreaster_mme_df$std_count[1:64])),
        col = adjustcolor("#88225560", alpha.f = 0.3), border = NA)

lines(weak_crcm5_noreaster_frequencies$era5, col = "#CC79A7", lwd = 3)
lines(era5_weak_noreaster_freq$year, era5_weak_noreaster_freq$count, col = "#00A9A5", lwd = 3)
lines(cfsr_weak_noreaster_freq$year, cfsr_weak_noreaster_freq$count, col = "#E69F00", lwd = 3)

mtext("Storm Occurrences (# per storm season)", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 2025, y = 33.5, labels = "b", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")

axis(side = 1, at = seq(1935, 2030, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.5, las = 1, 
     at = seq(from = 1940, to = 2020, by = 10), 
     labels = c("1940", "1950", "1960", "1970", "1979", "1990", "2000", "2010", "2020"))
axis(side = 2, at = seq(0, 55, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp=c(1, 0.2, 0), las = 1, cex.axis = 1.5, 
     at = seq(from = 0, to = 35, by = 5), labels = c("0", "5", "10", "15", "20", "25", "30", "35"))
axis(side = 4, tck = 0.007, at = seq(from = 0, to = 40, by = 1), labels = FALSE, col.ticks = "#707070")

box(lwd = 1.5)
#
#
# Strong storms
plot(strong_crcm5_noreaster_mme_df$year, strong_crcm5_noreaster_mme_df$mean_count, type = "l", lwd = 3, 
     col = "#882255", xlab = "Year", ylab = "Annual Storm Frequency", 
     xlim = c(1940, 2025), ylim = c(0, 30), yaxt = "n", xaxt = "n")
# Add shading for ±1 std dev
polygon(c(strong_crcm5_noreaster_mme_df$year[1:64], rev(strong_crcm5_noreaster_mme_df$year[1:64])),
        c(strong_crcm5_noreaster_mme_df$mean_count[1:64] - strong_crcm5_noreaster_mme_df$std_count[1:64], 
          rev(strong_crcm5_noreaster_mme_df$mean_count[1:64] + strong_crcm5_noreaster_mme_df$std_count[1:64])),
        col = adjustcolor("#88225560", alpha.f = 0.3), border = NA)

lines(strong_crcm5_noreaster_frequencies$era5, col = "#CC79A7", lwd = 3)
lines(era5_strong_noreaster_freq$year, era5_strong_noreaster_freq$count, col = "#00A9A5", lwd = 3)
lines(cfsr_strong_noreaster_freq$year, cfsr_strong_noreaster_freq$count, col = "#E69F00", lwd = 3)

text(x = 2025, y = 29.25, labels = "c", cex = 3, font = 4, col = "black", family = "Helvetica-Narrow")
mtext("Year", 1, 4, family = 'Helvetica', font = 2, cex = 2)

axis(side = 1, at = seq(1935, 2030, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.5, las = 1, 
     at = seq(from = 1940, to = 2020, by = 10), 
     labels = c("1940", "1950", "1960", "1970", "1979", "1990", "2000", "2010", "2020"))
axis(side = 2, at = seq(0, 30, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp=c(1, 0.2, 0), las = 1, cex.axis = 1.5, 
     at = seq(from = 0, to = 30, by = 5), labels = c("0", "5", "10", "15", "20", "25", "30"))
axis(side = 4, tck = 0.007, at = seq(from = 0, to = 40, by = 1), labels = FALSE, col.ticks = "#707070")

box(lwd = 1.5)

dev.off()

# ==============================================================================

# FIGURE S21

# Strong storm seasonal frequency

pdf("/Path/to/output/figures/fig_s21_seasonal_frequency_time_series.pdf", 
    family = "Helvetica-Narrow", width = 8, height = 12)
par(mfrow = c(3, 1), mar = c(2.5, 3, 1, 1), oma = c(5, 5, 5, 1))

# SON
plot(strong_crcm5_noreaster_SON_mme_df$year, 
     strong_crcm5_noreaster_SON_mme_df$mean_count, type = "l", 
     lwd = 3, col = "#882255", xlab = "Year", ylab = "Annual Storm Frequency", 
     xlim = c(1940, 2025), ylim = c(0, 9), yaxt = "n", xaxt = "n")
# Add shading for ±1 std dev
polygon(c(strong_crcm5_noreaster_SON_mme_df$year[1:64], rev(strong_crcm5_noreaster_SON_mme_df$year[1:64])),
        c(strong_crcm5_noreaster_SON_mme_df$mean_count[1:64] - strong_crcm5_noreaster_SON_mme_df$std_count[1:64], 
          rev(strong_crcm5_noreaster_SON_mme_df$mean_count[1:64] + strong_crcm5_noreaster_SON_mme_df$std_count[1:64])),
        col = adjustcolor("#88225560", alpha.f = 0.3), border = NA)

lines(strong_crcm5_noreaster_SON_freq$era5, col = "#CC79A7", lwd = 3)
lines(era5_strong_noreaster_SON_freq$year, 
      era5_strong_noreaster_SON_freq$count, col = "#00A9A5", lwd = 3)
lines(cfsr_strong_noreaster_SON_freq$year, 
      cfsr_strong_noreaster_SON_freq$count, col = "#E69F00", lwd = 3)

legend("topright", inset = c(0.225, -0.2), xpd = "True", 
       legend = c("CRCM5-CMIP6", "CRCM5-CMIP6 ±1 std dev"), 
       fill = c("#882255", "#88225560"), cex = 1.5, ncol = 1, bty = "n")
legend("topright", inset = c(0, -0.25), xpd = "True", 
       legend = c("CFSR", "ERA5", "CRCM5-ERA5"), 
       fill = c("#E69F00", "#00A9A5", "#CC79A7"), cex = 1.5, ncol = 1, bty = "n")

text(x = 2023, y = 8.75, labels = "a) SON", cex = 2, font = 4, col = "black", family = "Helvetica-Narrow")

axis(side = 1, at = seq(1935, 2030, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.5, las = 1, 
     at = seq(from = 1940, to = 2020, by = 10), 
     labels = c("1940", "1950", "1960", "1970", "1979", "1990", "2000", "2010", "2020"))
axis(side = 2, at = seq(0, 55, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp=c(1, 0.2, 0), las = 1, cex.axis = 1.5, 
     at = seq(from = 0, to = 10, by = 2), labels = c("0", "2", "4", "6", "8", "10"))
axis(side = 4, tck = 0.007, at = seq(from = 0, to = 40, by = 1), labels = FALSE, col.ticks = "#707070")

box(lwd = 1.5)
#
#
# DJF
plot(strong_crcm5_noreaster_DJF_mme_df$year, 
     strong_crcm5_noreaster_DJF_mme_df$mean_count, type = "l", lwd = 3, 
     col = "#882255", xlab = "Year", ylab = "Annual Storm Frequency", 
     xlim = c(1940, 2025), ylim = c(0, 17), yaxt = "n", xaxt = "n")
# Add shading for ±1 std dev
polygon(c(strong_crcm5_noreaster_DJF_mme_df$year[1:64], rev(strong_crcm5_noreaster_DJF_mme_df$year[1:64])),
        c(strong_crcm5_noreaster_DJF_mme_df$mean_count[1:64] - strong_crcm5_noreaster_DJF_mme_df$std_count[1:64], 
          rev(strong_crcm5_noreaster_DJF_mme_df$mean_count[1:64] + strong_crcm5_noreaster_DJF_mme_df$std_count[1:64])),
        col = adjustcolor("#88225560", alpha.f = 0.3), border = NA)

lines(strong_crcm5_noreaster_DJF_freq$era5, col = "#CC79A7", lwd = 3)
lines(era5_strong_noreaster_DJF_freq$year, 
      era5_strong_noreaster_DJF_freq$count, col = "#00A9A5", lwd = 3)
lines(cfsr_strong_noreaster_DJF_freq$year, 
      cfsr_strong_noreaster_DJF_freq$count, col = "#E69F00", lwd = 3)

mtext("Storm Occurrences (# per storm season)", 2, 3.5, family = 'Helvetica', font = 2, cex = 2)
text(x = 2023, y = 16.5, labels = "b) DJF", cex = 2, font = 4, col = "black", family = "Helvetica-Narrow")

axis(side = 1, at = seq(1935, 2030, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.5, las = 1, 
     at = seq(from = 1940, to = 2020, by = 10), labels = c("1940", "1950", "1960", "1970", "1979", "1990", "2000", "2010", "2020"))
axis(side = 2, at = seq(0, 55, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp=c(1, 0.2, 0), las = 1, cex.axis = 1.5, 
     at = seq(from = 0, to = 35, by = 5), labels = c("0", "5", "10", "15", "20", "25", "30", "35"))
axis(side = 4, tck = 0.007, at = seq(from = 0, to = 40, by = 1), labels = FALSE, col.ticks = "#707070")

box(lwd = 1.5)
#
#
# MA
plot(strong_crcm5_noreaster_MA_mme_df$year, 
     strong_crcm5_noreaster_MA_mme_df$mean_count, type = "l", lwd = 3, 
     col = "#882255", xlab = "Year", ylab = "Annual Storm Frequency", 
     xlim = c(1940, 2025), ylim = c(0, 13), yaxt = "n", xaxt = "n")
# Add shading for ±1 std dev
polygon(c(strong_crcm5_noreaster_MA_mme_df$year[1:64], rev(strong_crcm5_noreaster_MA_mme_df$year[1:64])),
        c(strong_crcm5_noreaster_MA_mme_df$mean_count[1:64] - strong_crcm5_noreaster_MA_mme_df$std_count[1:64], 
          rev(strong_crcm5_noreaster_MA_mme_df$mean_count[1:64] + strong_crcm5_noreaster_MA_mme_df$std_count[1:64])),
        col = adjustcolor("#88225560", alpha.f = 0.3), border = NA)

lines(strong_crcm5_noreaster_MA_freq$era5, col = "#CC79A7", lwd = 3)
lines(era5_strong_noreaster_MA_freq$year, 
      era5_strong_noreaster_MA_freq$count, col = "#00A9A5", lwd = 3)
lines(cfsr_strong_noreaster_MA_freq$year, 
      cfsr_strong_noreaster_MA_freq$count, col = "#E69F00", lwd = 3)

text(x = 2023, y = 12.5, labels = "c) MA", cex = 2, font = 4, col = "black", family = "Helvetica-Narrow")
mtext("Year", 1, 4, family = 'Helvetica', font = 2, cex = 2)

axis(side = 1, at = seq(1935, 2030, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 1, tck = 0.015, mgp = c(1, 0.6, 0), cex.axis = 1.5, las = 1, 
     at = seq(from = 1940, to = 2020, by = 10), labels = c("1940", "1950", "1960", "1970", "1979", "1990", "2000", "2010", "2020"))
axis(side = 2, at = seq(0, 30, 1), labels = FALSE, tck = 0.01, col.ticks = "#707070")
axis(side = 2, tck = 0.015, mgp=c(1, 0.2, 0), las = 1, cex.axis = 1.5, 
     at = seq(from = 0, to = 12, by = 3), labels = c("0", "3", "6", "9", "12"))
axis(side = 4, tck = 0.007, at = seq(from = 0, to = 40, by = 1), labels = FALSE, col.ticks = "#707070")

box(lwd = 1.5)

dev.off()

# ==============================================================================

# FIGURE S23

# Plot intensity quantiles vs. p-values by model: 1950-2014

pdf("/Path/to/output/figures/fig_s23_crcm5_cmip6_pilots_intensity_p_values_1950_2014.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 5, 2))

# CanESM5 min SLP
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_1950_2014$p_values$canesm5, 
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
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_1950_2014$p_values$cnrm, 
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
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_1950_2014$p_values$mpi, 
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
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_1950_2014$p_values$noresm2, 
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
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_1950_2014$p_values$canesm5, 
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
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_1950_2014$p_values$cnrm, 
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
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_1950_2014$p_values$mpi, 
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
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_1950_2014$p_values$noresm2, 
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

# ==============================================================================

# FIGURE S24

# Plot intensity quantiles vs. p-values by model: 1979-2014

pdf("/Path/to/output/figures/fig_s24_crcm5_cmip6_pilots_intensity_p_values_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 5, 2))

# CanESM5 min SLP
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_1979_2014$p_values$canesm5, 
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
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_1979_2014$p_values$cnrm, 
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
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_1979_2014$p_values$mpi, 
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
plot(rev(my.probs), crcm5_cmip6_models_slp_quant_trends_1979_2014$p_values$noresm2, 
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
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_1979_2014$p_values$canesm5, 
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
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_1979_2014$p_values$cnrm, 
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
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_1979_2014$p_values$mpi, 
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
plot(my.probs, crcm5_cmip6_models_wind_quant_trends_1979_2014$p_values$noresm2, 
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

# ==============================================================================

# FIGURE S25

# Plot intensification rate quantiles vs. p-values by model: 1950-2014

pdf("/Path/to/output/figures/fig_s25_crcm5_cmip6_pilots_intensification_p_values_1950_2014.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 5, 2))

# CanESM5 min SLP
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_1950_2014$p_values$canesm5, 
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
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_1950_2014$p_values$cnrm, 
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
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_1950_2014$p_values$mpi, 
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
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_1950_2014$p_values$noresm2, 
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
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_1950_2014$p_values$canesm5, 
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
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_1950_2014$p_values$cnrm, 
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
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_1950_2014$p_values$mpi, 
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
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_1950_2014$p_values$noresm2, 
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

# ==============================================================================

# FIGURE S26

# Plot intensification rate quantiles vs. p-values by model: 1979-2014

pdf("/Path/to/output/figures/fig_s26_crcm5_cmip6_pilots_intensification_p_values_1979_2014.pdf", 
    family = "Helvetica-Narrow", width = 20, height = 10)
par(mfrow = c(2, 4), mar = c(3, 6, 1, 1), oma = c(3, 6, 5, 2))

# CanESM5 min SLP
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_1979_2014$p_values$canesm5, 
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
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_1979_2014$p_values$cnrm, 
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
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_1979_2014$p_values$mpi, 
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
plot(my.probs, crcm5_cmip6_models_slp_intens_quant_trends_1979_2014$p_values$noresm2, 
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
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_1979_2014$p_values$canesm5, 
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
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_1979_2014$p_values$cnrm, 
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
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_1979_2014$p_values$mpi, 
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
plot(my.probs, crcm5_cmip6_models_wind_intens_quant_trends_1979_2014$p_values$noresm2, 
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

# ==============================================================================

# FIGURE S27

# 1950-2014 seasonal intensity quantiles vs. p-values

pdf("/Path/to/output/figures/fig_s27_seasonal_intensity_p_values_1950_2014.pdf", 
    family = "Helvetica-Narrow", width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(3, 3, 1, 1), oma = c(3, 7, 7, 1))

# min. SLP SON
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_1950_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_1950_2014_SON$p_values, lwd = 2, col = "#00A9A5")
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
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_1950_2014_DJF$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_1950_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
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
plot(rev(my.probs), crcm5_cmip6_ensemble_slp_quant_trends_1950_2014_MA$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(rev(my.probs), era5_slp_quant_trends_1950_2014_MA$p_values, lwd = 2, col = "#00A9A5")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_1950_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_1950_2014_SON$p_values, lwd = 2, col = "#00A9A5")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_1950_2014_DJF$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_1950_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_quant_trends_1950_2014_MA$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_quant_trends_1950_2014_MA$p_values, lwd = 2, col = "#00A9A5")
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

# FIGURE S28

# 1950-2014 seasonal intensification rate quantiles vs. p-values

pdf("/Path/to/output/figures/fig_s28_seasonal_intensification_p_values_1950_2014.pdf", 
    family = "Helvetica-Narrow", width = 15, height = 10)
par(mfrow = c(2, 3), mar = c(3, 3, 1, 1), oma = c(3, 7, 7, 1))

# min. SLP SON
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_1950_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_1950_2014_SON$p_values, lwd = 2, col = "#00A9A5")
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
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_1950_2014_DJF$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_1950_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
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
plot(my.probs, crcm5_cmip6_ensemble_slp_intens_quant_trends_1950_2014_MA$p_values, type = "l", 
     col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_slp_intens_quant_trends_1950_2014_MA$p_values, lwd = 2, col = "#00A9A5")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_1950_2014_SON$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_1950_2014_SON$p_values, lwd = 2, col = "#00A9A5")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_1950_2014_DJF$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_1950_2014_DJF$p_values, lwd = 2, col = "#00A9A5")
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
plot(my.probs, crcm5_cmip6_ensemble_wind_intens_quant_trends_1950_2014_MA$p_values, 
     type = "l", col = ("#882255"), lwd = 2, xlim = c(0, 1), ylim = c(0, 1.05),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(my.probs, era5_wind_intens_quant_trends_1950_2014_MA$p_values, lwd = 2, col = "#00A9A5")
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

















