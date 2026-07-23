# MMW-JGRAtmo
This repository contains the tracked cyclone datasets and analysis code necessary to re-create the results presented in Weaver, M.M. et al. (2026), Historical nor'easter intensity trends in reanalysis and regional climate simulations: A cross-dataset perspective.
## Related Publication
Weaver, M.M., D.K Apoznanski, O. Asselin, A.J. Broccoli, A.J. Garner, A.L. Horton, M.P. Labonte, C.J. Larson, M. Leduc, D. Paquin, J.L. Seidenstein, & M.E. Mann (2026), Historical nor'easter intensity trends in reanalysis and regional climate simulations: A cross-dataset perspective.
## Files Descriptions
### Storm Track Data (located in the /Data/storm_tracks directory):
* Contains subdirectories for each dataset (CFSR, ERA5, and CRCM5).
* The CRCM5 subdirectory contains folders for each dynamically downscaled historical simulation (CRCM5-ERA5, CRCM5-CanESM5, CRCM5-CNRM, CRCM5-MPI, and CRCM5-NorESM2).
* The final folder in each directory contains `YYYY_stormtracks_unfiltered.csv`, which contains storm data for all cyclones tracked in a given nor'easter season for the specified simulation (where the season is September to April, and YYYY is the end year of that nor'easter season).
### Wind Data (located in the /Data/storm_wind directory):
* Contains subdirectories for each dataset that requires separate wind data files (CFSR and CRCM5 only — ERA5 wind data is included in the storm track CSVs).
* The CRCM5 subdirectory contains folders for each dynamically downscaled historical simulation (CRCM5-ERA5, CRCM5-CanESM5, CRCM5-CNRM, CRCM5-MPI, and CRCM5-NorESM2).
* The final folder in each directory contains `YYYY_stormwinds_MODEL.csv`, which contains storm data for all cyclones tracked in a given nor'easter season for the specified simulation (where YYYY is the end year of that nor'easter season, and MODEL is the simulation name (CFSR, ERA5, CanESM5, CNRM, MPI, or NorESM2).
### Code (located in the /Scripts directory):
* `cyclone_tracking_algorithm.py` contains all functions and usage example for the cyclone tracking algorithm adapted from Chen et al. (2025; see the related publication's references)
* `historical_noreaster_analysis_functions.R` contains all functions used in hist_noreaster_intensity_analysis.R` and `TC_filter_hist_noreaster_intensity.R` (see below).
* `hist_noreaster_intensity_analysis.R` contains all analysis code necessary to generate main text Figures 2-8 and Supplementary Figures S1-2 & S7-16, which includes analyses of historical nor'easter intensity and intensification rate quantiles, track locations, and peak intensity/intensification rate locations.
* `figures_hist_nor'easter_intensity.R` contains the code that creates main text Figures 2-8 and Supplementary Figures S1-2 & S7-16.
* `TC_filter_hist_noreaster_intensity.R` contains all analysis code necessary to generate Table 1, main text Figure 1, and Supplementary Figures S3-6 & S10, which includes functionality testing of the tropical cyclone filter applied to the tracked cyclone datasets.
