# Clownfish_persistence

This repository provides the data, scripts, and figures for the analyses in "Persistence of a reef fish metapopulation via network connectivity: theory and data," which assesses metapopulation persistence for yellowtail anemonefish (*Amphiprion clarkii*) on a set of patch reefs along the coast of Leyte, Philippines.

Authors: 
Allison G. Dedrick<sup>a,b</sup>, Katrina A. Catalano<sup>a</sup>, Michelle R. Stuart<sup>a</sup>, J. Wilson White<sup>c</sup>, Humberto R. Montes, Jr.<sup>d</sup>, Malin L. Pinsky<sup>a</sup>

<sup>a</sup>Department of Ecology, Evolution, and Natural Resources, Rutgers University, 14 College Farm Road, New Brunswick, NJ 08901 USA<br>
<sup>b</sup>Corresponding author: adedrick@stanford.edu. Current address: Stanford Woods Institute for the Environment, Stanford University, Stanford, CA 94305 USA<br>
<sup>c</sup>Department of Fisheries and Wildlife, Coastal Oregon Marine Experiment Station, Oregon State University, Newport, OR 97365 USA<br>
<sup>d</sup>Visayas State University, Pangasugan, Baybay City, 6521 Leyte, Philippines<br>

All data were collected by members of the Pinsky lab at sites off Leyte, [Philippines from 2012-2018](https://github.com/pinskylab/PhilippinesSurveyData). The raw data are housed in the folder Data/Data_from_database, with a brief description in Metadata.md in the Data folder. Access to the BCO-DMO data can be found [here](https://www.bco-dmo.org/dataset/862415) or with the citation found below.

To run the analyses, use the script Metrics_with_uncertainty_site_specific.R, which will source outputs from analyses done in the other scripts. To recreate the figures, use the script Create_results_figures.R, which will source outputs from Metrics_with_uncertainty_site_specific.R.

The figures in the manuscript are located in Plots/
* Figure 1: LifeCycleSchematic/metrics_life_cycle_schematics.pdf
* Figure 2: FigureDrafts/Map_and_photo_2.pdf
* Figure 3: FigureDrafts/Parameter_inputs.pdf
* Figure 4: FigureDrafts/Abundance_LEP_LRP_LocalReplacement_FreqPlots.pdf
* Figure 5: FigureDrafts/SP_NP_connMatrixR_freq.pdf
* Figure 6: Figure/Drafts/What_if_4_panels_3D.pdf
* Appendix Figure D.1: Schematic/Schematic.pdf
* Appendix Figure D.2: UndercountingRecruits/RecruitScalingSchematic.pdf
* Appendix Figure D.3: FigureDrafts/APP_FIG_Parameter_inputs.pdf
* Appendix Figure D.4: FigureDrafts/Prop_of_kernel_sampled_by_site.pdf
* Appendix Figure D.5: FigureDrafts/APP_FIG_surv_by_size_and_site_Phisiteplussize_psizeplusdist.pdf
* Appendix Figure D.6: FigureDrafts/APP_FIG_recap_effects_Phisiteplussize_psizeplusdist.pdf
* Appendix Figure D.7: FigureDrafts/Time_series_scaled_F_by_site_with_lines.pdf
* Appendix Figure D.8: FigureDrafts/LRP_by_site.pdf
* Appendix Figure D.9: FigureDrafts/APP_FIG_LRP_LocalReplacement_withoutDDconsidered_freq.pdf
* Appendix Figure D.10: FigureDrafts/APP_FIG_SP_NP_connMatrixR_withoutDDcompensation_freq.pdf
* Appendix Figure D.11: FigureDrafts/LEP_uncertainty_by_param.pdf
* Appendix Figure D.12: FigureDrafts/LRP_uncertainty_by_param.pdf
* Appendix Figure D.13: FigureDrafts/RperE_uncertainty_by_param.pdf
* Appendix Figure D.14: FigureDrafts/NP_uncertainty_by_param.pdf

Repository contents and explanation:

**Code** - holds scripts
* *Pull_data_from_database.R*: pulls data from database and saves as RData files in Data/Data_from_database, pulls output from other analyses (like the parentage analysis and dispersal kernel fitting) from GitHub and saves as RData files in Data/From_other_analyses
* *Constants_database_common_functions.R*: processes data from database into useable dataframes and saves as RData in Data/Script_outputs, sets constants
  * When run the script, get two errors saying "2 failed to parse" and "3 failed to parse". Nothing to worry about, just means that the lat/lon coordinates couldn't be matched to two anemone observations because they don't have times associated with them (2 anems AGD found while snorkeling after running low on air in 2018) (and then one more when all anems, including those without ids, are included)
* *Total_anems_proportion_hab_sampled.R*: finds the total number anems at a site (we used metal tagged anems as best estimate of total), calculates proportion habitat sampled at each site in each sampling year, calculates total area sampled across time (for use in egg-recruit survival estimate and as an input into the parentage analyses)
* *Site_widths_and_distances.R*: finds the width of each site, the distance between sites (including possible buffer to site edges to account for larval navigation), and the distance from each site to the northern and southern edges of the sampling area. 
* *Density_dependence_scaling.R*: finds the number of anemones occupied by *A. clarkii* (APCL), occupied by other clownfish, and unoccupied during the two seasons when anemones were reasonably comprehensively surveyed, 2015 winter and 2017, for use in scaling up recruits to account for effects of density-dependence in settlement
* *Growth_analysis.R*: estimates parameters for a von-Bertalanffy growth curve using marked fish captured approximately one year apart
* *Clownfish_encounters.R*: makes encounter histories for all fish marked, fills in missing sizes by projecting using growth curve (and adding mean size or 0 for years pre-capture)
* *AnemDistFromDiveTrack.R*: finds the lat and lon coordinates of the first capture anemone for each marked fish and the minimum distance from that anemone to the sampling dive tracks each year, in preparation for mark-recapture analyses in ClownfishMarkRecap.R
* *ClownfishMarkRecap.R*: runs mark-recapture analysis using MARK (via RMark) to estimate survival and recapture probabilities
* *TimeSeriesPersistence.R*: estimates abundance of females at each site through time, fits a mixed-effects model to assess population trend over sampling time period
* *Metrics_with_uncertainty_site_specific.R*: estimates the persistence metrics, including best estimates and distributions with uncertainty
* *Create_results_figures.R*: uses saved outputs from Metrics_with_uncertainty_site_specific and other scripts to create figures

**Data** - holds raw and processed data files and outputs of analyses (any files not within a subfolder in the Data folder are old unused versions or files from abandoned analyses and should be ignored).
* *Data_from_database* (all loaded in Pull_data_from_database.R)
  * anem_db.RData - anemones table from leyte database
  * dives_db.RData - diveinfo table from leyte database
  * fish_db.RData - clownfish table from leyte database
  * fish_seen_db.RData - clown_sightings table from leyte database
  * gps_db.RData - GPX table from leyte database
* *Database_backups* - holds backups of older versions of database data
* *From_other_analyses* (all loaded in Pull_data_from_database.R, either from GitHub or downloaded file)
  * all_offspring.RData - info on all offspring included in parentage analysis (from [Catalano et al. 2020](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.15732))
  * all_parents.RData - info on all parents included in parentage analysis (from [Catalano et al. 2020](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.15732))
  * fish-obs.RData - info on individual fish, tying their genetic and tag ids together 
  * k_theta_allyear_95CI_values.RData - 95% confidence bounds on k and theta distributions for dispersal kernel (from [Catalano et al. 2020](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.15732))
  * kernel_summary.RData - dispersal kernel fit information (from [Catalano et al. 2020](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.15732))
  * length_count8llEA.RData - egg count and female length model fit (work done by Adam Yawdoszyn)
* *Map_data/Site_hulls* - shapefiles of the sites, used in SiteMap.R to make maps
* *Script_outputs*
  * site_vec_order.RData - vector of site names and their alphabetical order and geographic order (north to south) (produced by Constants_database_common_functions.R)
  * allfish_caught.RData - all APCL caught, with dive and anemone info appended (produced by Constants_database_common_functions.R)
  * dives_processed.RData - number of dives at each site in each year (produced by Constants_database_common_functions.R)
  * dives_processed_clownfish.RData - number of dives during clownfish sampling season (produced by Constants_database_common_functions.R)
  * site_visits.RData - gives a 1 or NA for each site and year, indicating whether it was visited during the clownfish sampling field season (produced by Constants_database_common_functions.R)
  * anems_Processed.RData - anems with anem_ids, dive info and times processed and such so can link to gps (produced by Constants_database_common_functions.R)
  * anems_Processed_all.RData - times processed for all anems, even without ids (produced by Constants_database_common_functions.R)
  * all_parents_by_site.RData - number of parents by site (produced by Constants_database_common_functions.R)
  * site_edge_anems.RData - anem_ids for anemones at the northern edge, southern edge, and middle of each site (eyeballed from QGIS) (produced by Constants_database_common_functions.R)
  * anems_visited_by_year.RData - number of anemones visited at each site each year and total anems at each site (produced by Total_anems_proportion_hab_sampled.R)
  * cumulative_prop_hab_sampled_by_site.RData - cumulative proportion habitat sampled at each site through time (produced by Total_anems_proportion_hab_sampled.R), just considering metal-tagged anemones as anemones sampled and total anemones at a site
  * cumulative_prop_hab_sampled_by_site_v2.RData - same as cumulative_prop_hab_sampled_by_site but includes plastic-tagged anemones for the sites Caridad Proper, Sitio Lonas, and Sitio Tugas (which don't have metal-tagged anemones)
  * prop_hab_appendix_table.RData - has total anemones and proportion sampled each year for each site and overall region info for appendix table in paper (produced by Total_anems_proportion_hab_sampled.R)
  * site_width_info.RData - width of each site and distance to northern and southern edges of sampling area (produced by Site_widths_and_distances.R)
  * site_dist_info.RData - distances between each pair of sites (produced by Site_widths_and_distances.R)
  * sampling_area_edges.RData - lat and lon coordinates of northern and southern limits of sampling area (produced by Site_widths_and_distances.R)
  * site_buffer_info.RData - provides maximum buffer to N and S of each site without overlapping shadows of other sites, for use in accounting for larval navigation in distances for kernel integration (produced by Site_widths_and_distances.R)
  * anems_APCL_and_not.RData - gives proportion anemones occupied by APCL and unoccupied by any clownfish averaged between 2015 and 2017 (produced by Density_dependent_scaling.R)
  * anems_APCL_and_not_by_year.RData - gives numbers and proportion of anemones occupied by APCL, other clownfish, and unoccupied in 2015 and 2017 (produced by Density_dependent_scaling.R)
   * growth_info_estimate.RData - parameter estimates for von Bertalanffy growth curve (produced by Growth_analysis.R)
  * recap_pairs_year.RData - pairs of sizes for fish recaptured after about a year (produced by Growth_analysis.R)
  * marked_fish.RData - all fish that have been marked either via PIT tag or genotype (so have a fish_indiv) with their anem, dive, etc. info (produced by Clownfish_encounters.R)
  * encounters_list.RData - marked fish re-organized into encounter histories for MARK (produced by Clownfish_encounters.R)
  * encounters_size_means_by_year.RData - encounter histories with missing sizes filled in via growth curve projection or mean for pre-capture years (produced by Clownfish_encounters.R)
  * encounters_size_0s.RData - same as encounter_size_means_by_year but 0 instead of mean size in a year for the pre-capture years (produced by Clownfish_encounters.R)
  * min_survey_dist_to_anems.RData - minimum distance from dive GPS tracks each year to capture anems for tagged fish (produced by AnemDistFromDiveTrack.R)
  * encounters_dist.RData - encounter histories with distances to dive tracks (produced by AnemDistFromDiveTrack.R)
  * encounters_dist_mean.RData - same as encounter_dist but with overall mean distance filled in for NAs (produced by AnemDistFromDiveTrack.R)
  * encounters_dist_mean_by_year.RData - same as encounter_dist but with mean distance by year filled in for NAs (produced by AnemDistFromDiveTrack.R)
  * best_fit_model_dfs - best fit model from mark-recapture analysis and data frames for plotting (produced by ClownfishMarkRecap.R)
  * model_comp_meanYsize_meanYdist - model comparison for mark-recapture models (produced by ClownfishMarkRecap.R)
  * females_df_F.RData - raw and estimated scaled number of females at each site through time (produced by TimeSeriesPersistence.R)
  * output_uncert files - output from 1000 metric calculations with uncertainty, list the uncertainty considered (output_uncert_all includes all uncertainty, output_uncert_growth includes uncertainty only in growth), DD indicates that density dependence was compensated for in the early life stages (produced by Metrics_with_uncertainty_site_specific.R)
  * param_set_full.RData - full set of parameters for 1000 metric calculations (produced by Metrics_with_uncertainty_site_specific.R)
  * params_summary.RData - parameter summary to make writing easier (produced by Metrics_with_uncertainty_site_specific.R)
  * metrics_summary.RData - metrics summary to make writing easier (produced by Metrics_with_uncertainty_site_specific.R)
  * perc_hab and wider_region and larv_nav files - results of the simulations for sensitivity to percent of the region that is habitat, region width, and larval navigation (produced by Metrics_with_uncertainty_site_specific.R)
  
**Plots** - holds plots

List of manuscript figures:
* Figure 1: LifeCycleSchematic/metrics_life_cycle_schematics.pdf
* Figure 2: FigureDrafts/Map_and_photo_2.pdf
* Figure 3: FigureDrafts/Parameter_inputs.pdf
* Figure 4: FigureDrafts/Abundance_LEP_LRP_LocalReplacement_FreqPlots.pdf
* Figure 5: FigureDrafts/SP_NP_connMatrixR_freq.pdf
* Figure 6: Figure/Drafts/What_if_4_panels_3D.pdf
* Appendix Figure D.1: Schematic/Schematic.pdf
* Appendix Figure D.2: UndercountingRecruits/RecruitScalingSchematic.pdf
* Appendix Figure D.3: FigureDrafts/APP_FIG_Parameter_inputs.pdf
* Appendix Figure D.4: FigureDrafts/Prop_of_kernel_sampled_by_site.pdf
* Appendix Figure D.5: FigureDrafts/APP_FIG_surv_by_size_and_site_Phisiteplussize_psizeplusdist.pdf
* Appendix Figure D.6: FigureDrafts/APP_FIG_recap_effects_Phisiteplussize_psizeplusdist.pdf
* Appendix Figure D.7: FigureDrafts/Time_series_scaled_F_by_site_with_lines.pdf
* Appendix Figure D.8: FigureDrafts/LRP_by_site.pdf
* Appendix Figure D.9: FigureDrafts/APP_FIG_LRP_LocalReplacement_withoutDDconsidered_freq.pdf
* Appendix Figure D.10: FigureDrafts/APP_FIG_SP_NP_connMatrixR_withoutDDcompensation_freq.pdf
* Appendix Figure D.11: FigureDrafts/LEP_uncertainty_by_param.pdf
* Appendix Figure D.12: FigureDrafts/LRP_uncertainty_by_param.pdf
* Appendix Figure D.13: FigureDrafts/RperE_uncertainty_by_param.pdf
* Appendix Figure D.14: FigureDrafts/NP_uncertainty_by_param.pdf

Notes on other plots:
* *AnemLocations* - old plots looking at GPS locations of anems sampled multiple times
* *DataCharacteristics* - data summary plots, not necessarily up to date
* *DistFromTrackToAnems/Distance_to_marked_fish_capture_anems_allyears.pdf* - histograms by year of minimum distance from divers to capture anems for marked fish (produced by AnemDistFromDiveTrack.R)
* *FigureDrafts* - contains other versions of results plots
* *TotalAnemsandProbHabSampled* - plots that show anemone counts across years, sites, and methods of assessing anemones; also shows plots of proportion habitat sampled (produced by Total_anems_proportion_hab_sampled.R)
* *Growth* - plots that show sizes of fish recaptured after a year with the von Bertalanffy growth curve fit and distributions of von Bertalanffy parameters (produced by Growth_analysis.R)
* *PhiandpEstimates* - plots showing survival and recapture probabilities (ones starting with surv_ or recap_ are most recent) (produced by ClownfishMarkRecap.R)
* *LifeCycleSchematic*, *Schematic*, *UndercountingRecruits* - schematic figures

## BCO-DMO Data Accessibility
Pinsky, M., Stuart, M. (2023) Temperature loggers (HOBO) placed in two locations off the coast of the West coast of Leyte, the Philippines , 2012-2019. Biological and Chemical Oceanography Data Management Office (BCO-DMO). (Version 2) Version Date 2023-03-06. doi:10.26008/1912/bco-dmo.862415.2
