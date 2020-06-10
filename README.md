# Clownfish_persistence

Add an explanation about what this repository does at the top. Also add a "how to interact with this repository if you wanted to run the analyses" bit.

All data were collected by members of the Pinsky lab at sites off Leyte, Philippines from 2012-2018. The raw data are housed XXXX. 

Repository contents and explanation:

**Code** - holds scripts
* *Pull_data_from_database*: pulls data from database and saves as RData files in Data/Data_from_database, pulls output from other analyses from GitHub and saves as RData files in Data/From_other_analyses
* *Constants_database_common_functions*: processes data from database into useable dataframes and saves as RData in Data/Script_outputs, sets constants
  * When run the script, get an error saying "2 failed to parse". Nothing to worry about, just means that the lat/lon coordinates couldn't be matched to two anemone observations because they don't have times associated with them (2 anems AGD found while snorkeling after running out of air in 2018)
* *Total_anems_proportion_hab_sampled*: finds the total number anems at a site (we used metal tagged anems as best estimate of total), calculates proportion habitat sampled at each site in each sampling year, calculates total area sampled across time (for use in egg-recruit survival estimate and as an input into the parentage analyses)
* *Site_widths_and_distances.R*: finds the width of each site, the distance between sites (including possible buffer to site edges to account for larval navigation), and the distance from each site to the northern and southern edges of the sampling area. 
* *TimeSeriesPersistence*: estimates abundance of females at each site through time, fits a mixed-effects model to assess population trend over sampling time period
* *Growth_analysis*: estimates parameters for a von-Bertalanffy growth curve using marked fish captured approximately one year apart
* *Clownfish_encounters*: makes encounter histories for all fish marked, fills in missing sizes by projecting using growth curve (and adding mean size or 0 for years pre-capture)
* *AnemDistFromDiveTrack.R*: finds the lat and lon coordinates of the first capture anemone for each marked fish and the minimum distance from that anemone to the sampling dive tracks each year, in preparation for mark-recapture analyses in ClownfishMarkRecap.R
* *ClownfishMarkRecap.R*: runs mark-recapture analysis using MARK (via RMark) to estimate survival and recapture probabilities
* Metrics_with_uncertainty: estimates the persistence metrics, including best estimates and distributions with uncertainty. Makes plots.

**Data** - holds raw and processed data files and outputs of analyses
* *Data_from_database* (all loaded in Pull_data_from_database.R)
  * anem_db.RData - anemones table from leyte database
  * dives_db.RData - diveinfo table from leyte database
  * fish_db.RData - clownfish table from leyte database
  * fish_seen_db.RData - clown_sightings table from leyte database
  * gps_db.RData - GPX table from leyte database
* *From_other_analyses* (all loaded in Pull_data_from_database.R, either from GitHub or downloaded file)
  * all_offspring.RData - info on all offspring included in parentage analysis 
  * all_parents.RData - info on all parents included in parentage analysis 
  * fish-obs.RData - info on individual fish, tying their genetic and tag ids together 
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
  * females_df_F.RData - raw and estimated scaled number of females at each site through time (produced by TimeSeriesPersistence.R)
  * growth_info_estimate - parameter estimates for von Bertalanffy growth curve (produced by Growth_analysis.R)
  * recap_pairs_year - pairs of sizes for fish recaptured after about a year (produced by Growth_analysis.R)
  * marked_fish - all fish that have been marked either via PIT tag or genotype (so have a fish_indiv) with their anem, dive, etc. info (produced by Clownfish_encounters.R)
  * encounters_list - marked fish re-organized into encounter histories for MARK (produced by Clownfish_encounters.R)
  * encounters_size_means - encounter histories with missing sizes filled in via growth curve projection or mean for pre-capture years (produced by Clownfish_encounters.R)
  * encounters_size_0s - same as encounter_size_means but 0 instead of mean size in a year for the pre-capture years (produced by Clownfish_encounters.R)
  * best_fit_model_dfs - best fit model from mark-recapture analysis and data frames for plotting (produced by ClownfishMarkRecap.R)
  * model_comp_meanYsize_meanYdist - model comparison for mark-recapture models (produced by ClownfishMarkRecap.R)

**Plots** - holds plots
* *TotalAnemsandProbHabSampled* - plots that show anemone counts across years, sites, and methods of assessing anemones; also shows plots of proportion habitat sampled (produced by Total_anems_proportion_hab_sampled.R)
* *Growth* - plots that show sizes of fish recaptured after a year with the von Bertalanffy growth curve fit and distributions of von Bertalanffy parameters (produced by Growth_analysis.R)
* *PhiandpEstimates* - plots showing survival and recapture probabilities (ones starting with surv_ or recap_ are most recent) (produced by ClownfishMarkRecap.R)

