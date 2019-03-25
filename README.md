# Clownfish_persistence

Draft list of cleaned-up scripts in progress:

* Pull_data_from_database: pulls data from database, saves as RData files for other users
* Constants_database_common_functions: processes raw data from database into useable dataframes, sets constants, reads in data from other analyses
  * When run the script, get an error saying "2 failed to parse". Nothing to worry about, just means that the lat/lon coordinates couldn't be matched to two anemone observations because they don't have times associated with them (2 anems AGD found while snorkeling after running out of air in 2018)
* Total_anems_proportion_hab_sampled: finds the total number anems at a site (we used tagged anems as best estimate of total), calculates proportion habitat sampled at each site in each sampling year, calculates total area sampled across time (for use in egg-recruit survival estimate)
* NEED TO INSERT INFO ABOUT MARK-RECAP SURVIVAL SCRIPT AND GROWTH SCRIPT HERE!
* Metrics_with_uncertainty: estimates the persistence metrics, including best estimates and distributions with uncertainty. Makes plots.
