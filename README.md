# Clownfish_persistence

Add an explanation about what this repository does at the top. Also add a "how to interact with this repository if you wanted to run the analyses" bit.

All data were collected by members of the Pinsky lab at sites off Leyte, Philippines from 2012-2018. The raw data are housed XXXX. 

Repository contents and explanation:

**Code** - holds scripts
* Pull_data_from_database: pulls data from database and saves as RData files in Data/Data_from_database, pulls output from other analyses from GitHub and saves as RData files in Data/From_other_analyses
* Constants_database_common_functions: processes data from database into useable dataframes and saves as RData in Data/Script_outputs, sets constants
  * When run the script, get an error saying "2 failed to parse". Nothing to worry about, just means that the lat/lon coordinates couldn't be matched to two anemone observations because they don't have times associated with them (2 anems AGD found while snorkeling after running out of air in 2018)
* Total_anems_proportion_hab_sampled: finds the total number anems at a site (we used tagged anems as best estimate of total), calculates proportion habitat sampled at each site in each sampling year, calculates total area sampled across time (for use in egg-recruit survival estimate)
* NEED TO INSERT INFO ABOUT MARK-RECAP SURVIVAL SCRIPT AND GROWTH SCRIPT HERE!
* Metrics_with_uncertainty: estimates the persistence metrics, including best estimates and distributions with uncertainty. Makes plots.

**Data** - holds raw and processed data files and outputs of analyses
* *Data_from_database* (all loaded in Pull_data_from_database.R)
  * anem_db.RData - anemones table from leyte database
  * dives_db.RData - diveinfo table from leyte database
  * fish_db.RData - clownfish table from leyte database
  * fish_seen_db.RData - clown_sightings table from leyte database
  * gps_db - GPX table from leyte database
* *From_other_analyses* (all loaded in Pull_data_from_database.R, either from GitHub or downloaded file)
  * all_offspring - info on all offspring included in parentage analysis 
  * all_parents - info on all parents included in parentage analysis 
  * fish-obs - info on individual fish, tying their genetic and tag ids together 
* *Script_outputs*

