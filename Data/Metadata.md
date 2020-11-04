# Metadata descriptions for data files in Data_from_database

## Overview of sampling methods

Members of the Pinsky lab at Rutgers University and colleagues sampled a metapopulation of yellowtail anemonefish (*Amphiprion clarkii*) annually from 2012-2018. Sampling focused
on a set of nineteen reef patches spanning 30 km along the western coast of Leyte island facing the Camotes Sea in the Philippines. Divers using SCUBA and tethered to GPS readers swam the extent of each patch to visit anemones occupied by yellowtail anemonefish, sampling fish and habitat in most patches each year. At each anemone, the divers caught fish 3.5 cm and larger, took a tissue sample, measured fork length, and noted tail color as an indicator of life stage. Starting in 2015, fish 6.0 cm and larger were also tagged with a passive integrated transponder (PIT) tag unless already tagged. 

## File descriptions
* anem_db.RData - holds data on visits to anemones, including info on anemone species, depth, size, and who took the data. 
    * links to the dive data table through dive_table_id, which contains info on site and date
    * one entry per visit to an anemone, each with a unique anem_table_id
* dives_db.RData - data on each dive, including site, number of the gps units used, the start and stop time, and conditions
    * links to the anemone data table through dive_table_id
    * links to the gps data table through the gps unit number and the time
    * one entry per diver with gps unit per dive, each with a unique dive_table_id
* fish_db.RData - data on yellowtail anemonefish that were captured, including fish size, color, sex, sample id for genetic info, and tag id for PIT tag number
    * links to the anemone data table through anem_table_id, must link to both anemone table and dive table to get info on site and date of capture
    * fish-obs.RData in the folder Data/From_other_analyses links multiple captures of the same fish and assigns each fish a unique id, fish_indiv
* fish_seen_db.RData - data on anemonefish that were seen but not captured during sampling dives
    * links to the anemone data table through anem_table_id
    * includes sightings of yellowtail anemonefish and other species of anemonefish
    * includes size of fish estimated visually by the divers
* gps_db.RData - gps positions collected by the gps units worn by divers during sampling dives
    * links to dive data table through unit and time
    * once linked to the correct dive, time can link to anemone and fish data tables
