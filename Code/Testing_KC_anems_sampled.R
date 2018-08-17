# Testing how KC pulls anems to be sampled (see emails 7/17/18 with subject "anems visited in 2012")

leyte <- read_db("Leyte")

allfish_fish <- leyte %>% 
  tbl("clownfish") %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, cap_id, anem_table_id, recap, tag_id, color, size) %>%
  collect() %>%
  filter(fish_spp == "APCL")

######## this is the code in her email...
allfish_anems <- leyte %>%
  tbl("anemones") %>%
  select(anem_table_id, dive_table_id, anem_obs, obs_time, anem_id, anem_obs, old_anem_id) %>%
  collect() %>%
  filter(anem_table_id %in% allfish_fish$anem_table_id)

allfish_dives <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, dive_type, date, site, gps) %>%
  collect() %>%
  filter(dive_table_id %in% allfish_anems$dive_table_id)

all_anems <- left_join(allfish_anems, allfish_dives, by="dive_table_id") 
all_anems <- all_anems %>% mutate(year = as.integer(substring(all_anems$date,1,4))) # this line not in the code KC sent me

n_anemones_sampled <- all_anems %>%
  group_by(site, year) %>%
  filter(anem_table_id %in% all_anems$anem_table_id) %>%
  summarise(n_anem_sampled= n())


