# Look at egg-recruit relationship (in attempts to estimate survival from eggs to recruits)

#################### Set-up: ####################
# Load relevant libraries
library(RCurl) #allows running R scripts from GitHub
library(RMySQL) #might need to load this to connect to the database?
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(here)

# Set constants
min_size_breedingF <- 6.0 # minimum size for filtering out breeding females
min_size_breedingM <- 5.0 # minumum size for filtering out breeding males (totally made up right now)
eggs_per_clutch = 1763 #from LEP_calc_WillWhite.R
clutch_per_year = 11.9 #from LEP_calc_WillWhite.R

min_size_recruit <- 3.5 #minimum size to be considered a recruit from an egg produced the previous year
max_size_recruit <- 6.0 #maximum size to be considered a recruit from an egg produced the previous year

sample_months <- c(3,4,5,6,7,8) #set months to avoid winter 2015 anem surveys
dive_list <- c("0","C","D","E")

# Load other files
# load file with proportion habitat sampled estimates (from TotalAnemsAtSite.R)
load(file=here("Data",'anem_sampling_table.RData')) #file with anems, after 2018 anems matched

# prob of catching a fish by site, from KC script: https://github.com/katcatalano/parentage/blob/master/notebooks/proportion_sampled_allison.ipynb
prob_r <- c(0.5555556, 0.2647059, 0.8888889, 0.6666667, 0.2000000, #2016 recapture dives
            0.8333333, 0.4666667, 0.2000000, 0.8333333, 1.0000000, #2017 recapture dives
            0.3333333, 0.5789474, 0.6250000, 0.4090909) #2018 recapture dives

#################### Functions: ####################
# Functions and constants from my GitHub function/constant collection
# script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQAXExRItr2tnepmkRr7NRt4hylZrks5aciBtwA%3D%3D", ssl.verifypeer = FALSE)
# eval(parse(text = script))
script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQJT5uCEwjgDGYOneY0W6Zdjol5axks5alHmBwA%3D%3D", ssl.verifypeer = FALSE)
eval(parse(text = script))

# Functions from Michelle's GitHub helpers script
#field_helpers (similar idea to helpers, in field repository) - this might be the newer version of the helpers collection?
script <- getURL("https://raw.githubusercontent.com/pinskylab/field/master/scripts/field_helpers.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#################### Running things: ####################
########## Pull info from database
leyte <- read_db("Leyte")

allfish_fish <- leyte %>% 
  tbl("clownfish") %>%
  select(fish_table_id, anem_table_id, fish_spp, sample_id, cap_id, anem_table_id, recap, tag_id, color, size) %>%
  collect() %>%
  filter(fish_spp == "APCL")

allfish_anems <- leyte %>%
  tbl("anemones") %>%
  select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id) %>%
  collect() %>%
  filter(anem_table_id %in% allfish_fish$anem_table_id)

allfish_dives <- leyte %>%
  tbl("diveinfo") %>%
  select(dive_table_id, dive_type, date, site, gps) %>%
  collect() %>%
  filter(dive_table_id %in% allfish_anems$dive_table_id)

# pull out just the year and put that in a separate column
allfish_dives$year <- as.integer(substring(allfish_dives$date,1,4))
allfish_dives$month <- as.integer(substring(allfish_dives$date,6,7))

# join together
allfish <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
allfish <- left_join(allfish, allfish_dives, by="dive_table_id")

allfish$size <- as.numeric(allfish$size) #make size numeric (rather than a chr) so can do means and such

# try to reduce multiple observations of same fish - filter out Jan-Feb 2015 anem surveys
allfish_abbrev <- allfish %>% filter(month %in% sample_months)
allfish_abbrev_divetype <- allfish %>% filter(month %in% sample_months) %>% filter(dive_type %in% dive_list)

#pull out just tagged fish
taggedfish <- allfish %>% filter(!is.na(tag_id))

########## Estimate number of eggs coming from each site and the whole population each year
##### Find number of breeding females
# Find number of breeding females seen or captured
breedingF <- allfish %>% 
  group_by(site,year) %>%
  filter(color %in% c("Y","YP")) %>% #just filter out YP fish (and Y for 2012, when there was no YP)
  filter(size >= min_size_breedingF) %>% #but also make sure they're big enough to be a breeder (sometimes small fish get marked YP too)
  summarise(NbreedingF = n())

# Find number of breeding females tagged
breedingF_tagged <- allfish %>%
  group_by(year,site) %>%
  filter(color %in% c("Y","YP")) %>%
  filter(!is.na(tag_id)) %>%
  summarise(Ntagged_breedingF = n())

# Combine - since didn't tag all years, might not have any females for some site/years
breedingF_info <- left_join(breedingF, breedingF_tagged, by = c("site" = "site", "year" = "year")) # no 2012 obs here - why?

# Find average prob_r from KC code
prob_r_avg <- mean(prob_r)

# Join proportion hab sampled with breeding F info
breedingF_info <- left_join(breedingF_info, (anems_table %>% select(site, year, prop_hab_sampled, prop_hab_sampled_2)), by = c("site", "year"))

# Add one more column that is prop_hab_sampled_combo - prop_hab_sampled_2 for 2012, prop_hab_sampled (using anem_ids) for 2013-2018
breedingF_info$prop_hab_sampled_combo <- rep(NA, length(breedingF_info$site))

for(i in 1:length(breedingF_info$site)) {
  if(breedingF_info$year[i] == 2012) {
    breedingF_info$prop_hab_sampled_combo[i] = breedingF_info$prop_hab_sampled_2[i] #use the non-anem_id method for 2012
  } else {
    breedingF_info$prop_hab_sampled_combo[i] = breedingF_info$prop_hab_sampled[i] #and the anem_id method for the other years
  } 
}

# Add a breedingF combo that uses the tagged fish for 2015-2018 and recorded fish of the right size/color before that
breedingF_info$NbreedingF_combo <- rep(NA, length(breedingF_info$site))

for(i in 1:length(breedingF_info$site)) {
  if(breedingF_info$year[i] == 2012 | breedingF_info$year[i] == 2013 | breedingF_info$year[i] == 2014) {
    breedingF_info$NbreedingF_combo[i] = breedingF_info$NbreedingF[i] #use fish w/right color/size recorded pre-tags
  } else {
    breedingF_info$NbreedingF_combo[i] = breedingF_info$Ntagged_breedingF[i] #and tagged females 2015-2018
  }
}

# Estimate number of breeding females in each pop - scale up N captured or tagged by proportion of habitat sampled and prob of capture
breedingF_info$totalF_est <- breedingF_info$NbreedingF_combo/(breedingF_info$prop_hab_sampled_combo*prob_r_avg)

# Estimate total number of eggs produced
breedingF_info$est_eggs <- breedingF_info$totalF_est*eggs_per_clutch*clutch_per_year #using scaled-up estimates of females
breedingF_info$rawF_egg_est <- breedingF_info$NbreedingF_combo*eggs_per_clutch*clutch_per_year #using raw numbers of captured/tagged breeding females

##### Find number of recruits in a pop in each year
# Find number of recruits (3.5cm - 6.0cm) seen or captured
recruits <- allfish %>% 
  group_by(site,year) %>%
  filter(size >= min_size_recruit & size < max_size_recruit) %>% #define "recruits" as fish between 3.5-6.0cm
  filter(dive_type %in% dive_list) %>% #avoid anem and recapture dives to avoid double-counting?
  summarise(Nrecruits = n())

# Find number of recruits clipped - way fewer numbers, what is happening here?
recruits_clipped <- allfish %>%
  group_by(year,site) %>%
  filter(size < min_tag_size) %>%
  filter(!is.na(sample_id)) %>%
  summarise(Nrecruits_clipped = n())

# Combine - looks like really not very many clipped so going to go with the other column but should check with Michelle/Malin what it is!
recruits_info <- left_join(recruits, recruits_clipped, by = c("site" = "site", "year" = "year")) 

# And join with habitat info, using the combo column in breedingF_info
recruits_info <- left_join(recruits_info, (breedingF_info %>% select(year, site, prop_hab_sampled_combo)), by = c("site" = "site", "year" = "year"))

# Now scale up, like did with females - scale up N captured or tagged by proportion of habitat sampled and prob of capture
recruits_info$totalR_est <- recruits_info$Nrecruits/(recruits_info$prop_hab_sampled_combo*prob_r_avg)
recruits_info$egg_year <- recruits_info$year - 1

# Create a new data frame with both females/eggs and recruits a year later, for easier plotting
females_recruits_summary <- breedingF_info %>% select(site, year, NbreedingF_combo, totalF_est, est_eggs, rawF_egg_est, prop_hab_sampled_combo)
females_recruits_summary <- left_join((breedingF_info %>% select(site, year, NbreedingF_combo, totalF_est, est_eggs, rawF_egg_est, prop_hab_sampled_combo)), 
                                      (recruits_info %>% select(site, egg_year, Nrecruits, Nrecruits_clipped, totalR_est)), 
                                      by = c("site" = "site", "year" = "egg_year"))
females_recruits_summary <- females_recruits_summary %>% filter(year %in% c(2012,2013,2014,2015,2016,2017)) #since don't have recruits for 2018 female output...

# And one more data frame where each year all females, eggs, and recruits are summed up across sites
metapop_level <- females_recruits_summary %>%
  group_by(year) %>% 
  summarize(metapop_rawF = sum(NbreedingF_combo, na.rm =TRUE), 
            metapop_estF = sum(totalF_est, na.rm = TRUE),
            metapop_est_eggs = sum(est_eggs, na.rm = TRUE),
            metapop_rawR = sum(Nrecruits, na.rm = TRUE),
            metapop_estR = sum(totalR_est, na.rm = TRUE))

########## Find linear relationships
### Estimated eggs and recruits each site individually
egg_recruits_est_mod1 <- lm(totalR_est ~ est_eggs, data=females_recruits_summary) #estimate linear relationship
egg_recruits_est_mod1_predicted <- data.frame(recruits_pred = predict(egg_recruits_est_mod1, females_recruits_summary), est_eggs = females_recruits_summary$est_eggs) #predict values
females_recruits_summary_mod1 <- left_join(females_recruits_summary, egg_recruits_est_mod1_predicted, by = "est_eggs") #add predictions to data frame to plot

### Estimated eggs and recruits summed up to metapop level
meta_eggs_recruits_est_mod1 <- lm(metapop_estR ~ metapop_est_eggs, data=metapop_level) #estimate linear relationship
meta_egg_recruits_est_mod1_predicted <- data.frame(recruits_pred = predict(meta_eggs_recruits_est_mod1, metapop_level), metapop_est_eggs = metapop_level$metapop_est_eggs) #predict values
metapop_level_mod1 <- left_join(metapop_level, meta_egg_recruits_est_mod1_predicted, by = "metapop_est_eggs") #add predictions to data frame to plot


#################### Plots: ####################
# Number of eggs produced by site with number of recruits there a year later, at site level
pdf(file = here("Plots/EggRecruitRelationship", "Eggs_recruits_by_site.pdf"))
ggplot(data = females_recruits_summary_mod1, aes(x=est_eggs, y=totalR_est, color=year, shape=site)) +
  geom_point(size=3) +
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) +
  geom_line(color = "black", data=females_recruits_summary_mod1, aes(x=est_eggs, y=recruits_pred)) +
  xlab("estimated egg output") + ylab("estimated number of recruits following year") +
  ggtitle("Egg-recruit scatter plot for each site and year individually") +
  theme(text = element_text(size=25)) +
  theme(axis.text.x = element_text(size=20)) + theme(axis.text.y = element_text(size=20)) +
  theme_bw() #can do theme_bw(base_size = 18) to change size of all elements at once while keeping them proportional to each other
dev.off()

# Number of eggs produced by site with number of recruits there a year later, at site level, ESA presentation
pdf(file = here("Plots/EggRecruitRelationship", "Eggs_recruits_by_site_ESA.pdf"))
ggplot(data = females_recruits_summary_mod1, aes(x=est_eggs, y=totalR_est, color=year)) +
  geom_point(size=4) +
  #scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) +
  geom_line(color = "black", data=females_recruits_summary_mod1, aes(x=est_eggs, y=recruits_pred)) +
  xlab("estimated egg output") + ylab("estimated recruits") +
  #ggtitle("Egg-recruit scatter plot for each site and year individually") +
  theme_bw() +
  theme(text = element_text(size=40)) +
  theme(axis.text.x = element_text(size=30, angle=90, vjust=0.5)) + theme(axis.text.y = element_text(size=30)) +
  theme(legend.text=element_text(size=20)) + theme(legend.title=element_text(size=25))
  #theme_bw(base_size = 20)
  #theme(text = element_text(size=25)) +
  #theme(axis.text.x = element_text(size=20)) + theme(axis.text.y = element_text(size=20)) +
  #theme_bw() #can do theme_bw(base_size = 18) to change size of all elements at once while keeping them proportional to each other
dev.off()

# Number of eggs produced by site with number of recruits there a year later, at site level, same as above but for MPE poster
pdf(file = here("Plots/EggRecruitRelationship", "Eggs_recruits_by_site_MPEposter.pdf"))
ggplot(data = females_recruits_summary, aes(x=est_eggs, y=totalR_est, color=year)) +
  geom_point(size=4) +
  xlab("egg output") + ylab("recruits") +
  #ggtitle("Egg-recruits") +
  theme_bw() +
  theme(text = element_text(size=40)) +
  theme(axis.text.x = element_text(size=30, angle=90, vjust=0.5)) + theme(axis.text.y = element_text(size=30)) +
  theme(legend.text=element_text(size=20)) + theme(legend.title=element_text(size=25))
dev.off()

# Number of eggs produced by site with number of recruits there a year later, years broken out individually
pdf(file = here("Plots/EggRecruitRelationship", "Eggs_recruits_by_site_years_as_subplots.pdf"))
ggplot(data = females_recruits_summary %>% filter(year %in% c(2012,2013,2014,2015,2016,2017)), aes(x=est_eggs, y=totalR_est, shape=site)) +
  geom_point(size=3) +
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) +
  xlab("estimated egg output") + ylab("estimated number of recruits following year") +
  ggtitle("Egg-recruit scatter plot for each site, years as subplots") +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~year) +
  theme(text = element_text(size=25)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + #not sure why this isn't working right now...
  #theme(axis.text.x = element_text(size=20, angle=90, hjust=1)) + theme(axis.text.y = element_text(size=20)) +
  theme_bw()
dev.off()

# Same as above but w/estimated number of females instead of eggs (which looks exactly the same right now b/c eggs aren't dependent on female size yet)
pdf(file = here("Plots/EggRecruitRelationship", "Females_recruits_by_site.pdf"))
ggplot(data = females_recruits_summary %>% filter(year %in% c(2012,2013,2014,2015,2016,2017)), aes(x=totalF_est, y=totalR_est, color=year, shape=site)) +
  geom_point(size=3) +
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) +
  scale_fill_brewer(palette = "Set1") + #not sure why this is having no effect...
  xlab("estimated # breeding F") + ylab("estimated number of recruits following year") +
  ggtitle("Female-recruit scatter plot for each site and year individually") +
  theme(text = element_text(size=25)) +
  theme(axis.text.x = element_text(size=20)) + theme(axis.text.y = element_text(size=20)) +
  theme_bw()
dev.off()

# Truncated w/estimated number of females instead of eggs (which looks exactly the same right now b/c eggs aren't dependent on female size yet)
pdf(file = here("Plots/EggRecruitRelationship", "Females_recruits_by_site_truncated.pdf"))
ggplot(data = females_recruits_summary %>% filter(year %in% c(2012,2013,2014,2015,2016,2017)), aes(x=totalF_est, y=totalR_est, color=year, shape=site)) +
  geom_point(size=3) +
  scale_x_continuous(limits=c(0,150)) +
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) +
  scale_fill_brewer(palette = "Set1") + #not sure why this is having no effect...
  xlab("estimated # breeding F") + ylab("estimated number of recruits following year") +
  ggtitle("Truncated female-recruit scatter plot for each site and year individually") +
  theme(text = element_text(size=25)) +
  theme(axis.text.x = element_text(size=20)) + theme(axis.text.y = element_text(size=20)) +
  theme_bw()
dev.off()

# Raw females and raw recruits (not scaled up by catchability or proportion habitat sampled)
pdf(file = here("Plots/EggRecruitRelationship", "Females_recruits_by_site_rawcounts.pdf"))
ggplot(data = females_recruits_summary %>% filter(year %in% c(2012,2013,2014,2015,2016,2017)), aes(x=NbreedingF_combo, y=Nrecruits, color=year, shape=site)) +
  geom_point(size=3) +
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) +
  scale_fill_brewer(palette = "Set1") + #not sure why this is having no effect...
  xlab("raw # F") + ylab("raw # recruits following year") +
  ggtitle("Raw female-recruit scatter plot for each site and year individually") +
  theme(text = element_text(size=25)) +
  theme(axis.text.x = element_text(size=20)) + theme(axis.text.y = element_text(size=20)) +
  theme_bw()
dev.off()

# Now sum it up by year across all sites, plot females and recruits at metapop level
pdf(file = here("Plots/EggRecruitRelationship", "Metapop_females_recruits_estimated.pdf"))
ggplot(data = metapop_level %>% filter(year %in% c(2012,2013,2014,2015,2016,2017)), aes(x=metapop_estF, y=metapop_estR, color=year)) +
  geom_point(size=3) +
  xlab("estimated # F") + ylab("estimated # R following year") +
  ggtitle("Estimated F and R summed across all sites sampled in each year") +
  theme(text = element_text(size=25)) +
  theme(axis.text.x = element_text(size=20)) + theme(axis.text.y = element_text(size=20)) +
  theme_bw()
dev.off()

# Now sum it up by year across all sites, plot est eggs and recruits at metapop level
pdf(file = here("Plots/EggRecruitRelationship", "Metapop_eggs_recruits_estimated.pdf"))
ggplot(data = metapop_level_mod1 %>% filter(year %in% c(2012,2013,2014,2015,2016,2017)), aes(x=metapop_est_eggs, y=metapop_estR, color=year)) +
  geom_point(size=3) +
  geom_line(color = "black", data=metapop_level_mod1, aes(x=metapop_est_eggs, y=recruits_pred)) +
  xlab("estimated # eggs") + ylab("estimated # R following year") +
  ggtitle("Estimated eggs and R summed across all sites sampled in each year") +
  theme(text = element_text(size=25)) +
  theme(axis.text.x = element_text(size=20)) + theme(axis.text.y = element_text(size=20)) +
  theme_bw()
dev.off()

# Now sum it up by year across all sites, plot raw females and recruits at metapop level
pdf(file = here("Plots/EggRecruitRelationship", "Metapop_females_recruits_raw.pdf"))
ggplot(data = metapop_level %>% filter(year %in% c(2012,2013,2014,2015,2016,2017)), aes(x=metapop_rawF, y=metapop_rawR, color=year)) +
  geom_point(size=3) +
  xlab("raw # F") + ylab("raw # R following year") +
  ggtitle("Raw F and R summed across all sites sampled in each year") +
  theme(text = element_text(size=25)) +
  theme(axis.text.x = element_text(size=20)) + theme(axis.text.y = element_text(size=20)) +
  theme_bw()
dev.off()


#################### Saving output: ####################
save(recruits_info, file=here("Data", "recruits_info.RData"))
save(females_recruits_summary, file=here("Data", "females_recruits_summary.RData"))
save(metapop_level, file=here("Data", "metapop_level.RData"))
save(meta_eggs_recruits_est_mod1, file=here("Data","metapop_lm.RData"))
save(egg_recruits_est_mod1, file=here("Data","egg_recruit_lm.RData"))

#Thinking through the plan
#For each year, try to estimate the number of eggs released (either at a site level or at a whole-population level)
#To do that, count up number of breeding females, scale by proportion habitat sampled, scale by eggs/year (by size if want to)
#For recruits, assess number of fish between 3.5-6.0 (or whatever number it was that Michelle estimated they get to in one year), scale by hab sampled, catchability
#Then plot, eggs for previous year against recruits for this year


