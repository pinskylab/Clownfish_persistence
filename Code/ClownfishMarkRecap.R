# Takes in encounter histories and distances to capture anemone for each marked fish and runs mark-recap analysis to estimate survival and recapture probability

#################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))

# Load libraries
library(RMark)
library(ggplot2)
library(cowplot)

# Read in encounter history data frames for marked fish with sizes (from Clownfish_encounters.R) and distances to capture anemone (from AnemDistFromDiveTrack.R)
load(file = here::here("Data/Script_outputs", "encounters_dist_mean.RData"))  # encounter histories with overall mean dist filled in for NAs (pre-capture years, plus about 5 missing coords)
load(file = here::here("Data/Script_outputs", "encounters_dist_mean_by_year.RData"))  # encounter histories with mean dist by year filled in for NAs
#load(file = here::here("Data/Script_outputs", "encounters_size_means.RData"))  # encounter histories with sizes and missing ones filled in with projections or means
load(file = here::here("Data/Script_outputs", "encounters_size_0.RData"))  # encounter histories with sizes and missing ones filled in with projections or 0s
load(file = here::here("Data/Script_outputs", "encounters_size_means_by_year.RData"))

# Size vector for plotting
n_size_steps <- 30
size_values <- min_size+(0:n_size_steps)*(max_size-min_size)/n_size_steps  # min and max size are in Constants_database_common_functions.R

# Distance vector for plotting
n_dist_steps <- 100
min_dist = 0
max_dist = 200
# max_dist = max(encounters_dist_mean$dist2012, encounters_dist_mean$dist2013, encounters_dist_mean$dist2014, encounters_dist_mean$dist2015,
#                encounters_dist_mean$dist2016, encounters_dist_mean$dist2017, encounters_dist_mean$dist2018, na.rm=TRUE)  # this is very large (766, recap prob goes essentially to zero at < 100)
dist_values = min_dist+(0:n_dist_steps)*(max_dist-min_dist)/n_dist_steps

# Start year
start_year = 2012

# Sites we didn't revisit
one_time_sites <- c("Sitio Lonas", "Sitio Tugas", "Caridad Proper")

#################### Functions: ####################

#################### Running things: ####################
##### Take out fish at sites we didn't re-visit
encounters_dist_mean <- encounters_dist_mean %>% filter(!site %in% one_time_sites)
encounters_dist_mean_by_year <- encounters_dist_mean_by_year %>% filter(!site %in% one_time_sites)
encounters_size_0 <- encounters_size_0 %>% filter(!site %in% one_time_sites)
#encounters_size_means <- encounters_size_means %>% filter(!site %in% one_time_sites)
encounters_size_means_by_year <- encounters_size_means_by_year %>% filter(!site %in% one_time_sites)

##### Create data frames for analysis, joining trait info and distance info
# Mean size (by year), mean distance (by year)
eall_meanYsize_meanYdist <- left_join(encounters_size_means_by_year, encounters_dist_mean_by_year %>% select(fish_indiv, dist2013, dist2014, dist2015, dist2016, dist2017, dist2018),
                                      by = "fish_indiv") %>%
  mutate(site = as.factor(site), cap_color = as.factor(cap_color), cap_stage = as.factor(cap_stage)) %>%
  dplyr::rename(ch = encounter.hist)
  #mutate(site = as.factor(site), capture_color = as.factor(cap_color), capture_stage = as.factor(capture_stage)) %>%
  #dplyr::rename(ch = encounter.hist, cap_color = capture_color, cap_stage = capture_stage)

eall_meanYsize_meanYdist <- eall_meanYsize_meanYdist[complete.cases(eall_meanYsize_meanYdist),]

eall_meanYsize_meanYdist <- eall_meanYsize_meanYdist %>%
  mutate(site = case_when(site == "Cabatoan" ~ "Cabatoan",
                          site == "Caridad Cemetery" ~ "CaridadCemetery",
                          #site == "Caridad Proper" ~ "CaridadProper",
                          site == "Elementary School" ~ "ElementarySchool",
                          site == "Gabas" ~ "Gabas",
                          site == "Haina" ~"Haina",
                          site == "Hicgop South" ~ "HicgopSouth",
                          site == "N. Magbangon" ~ "NMagbangon",
                          site == "Palanas" ~ "Palanas",
                          site == "Poroc Rose" ~ "PorocRose",
                          site == "Poroc San Flower" ~ "PorocSanFlower",
                          site == "S. Magbangon" ~ "SMagbangon",
                          site == "San Agustin" ~ "SanAgustin",
                          site == "Sitio Baybayon" ~ "SitioBaybayon",
                          #site == "Sitio Lonas" ~ "SitioLonas",
                          #site == "Sitio Tugas" ~ "SitioTugas",
                          site == "Tamakin Dacot" ~ "TamakinDacot",
                          site == "Visca" ~ "Visca",
                          site == "Wangag" ~ "Wangag"))

##### Prep for MARK runs
# Process data (both with and without site as a group)
eall_meanYsize_meanYdist.processed = process.data(eall_meanYsize_meanYdist, model='CJS', begin.time=start_year)
eall_meanYsize_meanYdist.processed_site = process.data(eall_meanYsize_meanYdist, model="CJS", begin.time=start_year, groups="site")

# Make ddl
eall_meanYsize_meanYdist.ddl = make.design.data(eall_meanYsize_meanYdist.processed)
eall_meanYsize_meanYdist.ddl_site = make.design.data(eall_meanYsize_meanYdist.processed_site)

# Set models for Phi
Phi.dot = list(formula=~1, link="logit")
Phi.time = list(formula=~time, link="logit")
Phi.site = list(formula=~site, link="logit")
Phi.size = list(formula=~size, link="logit")
Phi.site.plus.size = list(formula=~site+size, link="logit")

# Set models for p
p.dot = list(formula=~1, link="logit")
p.time = list(formula=~time, link="logit")
p.dist = list(formula=~dist, link="logit")
p.site = list(formula=~site, link="logit")
p.size = list(formula=~size, link="logit")
p.size.plus.dist = list(formula=~size+dist, link="logit")

# Run models
# using mean-dist-for-NA and mean-size-for-NA (mean within each year) dataset
eall_meanYsize_meanYdist.Phi.dot.p.dot = mark(eall_meanYsize_meanYdist.processed, eall_meanYsize_meanYdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot), prefix = "eall_meanYsize_meanYdist.Phi.dot.p.dot")  # constant survival and recap, no covariates
eall_meanYsize_meanYdist.Phi.dot.p.time = mark(eall_meanYsize_meanYdist.processed, eall_meanYsize_meanYdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.time), prefix = "eall_meanYsize_meanYdist.Phi.dot.p.time")  # constant survival, recapture time-dependent
eall_meanYsize_meanYdist.Phi.time.p.dot = mark(eall_meanYsize_meanYdist.processed, eall_meanYsize_meanYdist.ddl, model.parameters=list(Phi=Phi.time, p=p.dot), prefix = "eall_meanYsize_meanYdist.Phi.time.p.dot")  # time-varying survival, recapture constant
eall_meanYsize_meanYdist.Phi.dot.p.dist = mark(eall_meanYsize_meanYdist.processed, eall_meanYsize_meanYdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist), prefix = "eall_meanYsize_meanYdist.Phi.dot.p.dist")  # constant survival, recapture distance-dependent
eall_meanYsize_meanYdist.Phi.size.p.size = mark(eall_meanYsize_meanYdist.processed, eall_meanYsize_meanYdist.ddl, model.parameters=list(Phi=Phi.size, p=p.size), prefix = "eall_meanYsize_meanYdist.Phi.size.p.size")  # capture-size-dependent survival and recapture
eall_meanYsize_meanYdist.Phi.size.p.dist = mark(eall_meanYsize_meanYdist.processed, eall_meanYsize_meanYdist.ddl, model.parameters=list(Phi=Phi.size, p=p.dist), prefix = "eall_meanYsize_meanYdist.Phi.size.p.dist")  # size-dependent survival, distance-dependent recapture
eall_meanYsize_meanYdist.Phi.size.p.size.plus.dist = mark(eall_meanYsize_meanYdist.processed, eall_meanYsize_meanYdist.ddl, model.parameters=list(Phi=Phi.size, p=p.size.plus.dist), prefix = "eall_meanYsize_meanYdist.Phi.size.p.size.plus.dist")  # size-dependent survival, size-and-distance-dependent recapture
eall_meanYsize_meanYdist.Phi.dot.p.size.plus.dist = mark(eall_meanYsize_meanYdist.processed, eall_meanYsize_meanYdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.size.plus.dist), prefix = "eall_meanYsize_meanYdist.Phi.dot.p.size.plus.dist")  # constant survival, size-and-distance-dependent recapture
eall_meanYsize_meanYdist.Phi.site.p.size.plus.dist = mark(eall_meanYsize_meanYdist.processed_site, eall_meanYsize_meanYdist.ddl_site, model.parameters=list(Phi=Phi.site, p=p.size.plus.dist), prefix = "eall_meanYsize_meanYdist.Phi.site.p.size.plus.dist")  # site-dependent survival, size-and-distance-dependent recapture
eall_meanYsize_meanYdist.Phi.site.p.dot = mark(eall_meanYsize_meanYdist.processed_site, eall_meanYsize_meanYdist.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dot), prefix = "eall_meanYsize_meanYdist.Phi.site.p.dot")  # site-dependent survival, constant recapture
eall_meanYsize_meanYdist.Phi.site.p.dist = mark(eall_meanYsize_meanYdist.processed_site, eall_meanYsize_meanYdist.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dist), prefix = "eall_meanYsize_meanYdist.Phi.site.p.dist")  # site-dependent survival, dist-related recapture
eall_meanYsize_meanYdist.Phi.site.plus.size.p.size.plus.dist = mark(eall_meanYsize_meanYdist.processed_site, eall_meanYsize_meanYdist.ddl_site, model.parameters=list(Phi=Phi.site.plus.size, p=p.size.plus.dist), prefix = "eall_meanYdist_meanYsize.Phi.site.plus.size.p.size.plus.dist")  # site- and size-dependent survival, size-and-distance-dependent recapture

# Compare model AICc
model_comp_meanYsize_meanYdist = data.frame(model = c('eall_meanYsize_meanYdist.Phi.dot.p.dot','eall_meanYsize_meanYdist.Phi.dot.p.time','eall_meanYsize_meanYdist.Phi.time.p.dot',
                                       'eall_meanYsize_meanYdist.Phi.dot.p.dist','eall_meanYsize_meanYdist.Phi.size.p.size','eall_meanYsize_meanYdist.Phi.size.p.dist',
                                       'eall_meanYsize_meanYdist.Phi.size.p.size.plus.dist','eall_meanYsize_meanYdist.Phi.dot.p.size.plus.dist',
                                       'eall_meanYsize_meanYdist.Phi.site.p.dot','eall_meanYsize_meanYdist.Phi.site.p.size.plus.dist', 
                                       'eall_meanYsize_meanYdist.Phi.site.p.dist','eall_meanYsize_meanYdist.Phi.site.plus.size.p.size.plus.dist'),
                             AICc = c(eall_meanYsize_meanYdist.Phi.dot.p.dot$results$AICc, eall_meanYsize_meanYdist.Phi.dot.p.time$results$AICc, eall_meanYsize_meanYdist.Phi.time.p.dot$results$AICc,
                                      eall_meanYsize_meanYdist.Phi.dot.p.dist$results$AICc, eall_meanYsize_meanYdist.Phi.size.p.size$results$AICc, eall_meanYsize_meanYdist.Phi.size.p.dist$results$AICc,
                                      eall_meanYsize_meanYdist.Phi.size.p.size.plus.dist$results$AICc, eall_meanYsize_meanYdist.Phi.dot.p.size.plus.dist$results$AICc,
                                      eall_meanYsize_meanYdist.Phi.site.p.dot$results$AICc, eall_meanYsize_meanYdist.Phi.site.p.size.plus.dist$results$AICc,
                                      eall_meanYsize_meanYdist.Phi.site.p.dist$results$AICc, eall_meanYsize_meanYdist.Phi.site.plus.size.p.size.plus.dist$results$AICc),
                             stringsAsFactors = FALSE)

# Arrange by increasing AICc
model_comp_meanYsize_meanYdist <- model_comp_meanYsize_meanYdist %>%
  mutate(dAICc = min(AICc) - AICc) %>%
  arrange(-dAICc)

#################### Data frames for plotting for best-fit model: ####################
###### survival by site and size, recapture by both size and distance (eall_meanYsize_meanYdist.Phi.site.plus.size.p.size.plus.dist)
# Pull out results to save
results_df_1 <- as.data.frame(eall_meanYsize_meanYdist.Phi.site.plus.size.p.size.plus.dist$results$beta)

# Set indices for effects within results data frame
surv_intercept_index = 1
surv_size_effect_index = 17
recap_intercept_index = 18
recap_size_effect_index = 19
recap_dist_effect_index = 20

### Look at surv by site for mean size (overall, not by site)
mean_size_MARK <- mean(c(eall_meanYsize_meanYdist$size2012, eall_meanYsize_meanYdist$size2013, eall_meanYsize_meanYdist$size2014,
                         eall_meanYsize_meanYdist$size2015, eall_meanYsize_meanYdist$size2016, eall_meanYsize_meanYdist$size2017, eall_meanYsize_meanYdist$size2018))
# start with Cabatoan (intercept is for Cabatoan)
surv_by_site_mean_size <- data.frame(site = no_space_sites_revisited[1], size = mean_size_MARK) %>%
  mutate(estimate = results_df_1$estimate[1] + results_df_1$estimate[surv_size_effect_index]*size,
         lcl = results_df_1$lcl[1] + results_df_1$lcl[surv_size_effect_index]*size,
         ucl = results_df_1$ucl[1] + results_df_1$ucl[surv_size_effect_index]*size,
         estimate_prob = logit_recip(estimate),
         lcl_prob = logit_recip(lcl),
         ucl_prob = logit_recip(ucl))
# fill in rest of sites
for(i in 2:length(no_space_sites_revisited)) {
  df_to_add <- data.frame(site = no_space_sites_revisited[i], size = mean_size_MARK) %>%
    mutate(estimate = results_df_1$estimate[1] + results_df_1$estimate[i] + results_df_1$estimate[surv_size_effect_index]*size,
           lcl = results_df_1$lcl[1] + results_df_1$lcl[i] + results_df_1$lcl[surv_size_effect_index]*size,
           ucl = results_df_1$ucl[1] + results_df_1$ucl[i] + results_df_1$ucl[surv_size_effect_index]*size,
           estimate_prob = logit_recip(estimate),
           lcl_prob = logit_recip(lcl),
           ucl_prob = logit_recip(ucl))
  
  surv_by_site_mean_size <- rbind(surv_by_site_mean_size, df_to_add)
}

### Look at surv by site for mean size by each site
mean_size_MARK_site <- eall_meanYsize_meanYdist %>%
  group_by(site) %>%
  summarize(mean_size = mean(c(size2012, size2013, size2014, size2015, size2016, size2017, size2018)))
                  
# start with Cabatoan
surv_by_site_mean_size_by_site <- data.frame(site = no_space_sites_revisited[1], size = mean_size_MARK_site$mean_size[1]) %>%
  mutate(estimate = results_df_1$estimate[1] + results_df_1$estimate[surv_size_effect_index]*size,
         lcl = results_df_1$lcl[1] + results_df_1$lcl[surv_size_effect_index]*size,
         ucl = results_df_1$ucl[1] + results_df_1$ucl[surv_size_effect_index]*size,
         estimate_prob = logit_recip(estimate),
         lcl_prob = logit_recip(lcl),
         ucl_prob = logit_recip(ucl))
# fill in rest of sites
for(i in 2:length(no_space_sites_revisited)) {
  df_to_add <- data.frame(site = no_space_sites_revisited[i], size = mean_size_MARK_site$mean_size[i]) %>%
    mutate(estimate = results_df_1$estimate[1] + results_df_1$estimate[i] + results_df_1$estimate[surv_size_effect_index]*size,
           lcl = results_df_1$lcl[1] + results_df_1$lcl[i] + results_df_1$lcl[surv_size_effect_index]*size,
           ucl = results_df_1$ucl[1] + results_df_1$ucl[i] + results_df_1$ucl[surv_size_effect_index]*size,
           estimate_prob = logit_recip(estimate),
           lcl_prob = logit_recip(lcl),
           ucl_prob = logit_recip(ucl))
  
  surv_by_site_mean_size_by_site <- rbind(surv_by_site_mean_size_by_site, df_to_add)
}

### Now find relationship with size by site
# start with Cabatoan
df_size_by_site <- data.frame(site = no_space_sites_revisited[1], size = size_values) %>%
  mutate(estimate = results_df_1$estimate[1] + results_df_1$estimate[surv_size_effect_index]*size,
         lcl = results_df_1$lcl[1] + results_df_1$lcl[surv_size_effect_index]*size,
         ucl = results_df_1$ucl[1] + results_df_1$ucl[surv_size_effect_index]*size,
         estimate_prob = logit_recip(estimate),
         lcl_prob = logit_recip(lcl),
         ucl_prob = logit_recip(ucl))
# fill in rest 
for(i in 2:length(no_space_sites_revisited)) {
  df_to_add <- data.frame(site = no_space_sites_revisited[i], size = size_values) %>%
    mutate(estimate = results_df_1$estimate[1] + results_df_1$estimate[i] + results_df_1$estimate[surv_size_effect_index]*size,
           lcl = results_df_1$lcl[1] + results_df_1$lcl[i] + results_df_1$lcl[surv_size_effect_index]*size,
           ucl = results_df_1$ucl[1] + results_df_1$ucl[i] + results_df_1$ucl[surv_size_effect_index]*size,
           estimate_prob = logit_recip(estimate),
           lcl_prob = logit_recip(lcl),
           ucl_prob = logit_recip(ucl))
  
  df_size_by_site <- rbind(df_size_by_site, df_to_add)
}

### Find size effect on recap
recap_by_size <- data.frame(size = size_values) %>%
  mutate(estimate = results_df_1$estimate[recap_intercept_index] + results_df_1$estimate[recap_size_effect_index]*size,
         lcl = results_df_1$lcl[recap_intercept_index] + results_df_1$lcl[recap_size_effect_index]*size,
         ucl = results_df_1$ucl[recap_intercept_index] + results_df_1$ucl[recap_size_effect_index]*size,
         estimate_prob = logit_recip(estimate),
         lcl_prob = logit_recip(lcl),
         ucl_prob = logit_recip(ucl))

### Find distance effect on recap
recap_by_dist <- data.frame(dist = dist_values) %>%
  mutate(estimate = results_df_1$estimate[recap_intercept_index] + results_df_1$estimate[recap_dist_effect_index]*dist,
         lcl = results_df_1$lcl[recap_intercept_index] + results_df_1$lcl[recap_dist_effect_index]*dist,
         ucl = results_df_1$ucl[recap_intercept_index] + results_df_1$ucl[recap_dist_effect_index]*dist,
         estimate_prob = logit_recip(estimate),
         lcl_prob = logit_recip(lcl),
         ucl_prob = logit_recip(ucl))

# Put best-fit model and data frames together into a list to save
best_fit_model_dfs <- list("model" = eall_meanYsize_meanYdist.Phi.site.plus.size.p.size.plus.dist, "results" = results_df_1, 
                           "surv_site_size" = df_size_by_site, "recap_size" = recap_by_size, "recap_dist" = recap_by_dist)


#################### Pull out model where site isn't in survival, in case want to use: ####################
# Phi.size.p.size.plus.dist, 4th best-fit

# Pull out results to save
results_df_4 <- as.data.frame(eall_meanYsize_meanYdist.Phi.size.p.size.plus.dist$results$beta)

best_fit_no_site_in_surv_model_dfs <- list("model" = eall_meanYsize_meanYdist.Phi.size.p.size.plus.dist, "results" = results_df_4)

#################### Check a few other models to see how they compare: ####################
###### 2nd best: survival by site, recapture by both size and distance (eall_meanYsize_meanYdist.Phi.site.p.size.plus.dist)
# Mean survival by site
results_df_2 <- as.data.frame(eall_meanYsize_meanYdist.Phi.site.p.size.plus.dist$results$beta)

site_specific_surv <- data.frame(site = no_space_sites_revisited, estimate = NA, se = NA, lcl = NA, ucl = NA)
# Cabatoan is the first
site_specific_surv$estimate[1] = results_df_2$estimate[1]
site_specific_surv$lcl[1] = results_df_2$lcl[1]
site_specific_surv$ucl[1] = results_df_2$ucl[1]

# Go through rest of sites - need to be added to Cabatoan estimate
for (i in 2:length(no_space_sites_revisited)) {
  site_specific_surv$estimate[i] = results_df_2$estimate[1] + results_df_2$estimate[i]
  site_specific_surv$lcl[i] = results_df_2$lcl[1] + results_df_2$lcl[i]
  site_specific_surv$ucl[i] = results_df_2$ucl[1] + results_df_2$ucl[i]
}

# Convert back from logit scale
site_specific_surv <- site_specific_surv %>%
  mutate(estimate_prob = logit_recip(estimate),
         lcl_prob = logit_recip(lcl),
         ucl_prob = logit_recip(ucl))

###### 3rd best: survival by site, recapture by distance (eall_meanYsize_meanYdist.Phi.site.p.dist)
# Mean survival by site
results_df_3 <- as.data.frame(eall_meanYsize_meanYdist.Phi.site.p.dist$results$beta)

site_specific_surv_nodist <- data.frame(site = no_space_sites_revisited, estimate = NA, se = NA, lcl = NA, ucl = NA)

# Cabatoan is the first
site_specific_surv_nodist$estimate[1] = results_df_3$estimate[1]
site_specific_surv_nodist$lcl[1] = results_df_3$lcl[1]
site_specific_surv_nodist$ucl[1] = results_df_3$ucl[1]

# Go through rest of sites - need to be added to Cabatoan estimate
for (i in 2:length(no_space_sites_revisited)) {
  site_specific_surv_nodist$estimate[i] = results_df_3$estimate[1] + results_df_3$estimate[i]
  site_specific_surv_nodist$lcl[i] = results_df_3$lcl[1] + results_df_3$lcl[i]
  site_specific_surv_nodist$ucl[i] = results_df_3$ucl[1] + results_df_3$ucl[i]
}

# Convert back from logit scale
site_specific_surv_nodist <- site_specific_surv_nodist %>%
  mutate(estimate_prob = logit_recip(estimate),
         lcl_prob = logit_recip(lcl),
         ucl_prob = logit_recip(ucl))


#################### Plots: ####################
##### Best fit: Phi by site+size, recap by size+dist (eall_meanYsize_meanYdist.Phi.site.plus.size.p.size.plus.dist)
# Phi (by site), mean size overall 
pdf(file = here::here("Plots/PhiandpEstimates", "surv_by_site_mean_size_Phisiteplussize_psizeplusdist.pdf"))
ggplot(data = surv_by_site_mean_size, aes(x = site, y = estimate_prob, color = site, fill = site)) +
  geom_point() + 
  geom_linerange(aes(ymin=lcl_prob, ymax=ucl_prob)) +
  ylab("survival probability") + ggtitle("Phi:site+size, p:size+dist, mean size overall") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Phi (by site), mean size by site
pdf(file = here::here("Plots/PhiandpEstimates", "surv_by_site_mean_size_by_site_Phisiteplussize_psizeplusdist.pdf"))
ggplot(data = surv_by_site_mean_size_by_site, aes(x = site, y = estimate_prob, color = site, fill = site)) +
  geom_point() + 
  geom_linerange(aes(ymin=lcl_prob, ymax=ucl_prob)) +
  ylab("survival probability") + ggtitle("Phi:site+size, p:size+dist, mean size by site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Phi (by size and site)
pdf(file = here::here("Plots/PhiandpEstimates", "surv_by_size_and_site_Phisiteplussize_psizeplusdist.pdf"))
ggplot(data = df_size_by_site, aes(x = size, y = estimate_prob, fill = site)) +
  geom_line(aes(color=site)) + 
  geom_ribbon(aes(ymin=lcl_prob, ymax=ucl_prob), alpha=0.5) +
  ylab("survival probability") + ggtitle("Phi:site+size, p:size+dist") + xlab("size (cm)") +
  facet_wrap(~site) +
  theme_bw() 
dev.off()

# p (by dist)
pdf(file = here::here("Plots/PhiandpEstimates", "recap_by_dist_Phisiteplussize_psizeplusdist.pdf"))
ggplot(data = recap_by_dist, aes(x = dist, y = estimate_prob)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl_prob, ymax=ucl_prob), alpha=0.5) +
  ylab("recapture probability") + ggtitle("Phi:site+size, p:size+dist") + xlab("distance (m)") +
  theme_bw() 
dev.off()

# p (by size)
pdf(file = here::here("Plots/PhiandpEstimates", "recap_by_size_Phisiteplussize_psizeplusdist.pdf"))
ggplot(data = recap_by_size, aes(x = size, y = estimate_prob)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl_prob, ymax=ucl_prob), alpha=0.5) +
  ylab("recapture probability") + ggtitle("Phi:site+size, p:size+dist") + xlab("size (cm)") +
  theme_bw() 
dev.off()

# Put recap effects together
recap_size_plot <- ggplot(data = recap_by_size, aes(x = size, y = estimate_prob)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl_prob, ymax=ucl_prob), alpha=0.5) +
  ylab("recapture probability") + ggtitle("Phi:site+size, p:size+dist") + xlab("size (cm)") +
  theme_bw() 

recap_dist_plot <- ggplot(data = recap_by_dist, aes(x = dist, y = estimate_prob)) +
  geom_line() + 
  geom_ribbon(aes(ymin=lcl_prob, ymax=ucl_prob), alpha=0.5) +
  ylab("recapture probability") + ggtitle("Phi:site+size, p:size+dist") + xlab("distance (m)") +
  theme_bw() 

pdf(file = here::here("Plots/PhiandpEstimates", "recap_effects_Phisiteplussize_psizeplusdist.pdf"))
plot_grid(recap_dist_plot, recap_size_plot, labels=c("a","b"), nrow=1)
dev.off()


##### Other models (second best fit), just to compare: Phi by site, recap by size+dist (eall_meanYsize_meanYdist.Phi.site.p.size.plus.dist)
# Phi (by site)
pdf(file = here::here("Plots/PhiandpEstimates", "surv_by_site_recap_by_sizeplusdist.pdf"))
ggplot(data = site_specific_surv, aes(x = site, y = estimate_prob, color = site, fill = site)) +
  geom_point() + 
  geom_linerange(aes(ymin=lcl_prob, ymax=ucl_prob)) +
  ylab("survival probability") + ggtitle("Phi:site, p:size+dist") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

##### 3rd best fit: Phi by site, recap by dist (eall_meanYsize_meanYdist.Phi.site.p.dist)
# Phi (by site)
pdf(file = here::here("Plots/PhiandpEstimates", "surv_by_site_recap_by_dist.pdf"))
ggplot(data = site_specific_surv_nodist, aes(x = site, y = estimate_prob, color = site, fill = site)) +
  geom_point() + 
  geom_linerange(aes(ymin=lcl_prob, ymax=ucl_prob)) +
  ylab("survival probability") + ggtitle("Phi:site, p:dist") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#################### Saving output: ####################
save(model_comp_meanYsize_meanYdist, file=here::here("Data/Script_outputs", "model_comp_meanYsize_meanYdist.RData"))
save(best_fit_model_dfs, file=here::here("Data/Script_outputs", "best_fit_model_dfs.RData"))
save(best_fit_no_site_in_surv_model_dfs, file=here::here("Data/Script_outputs", "best_fit_no_site_in_surv_model_dfs.RData"))

#################### Sensitivity to how fill in missing values tests: ####################
##### Do some sensitivity to test that the same models get picked when fill in missing values in other ways
### Mean size (by year), mean distance (overall)
# Set data frames
eall_meanYsize_meanOdist <- left_join(encounters_size_means_by_year, encounters_dist_mean %>% select(fish_indiv, dist2013, dist2014, dist2015, dist2016, dist2017, dist2018),
                                      by = "fish_indiv") %>%
  mutate(site = as.factor(site), cap_color = as.factor(cap_color), cap_stage = as.factor(cap_stage)) %>%
  dplyr::rename(ch = encounter.hist)

  #mutate(site = as.factor(site), capture_color = as.factor(capture_color), capture_stage = as.factor(capture_stage)) %>%
  #dplyr::rename(ch = encounter.hist, cap_color = capture_color, cap_stage = capture_stage)

eall_meanYsize_meanOdist <- eall_meanYsize_meanOdist[complete.cases(eall_meanYsize_meanOdist),]

eall_meanYsize_meanOdist <- eall_meanYsize_meanOdist %>%
  mutate(site = case_when(site == "Cabatoan" ~ "Cabatoan",
                          site == "Caridad Cemetery" ~ "CaridadCemetery",
                          site == "Caridad Proper" ~ "CaridadProper",
                          site == "Elementary School" ~ "ElementarySchool",
                          site == "Gabas" ~ "Gabas",
                          site == "Haina" ~"Haina",
                          site == "Hicgop South" ~ "HicgopSouth",
                          site == "N. Magbangon" ~ "NMagbangon",
                          site == "Palanas" ~ "Palanas",
                          site == "Poroc Rose" ~ "PorocRose",
                          site == "Poroc San Flower" ~ "PorocSanFlower",
                          site == "S. Magbangon" ~ "SMagbangon",
                          site == "San Agustin" ~ "SanAgustin",
                          site == "Sitio Baybayon" ~ "SitioBaybayon",
                          site == "Sitio Lonas" ~ "SitioLonas",
                          site == "Sitio Tugas" ~ "SitioTugas",
                          site == "Tamakin Dacot" ~ "TamakinDacot",
                          site == "Visca" ~ "Visca",
                          site == "Wangag" ~ "Wangag"))

# Process data and make ddl
eall_meanYsize_meanOdist.processed = process.data(eall_meanYsize_meanOdist, model='CJS', begin.time=2012)
eall_meanYsize_meanOdist.processed_site = process.data(eall_meanYsize_meanOdist, model="CJS", begin.time=2012, groups="site")

eall_meanYsize_meanOdist.ddl = make.design.data(eall_meanYsize_meanOdist.processed)
eall_meanYsize_meanOdist.ddl_site = make.design.data(eall_meanYsize_meanOdist.processed_site)

# Run models
eall_meanYsize_meanOdist.Phi.dot.p.dot = mark(eall_meanYsize_meanOdist.processed, eall_meanYsize_meanOdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot), prefix = "eall_meanYsize_meanOdist.Phi.dot.p.dot")  # constant survival and recap, no covariates
eall_meanYsize_meanOdist.Phi.dot.p.time = mark(eall_meanYsize_meanOdist.processed, eall_meanYsize_meanOdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.time), prefix = "eall_meanYsize_meanOdist.Phi.dot.p.time")  # constant survival, recapture time-dependent
eall_meanYsize_meanOdist.Phi.time.p.dot = mark(eall_meanYsize_meanOdist.processed, eall_meanYsize_meanOdist.ddl, model.parameters=list(Phi=Phi.time, p=p.dot), prefix = "eall_meanYsize_meanOdist.Phi.time.p.dot")  # time-varying survival, recapture constant
eall_meanYsize_meanOdist.Phi.dot.p.dist = mark(eall_meanYsize_meanOdist.processed, eall_meanYsize_meanOdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist), prefix = "eall_meanYsize_meanOdist.Phi.dot.p.dist")  # constant survival, recapture distance-dependent
eall_meanYsize_meanOdist.Phi.size.p.size = mark(eall_meanYsize_meanOdist.processed, eall_meanYsize_meanOdist.ddl, model.parameters=list(Phi=Phi.size, p=p.size), prefix = "eall_meanYsize_meanOdist.Phi.size.p.size")  # capture-size-dependent survival and recapture
eall_meanYsize_meanOdist.Phi.size.p.dist = mark(eall_meanYsize_meanOdist.processed, eall_meanYsize_meanOdist.ddl, model.parameters=list(Phi=Phi.size, p=p.dist), prefix = "eall_meanYsize_meanOdist.Phi.size.p.dist")  # size-dependent survival, distance-dependent recapture
eall_meanYsize_meanOdist.Phi.size.p.size.plus.dist = mark(eall_meanYsize_meanOdist.processed, eall_meanYsize_meanOdist.ddl, model.parameters=list(Phi=Phi.size, p=p.size.plus.dist), prefix = "eall_meanYsize_meanOdist.Phi.size.p.size.plus.dist")  # size-dependent survival, size-and-distance-dependent recapture
eall_meanYsize_meanOdist.Phi.dot.p.size.plus.dist = mark(eall_meanYsize_meanOdist.processed, eall_meanYsize_meanOdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.size.plus.dist), prefix = "eall_meanYsize_meanOdist.Phi.dot.p.size.plus.dist")  # constant survival, size-and-distance-dependent recapture
eall_meanYsize_meanOdist.Phi.site.p.size.plus.dist = mark(eall_meanYsize_meanOdist.processed_site, eall_meanYsize_meanOdist.ddl_site, model.parameters=list(Phi=Phi.site, p=p.size.plus.dist), prefix = "eall_meanYsize_meanOdist.Phi.site.p.size.plus.dist")  # site-dependent survival, size-and-distance-dependent recapture
eall_meanYsize_meanOdist.Phi.site.p.dot = mark(eall_meanYsize_meanOdist.processed_site, eall_meanYsize_meanOdist.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dot), prefix = "eall_meanYsize_meanOdist.Phi.site.p.dot")  # site-dependent survival, constant recapture
eall_meanYsize_meanOdist.Phi.site.p.dist = mark(eall_meanYsize_meanOdist.processed_site, eall_meanYsize_meanOdist.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dist), prefix = "eall_meanYsize_meanOdist.Phi.site.p.dist")  # site-dependent survival, dist-related recapture
eall_meanYsize_meanOdist.Phi.site.plus.size.p.size.plus.dist = mark(eall_meanYsize_meanOdist.processed_site, eall_meanYsize_meanOdist.ddl_site, model.parameters=list(Phi=Phi.site.plus.size, p=p.size.plus.dist), prefix = "eall_meanYdist_meanOsize.Phi.site.plus.size.p.size.plus.dist")  # site- and size-dependent survival, size-and-distance-dependent recapture

# Compare model AICc
model_comp_meanYsize_meanOdist = data.frame(model = c('eall_meanYsize_meanOdist.Phi.dot.p.dot','eall_meanYsize_meanOdist.Phi.dot.p.time','eall_meanYsize_meanOdist.Phi.time.p.dot',
                                                      'eall_meanYsize_meanOdist.Phi.dot.p.dist','eall_meanYsize_meanOdist.Phi.size.p.size','eall_meanYsize_meanOdist.Phi.size.p.dist',
                                                      'eall_meanYsize_meanOdist.Phi.size.p.size.plus.dist','eall_meanYsize_meanOdist.Phi.dot.p.size.plus.dist',
                                                      'eall_meanYsize_meanOdist.Phi.site.p.dot','eall_meanYsize_meanOdist.Phi.site.p.size.plus.dist',
                                                      "eall_meanYsize_meanOdist.Phi.site.p.dist","eall_meanYsize_meanOdist.Phi.site.plus.size.p.size.plus.dist"),
                                            AICc = c(eall_meanYsize_meanOdist.Phi.dot.p.dot$results$AICc, eall_meanYsize_meanOdist.Phi.dot.p.time$results$AICc, eall_meanYsize_meanOdist.Phi.time.p.dot$results$AICc,
                                                     eall_meanYsize_meanOdist.Phi.dot.p.dist$results$AICc, eall_meanYsize_meanOdist.Phi.size.p.size$results$AICc, eall_meanYsize_meanOdist.Phi.size.p.dist$results$AICc,
                                                     eall_meanYsize_meanOdist.Phi.size.p.size.plus.dist$results$AICc, eall_meanYsize_meanOdist.Phi.dot.p.size.plus.dist$results$AICc,
                                                     eall_meanYsize_meanOdist.Phi.site.p.dot$results$AICc, eall_meanYsize_meanOdist.Phi.site.p.size.plus.dist$results$AICc,
                                                     eall_meanYsize_meanOdist.Phi.site.p.dist$results$AICc, eall_meanYsize_meanOdist.Phi.site.plus.size.p.size.plus.dist$results$AICc),
                                            stringsAsFactors = FALSE)

# Arrange by increasing AICc
model_comp_meanYsize_meanOdist <- model_comp_meanYsize_meanOdist %>%
  mutate(dAICc = min(AICc) - AICc) %>%
  arrange(-dAICc)

### 0 size, mean distance (overall)
# Make data frames
eall_0size_meanOdist <- left_join(encounters_size_0, encounters_dist_mean %>% select(fish_indiv, dist2013, dist2014, dist2015, dist2016, dist2017, dist2018),
                                  by = "fish_indiv") %>%
  mutate(site = as.factor(site), cap_color = as.factor(cap_color), cap_stage = as.factor(cap_stage)) %>%
  dplyr::rename(ch = encounter.hist)
           
  #mutate(site = as.factor(site), capture_color = as.factor(capture_color), capture_stage = as.factor(capture_stage)) %>%
  #dplyr::rename(ch = encounter.hist, cap_color = capture_color, cap_stage = capture_stage)

eall_0size_meanOdist <- eall_0size_meanOdist[complete.cases(eall_0size_meanOdist),]

eall_0size_meanOdist <- eall_0size_meanOdist %>%
  mutate(site = case_when(site == "Cabatoan" ~ "Cabatoan",
                          site == "Caridad Cemetery" ~ "CaridadCemetery",
                          site == "Caridad Proper" ~ "CaridadProper",
                          site == "Elementary School" ~ "ElementarySchool",
                          site == "Gabas" ~ "Gabas",
                          site == "Haina" ~"Haina",
                          site == "Hicgop South" ~ "HicgopSouth",
                          site == "N. Magbangon" ~ "NMagbangon",
                          site == "Palanas" ~ "Palanas",
                          site == "Poroc Rose" ~ "PorocRose",
                          site == "Poroc San Flower" ~ "PorocSanFlower",
                          site == "S. Magbangon" ~ "SMagbangon",
                          site == "San Agustin" ~ "SanAgustin",
                          site == "Sitio Baybayon" ~ "SitioBaybayon",
                          site == "Sitio Lonas" ~ "SitioLonas",
                          site == "Sitio Tugas" ~ "SitioTugas",
                          site == "Tamakin Dacot" ~ "TamakinDacot",
                          site == "Visca" ~ "Visca",
                          site == "Wangag" ~ "Wangag"))


# Process data and make ddl
eall_0size_meanOdist.processed = process.data(eall_0size_meanOdist, model='CJS', begin.time=2012)
eall_0size_meanOdist.processed_site = process.data(eall_0size_meanOdist, model="CJS", begin.time=2012, groups="site")

eall_0size_meanOdist.ddl = make.design.data(eall_0size_meanOdist.processed)
eall_0size_meanOdist.ddl_site = make.design.data(eall_0size_meanOdist.processed_site)

# Run models
eall_0size_meanOdist.Phi.dot.p.dot = mark(eall_0size_meanOdist.processed, eall_0size_meanOdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot), prefix = "eall_0size_meanOdist.Phi.dot.p.dot")  # constant survival and recap, no covariates
eall_0size_meanOdist.Phi.dot.p.time = mark(eall_0size_meanOdist.processed, eall_0size_meanOdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.time), prefix = "eall_0size_meanOdist.Phi.dot.p.time")  # constant survival, recapture time-dependent
eall_0size_meanOdist.Phi.time.p.dot = mark(eall_0size_meanOdist.processed, eall_0size_meanOdist.ddl, model.parameters=list(Phi=Phi.time, p=p.dot), prefix = "eall_0size_meanOdist.Phi.time.p.dot")  # time-varying survival, recapture constant
eall_0size_meanOdist.Phi.dot.p.dist = mark(eall_0size_meanOdist.processed, eall_0size_meanOdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist), prefix = "eall_0size_meanOdist.Phi.dot.p.dist")  # constant survival, recapture distance-dependent
eall_0size_meanOdist.Phi.size.p.size = mark(eall_0size_meanOdist.processed, eall_0size_meanOdist.ddl, model.parameters=list(Phi=Phi.size, p=p.size), prefix = "eall_0size_meanOdist.Phi.size.p.size")  # capture-size-dependent survival and recapture
eall_0size_meanOdist.Phi.size.p.dist = mark(eall_0size_meanOdist.processed, eall_0size_meanOdist.ddl, model.parameters=list(Phi=Phi.size, p=p.dist), prefix = "eall_0size_meanOdist.Phi.size.p.dist")  # size-dependent survival, distance-dependent recapture
eall_0size_meanOdist.Phi.size.p.size.plus.dist = mark(eall_0size_meanOdist.processed, eall_0size_meanOdist.ddl, model.parameters=list(Phi=Phi.size, p=p.size.plus.dist), prefix = "eall_0size_meanOdist.Phi.size.p.size.plus.dist")  # size-dependent survival, size-and-distance-dependent recapture
eall_0size_meanOdist.Phi.dot.p.size.plus.dist = mark(eall_0size_meanOdist.processed, eall_0size_meanOdist.ddl, model.parameters=list(Phi=Phi.dot, p=p.size.plus.dist), prefix = "eall_0size_meanOdist.Phi.dot.p.size.plus.dist")  # constant survival, size-and-distance-dependent recapture
eall_0size_meanOdist.Phi.site.p.size.plus.dist = mark(eall_0size_meanOdist.processed_site, eall_0size_meanOdist.ddl_site, model.parameters=list(Phi=Phi.site, p=p.size.plus.dist), prefix = "eall_0size_meanOdist.Phi.site.p.size.plus.dist")  # site-dependent survival, size-and-distance-dependent recapture
eall_0size_meanOdist.Phi.site.p.dot = mark(eall_0size_meanOdist.processed_site, eall_0size_meanOdist.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dot), prefix = "eall_0size_meanOdist.Phi.site.p.dot")  # site-dependent survival, constant recapture
eall_0size_meanOdist.Phi.site.p.dist = mark(eall_0size_meanOdist.processed_site, eall_0size_meanOdist.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dist), prefix = "eall_0size_meanOdist.Phi.site.p.dist")  # site-dependent survival, dist-related recapture
eall_0size_meanOdist.Phi.site.plus.size.p.size.plus.dist = mark(eall_0size_meanOdist.processed_site, eall_0size_meanOdist.ddl_site, model.parameters=list(Phi=Phi.site.plus.size, p=p.size.plus.dist), prefix = "eall_0size_meanOdist.Phi.site.plus.size.p.size.plus.dist")  # site- and size-dependent survival, size-and-distance-dependent recapture

# Compare model AICc
model_comp_0size_meanOdist = data.frame(model = c('eall_0size_meanOdist.Phi.dot.p.dot','eall_0size_meanOdist.Phi.dot.p.time','eall_0size_meanOdist.Phi.time.p.dot',
                                                  'eall_0size_meanOdist.Phi.dot.p.dist','eall_0size_meanOdist.Phi.size.p.size','eall_0size_meanOdist.Phi.size.p.dist',
                                                  'eall_0size_meanOdist.Phi.size.p.size.plus.dist','eall_0size_meanOdist.Phi.dot.p.size.plus.dist',
                                                  'eall_0size_meanOdist.Phi.site.p.dot','eall_0size_meanOdist.Phi.site.p.size.plus.dist',
                                                  'eall_0size_meanOdist.Phi.site.p.dist','eall_0size_meanOdist.Phi.site.plus.size.p.size.plus.dist'),
                                        AICc = c(eall_0size_meanOdist.Phi.dot.p.dot$results$AICc, eall_0size_meanOdist.Phi.dot.p.time$results$AICc, eall_0size_meanOdist.Phi.time.p.dot$results$AICc,
                                                 eall_0size_meanOdist.Phi.dot.p.dist$results$AICc, eall_0size_meanOdist.Phi.size.p.size$results$AICc, eall_0size_meanOdist.Phi.size.p.dist$results$AICc,
                                                 eall_0size_meanOdist.Phi.size.p.size.plus.dist$results$AICc, eall_0size_meanOdist.Phi.dot.p.size.plus.dist$results$AICc,
                                                 eall_0size_meanOdist.Phi.site.p.dot$results$AICc, eall_0size_meanOdist.Phi.site.p.size.plus.dist$results$AICc,
                                                 eall_0size_meanOdist.Phi.site.p.dist$results$AICc, eall_0size_meanOdist.Phi.site.plus.size.p.size.plus.dist$results$AICc),
                                        stringsAsFactors = FALSE)

# Arrange by increasing AICc
model_comp_0size_meanOdist <- model_comp_0size_meanOdist %>%
  mutate(dAICc = min(AICc) - AICc) %>%
  arrange(-dAICc)

#################### Old code: ####################
##### Script clean-up

#save(site_specific_surv, file=here::here("Data/Script_outputs", "site_specific_surv.RData"))

# # Phi (by site) - need to figure out what this represents, surv for mean size of site? surv for mean size? surv at size 0?
# pdf(file = here::here("Plots/PhiandpEstimates", "surv_by_site_Phisiteplussize_psizeplusdist.pdf"))
# ggplot(data = site_specific_surv_size, aes(x = site, y = estimate_prob, color = site, fill = site)) +
#   geom_point() + 
#   geom_linerange(aes(ymin=lcl_prob, ymax=ucl_prob)) +
#   ylab("survival probability") + ggtitle("Phi:site+size, p:size+dist, no size") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# dev.off()

# # List of sites to remove for Phi site+size model
# sites_no_recaps <- c("Caridad Cemetery", "Caridad Proper", "Sitio Lonas", "Sitio Tugas")

# # Mean survival by site (assume this is for mean size?)
# results_df_1 <- as.data.frame(eall_meanYsize_meanYdist.Phi.site.plus.size.p.size.plus.dist$results$beta)
# 
# site_specific_surv_size <- data.frame(site = site_vec, estimate = NA, se = NA, lcl = NA, ucl = NA)
# # Cabatoan is the intercept
# site_specific_surv_size$estimate[1] = results_df_1$estimate[1]
# site_specific_surv_size$se[1] = results_df_1$se[1]
# site_specific_surv_size$lcl[1] = results_df_1$lcl[1]
# site_specific_surv_size$ucl[1] = results_df_1$ucl[1]
# 
# # Go through rest of sites - need to be added to Cabatoan estimate
# for (i in 2:length(site_vec)) {
#   site_specific_surv_size$estimate[i] = results_df_1$estimate[1] + results_df_1$estimate[i]
#   site_specific_surv_size$se[i] = results_df_1$se[1] + results_df_1$se[i]
#   site_specific_surv_size$lcl[i] = results_df_1$lcl[1] + results_df_1$lcl[i]
#   site_specific_surv_size$ucl[i] = results_df_1$ucl[1] + results_df_1$ucl[i]
# }
# 
# # Convert back from logit scale
# site_specific_surv_size <- site_specific_surv_size %>%
#   mutate(estimate_prob = logit_recip(estimate),
#          lcl_prob = logit_recip(lcl),
#          ucl_prob = logit_recip(ucl))
# 

# ##### Do runs where taking the sites with no recaptures (Caridad Cemetery, Caridad Proper, Sitio Lonas, Sitio Tugas) out so can see if site+size works for Phi
# # Make data frame
# eall_meanYsize_meanYdist_sitesubset <- eall_meanYsize_meanYdist %>% filter(site %in% c("Cabatoan","Palanas","Sitio Baybayon"))
# 
# # ###### MAKE A FUNCTION TO DO THIS??!?!?! DIDN'T FINISH EDITING MODELS BELOW!
# # run_MARK_models <- function(encounters_df, prefix_vec) {
# #   
# #   models_list = list()
# #   
# #   # Process data for MARK runs
# #   df.processed = process.data(encounters_df, model="CJS", begin.time=2012)
# #   df.processed_site = process.data(encounters_df, model="CJS", begin.time=2012, groups="site")
# #   
# #   # Make ddl
# #   df.ddl = make.design.data(df.processed)
# #   df.ddl_site = make.design.data(df.processed_site)
# #   
# #   # Run models
# #   models_list[[1]] = mark(df.processed, df.ddl, model.parameters=list(Phi))
# # }
# # Prep for MARK runs
# # Process data (both with and without site as a group)
# eall_meanYsize_meanYdist_sitesubset.processed = process.data(eall_meanYsize_meanYdist_sitesubset, model='CJS', begin.time=2012)
# eall_meanYsize_meanYdist_sitesubset.processed_site = process.data(eall_meanYsize_meanYdist_sitesubset, model="CJS", begin.time=2012, groups="site")
# 
# # Make ddl
# eall_meanYsize_meanYdist_sitesubset.ddl = make.design.data(eall_meanYsize_meanYdist_sitesubset.processed)
# eall_meanYsize_meanYdist_sitesubset.ddl_site = make.design.data(eall_meanYsize_meanYdist_sitesubset.processed_site)
# 
# # Run models
# eall_meanYsize_meanYdist_sitesubset.Phi.dot.p.dot = mark(eall_meanYsize_meanYdist_sitesubset.processed, eall_meanYsize_meanYdist_sitesubset.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.dot.p.dot")  # constant survival and recap, no covariates
# eall_meanYsize_meanYdist_sitesubset.Phi.dot.p.time = mark(eall_meanYsize_meanYdist_sitesubset.processed, eall_meanYsize_meanYdist_sitesubset.ddl, model.parameters=list(Phi=Phi.dot, p=p.time), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.dot.p.time")  # constant survival, recapture time-dependent
# eall_meanYsize_meanYdist_sitesubset.Phi.time.p.dot = mark(eall_meanYsize_meanYdist_sitesubset.processed, eall_meanYsize_meanYdist_sitesubset.ddl, model.parameters=list(Phi=Phi.time, p=p.dot), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.time.p.dot")  # time-varying survival, recapture constant
# eall_meanYsize_meanYdist_sitesubset.Phi.dot.p.dist = mark(eall_meanYsize_meanYdist_sitesubset.processed, eall_meanYsize_meanYdist_sitesubset.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.dot.p.dist")  # constant survival, recapture distance-dependent
# eall_meanYsize_meanYdist_sitesubset.Phi.size.p.size = mark(eall_meanYsize_meanYdist_sitesubset.processed, eall_meanYsize_meanYdist_sitesubset.ddl, model.parameters=list(Phi=Phi.size, p=p.size), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.size.p.size")  # capture-size-dependent survival and recapture
# eall_meanYsize_meanYdist_sitesubset.Phi.size.p.dist = mark(eall_meanYsize_meanYdist_sitesubset.processed, eall_meanYsize_meanYdist_sitesubset.ddl, model.parameters=list(Phi=Phi.size, p=p.dist), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.size.p.dist")  # size-dependent survival, distance-dependent recapture
# eall_meanYsize_meanYdist_sitesubset.Phi.size.p.size.plus.dist = mark(eall_meanYsize_meanYdist_sitesubset.processed, eall_meanYsize_meanYdist_sitesubset.ddl, model.parameters=list(Phi=Phi.size, p=p.size.plus.dist), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.size.p.size.plus.dist")  # size-dependent survival, size-and-distance-dependent recapture
# eall_meanYsize_meanYdist_sitesubset.Phi.dot.p.size.plus.dist = mark(eall_meanYsize_meanYdist_sitesubset.processed, eall_meanYsize_meanYdist_sitesubset.ddl, model.parameters=list(Phi=Phi.dot, p=p.size.plus.dist), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.dot.p.size.plus.dist")  # constant survival, size-and-distance-dependent recapture
# eall_meanYsize_meanYdist_sitesubset.Phi.site.p.size.plus.dist = mark(eall_meanYsize_meanYdist_sitesubset.processed_site, eall_meanYsize_meanYdist_sitesubset.ddl_site, model.parameters=list(Phi=Phi.site, p=p.size.plus.dist), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.site.p.size.plus.dist")  # site-dependent survival, size-and-distance-dependent recapture
# eall_meanYsize_meanYdist_sitesubset.Phi.site.p.dot = mark(eall_meanYsize_meanYdist_sitesubset.processed_site, eall_meanYsize_meanYdist_sitesubset.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dot), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.site.p.dot")  # site-dependent survival, constant recapture
# eall_meanYsize_meanYdist_sitesubset.Phi.site.p.dist = mark(eall_meanYsize_meanYdist_sitesubset.processed_site, eall_meanYsize_meanYdist_sitesubset.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dist), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.site.p.dist")  # site-dependent survival, dist-related recapture
# 
# eall_meanYsize_meanYdist_sitesubset.Phi.site.plus.size.p.size.plus.dist = mark(eall_meanYsize_meanYdist_sitesubset.processed_site, eall_meanYsize_meanYdist_sitesubset.ddl_site, model.parameters=list(Phi=Phi.site.plus.size, p=p.size.plus.dist), prefix = "eall_meanYsize_meanYdist_sitesubset.Phi.site.plus.size.p.size.plus.dist")  # site- and size-dependent survival, size-and-distance-dependent recapture
# 

# ###### Survival by size, recapture both by size and distance (eall_mean.Phi.size.p.size.plus.dist)
# minsize = 1
# maxsize = 15
# size.values <- minsize+(0:30)*(maxsize-minsize)/30
# 
# mindist = 0
# # maxdist1 = max(eall$dist2016,eall$dist2017,eall$dist2018) #seems kind of high and probably due to an error
# maxdist2 = 200 #encompasses highest values in 2016, 2018 
# # dist.values1 = mindist+(0:30)*(maxdist1-mindist)/30
# dist.values2 = mindist+(0:30)*(maxdist2-mindist)/30
# 
# # Means
# eall_mean.Phi.size.p.size.plus.dist.results = as.data.frame(eall_mean.Phi.size.p.size.plus.dist$results$beta)
# 
# Phibysize_Phi.size.p.size.plus.dist_means <- data.frame(size = size.values) %>%
#   mutate(Phi_logit = eall_mean.Phi.size.p.size.plus.dist.results$estimate[1] + eall_mean.Phi.size.p.size.plus.dist.results$estimate[2]*size,
#          Phi_lcl_logit = eall_mean.Phi.size.p.size.plus.dist.results$lcl[1] + eall_mean.Phi.size.p.size.plus.dist.results$lcl[2]*size,
#          Phi_ucl_logit = eall_mean.Phi.size.p.size.plus.dist.results$ucl[1] + eall_mean.Phi.size.p.size.plus.dist.results$ucl[2]*size,
#          Phi = logit_recip(Phi_logit),
#          Phi_lcl = logit_recip(Phi_lcl_logit),
#          Phi_ucl = logit_recip(Phi_ucl_logit))
# 
# pbysize_Phi.size.p.size.plus.dist_means <- data.frame(size = size.values) %>%
#   mutate(p_logit = eall_mean.Phi.size.p.size.plus.dist.results$estimate[3] + eall_mean.Phi.size.p.size.plus.dist.results$estimate[4]*size,
#          p_lcl_logit = eall_mean.Phi.size.p.size.plus.dist.results$lcl[3] + eall_mean.Phi.size.p.size.plus.dist.results$lcl[4]*size,
#          p_ucl_logit = eall_mean.Phi.size.p.size.plus.dist.results$ucl[3] + eall_mean.Phi.size.p.size.plus.dist.results$ucl[4]*size,
#          p = logit_recip(p_logit),
#          p_lcl = logit_recip(p_lcl_logit),
#          p_ucl = logit_recip(p_ucl_logit))
# 
# pbydist_Phi.size.p.size.plus.dist_means <- data.frame(dist = dist.values2) %>%
#   mutate(p_logit = eall_mean.Phi.size.p.size.plus.dist.results$estimate[3] + eall_mean.Phi.size.p.size.plus.dist.results$estimate[5]*dist,
#          p_lcl_logit = eall_mean.Phi.size.p.size.plus.dist.results$lcl[3] + eall_mean.Phi.size.p.size.plus.dist.results$lcl[5]*dist,
#          p_ucl_logit = eall_mean.Phi.size.p.size.plus.dist.results$ucl[3] + eall_mean.Phi.size.p.size.plus.dist.results$ucl[5]*dist,
#          p = logit_recip(p_logit),
#          p_lcl = logit_recip(p_lcl_logit),
#          p_ucl = logit_recip(p_ucl_logit))
# 
# # Zeros
# eall_0.Phi.size.p.size.plus.dist.results = as.data.frame(eall_0.Phi.size.p.size.plus.dist$results$beta)
# 
# Phibysize_Phi.size.p.size.plus.dist_0s <- data.frame(size = size.values) %>%
#   mutate(Phi_logit = eall_0.Phi.size.p.size.plus.dist.results$estimate[1] + eall_0.Phi.size.p.size.plus.dist.results$estimate[2]*size,
#          Phi_lcl_logit = eall_0.Phi.size.p.size.plus.dist.results$lcl[1] + eall_0.Phi.size.p.size.plus.dist.results$lcl[2]*size,
#          Phi_ucl_logit = eall_0.Phi.size.p.size.plus.dist.results$ucl[1] + eall_0.Phi.size.p.size.plus.dist.results$ucl[2]*size,
#          Phi = logit_recip(Phi_logit),
#          Phi_lcl = logit_recip(Phi_lcl_logit),
#          Phi_ucl = logit_recip(Phi_ucl_logit))
# 
# pbysize_Phi.size.p.size.plus.dist_0s <- data.frame(size = size.values) %>%
#   mutate(p_logit = eall_0.Phi.size.p.size.plus.dist.results$estimate[3] + eall_0.Phi.size.p.size.plus.dist.results$estimate[4]*size,
#          p_lcl_logit = eall_0.Phi.size.p.size.plus.dist.results$lcl[3] + eall_0.Phi.size.p.size.plus.dist.results$lcl[4]*size,
#          p_ucl_logit = eall_0.Phi.size.p.size.plus.dist.results$ucl[3] + eall_0.Phi.size.p.size.plus.dist.results$ucl[4]*size,
#          p = logit_recip(p_logit),
#          p_lcl = logit_recip(p_lcl_logit),
#          p_ucl = logit_recip(p_ucl_logit))
# 
# pbydist_Phi.size.p.size.plus.dist_0s <- data.frame(dist = dist.values2) %>%
#   mutate(p_logit = eall_0.Phi.size.p.size.plus.dist.results$estimate[3] + eall_0.Phi.size.p.size.plus.dist.results$estimate[5]*dist,
#          p_lcl_logit = eall_0.Phi.size.p.size.plus.dist.results$lcl[3] + eall_0.Phi.size.p.size.plus.dist.results$lcl[5]*dist,
#          p_ucl_logit = eall_0.Phi.size.p.size.plus.dist.results$ucl[3] + eall_0.Phi.size.p.size.plus.dist.results$ucl[5]*dist,
#          p = logit_recip(p_logit),
#          p_lcl = logit_recip(p_lcl_logit),
#          p_ucl = logit_recip(p_ucl_logit))


######## From testing, previous script
#' 
#' #### Models with distance (with 0s put in for distance pre-fish first caught), capture size, capture tail color, capture stage
#' # Data
#' # Using 0 for distances for years before fish are caught
#' eall_0 <- encounters_all_0s %>%
#'   select(ch, site, cap_size, cap_color, cap_stage, dist2012, dist2013, dist2014, dist2015, dist2016, dist2017, dist2018,
#'          size2012, size2013, size2014, size2015, size2016, size2017, size2018) %>%
#'   mutate(site = as.factor(site), cap_color = as.factor(cap_color), cap_stage = as.factor(cap_stage))
#' 
#' eall_0 <- eall_0[complete.cases(eall_0),] #need complete cases, sometimes NA in size
#' 
#' 
#' # Using mean dist in that year for years before fish are caught
#' eall_mean <- encounters_all_means %>%
#'   select(ch, site, cap_size, cap_color, cap_stage, dist2012, dist2013, dist2014, dist2015, dist2016, dist2017, dist2018,
#'          size2012, size2013, size2014, size2015, size2016, size2017, size2018) %>%
#'   mutate(site = as.factor(site), cap_color = as.factor(cap_color), cap_stage = as.factor(cap_stage))
#' 
#' eall_mean <- eall_mean[complete.cases(eall_mean),]
#' 
#' # Process data and make ddl
#' # first for data with mean filled in for NAs in both fish size and distance
#' eall_mean.processed = process.data(eall_mean, model='CJS', begin.time=2012)
#' eall_mean.processed_site = process.data(eall_mean, model="CJS", begin.time=2012, groups="site")
#' #eall_mean.processed_stage = process.data(eall_mean, model='CJS', begin.time=2012, groups='cap_stage')
#' #eall_mean.processed_color = process.data(eall_mean, model="CJS", begin.time=2012, groups="cap_color")
#' 
#' eall_mean.ddl = make.design.data(eall_mean.processed)
#' eall_mean.ddl_site = make.design.data(eall_mean.processed_site)
#' #eall_mean.ddl_stage = make.design.data(eall_mean.processed_stage)
#' #eall_mean.ddl_color = make.design.data(eall_mean.processed_color)
#' 
#' # then for data with 0s filled in for distance and size NAs
#' eall_0.processed = process.data(eall_0, model='CJS', begin.time=2012)
#' eall_0.processed_site = process.data(eall_0, model="CJS", begin.time=2012, groups="site")
#' #eall_0.processed_stage = process.data(eall_0, model='CJS', begin.time=2012, groups='cap_stage')
#' #eall_0.processed_color = process.data(eall_0, model="CJS", begin.time=2012, groups="cap_color")
#' 
#' eall_0.ddl = make.design.data(eall_0.processed)
#' eall_0.ddl_site = make.design.data(eall_0.processed_site)
#' #eall_0.ddl_stage = make.design.data(eall_0.processed_stage)
#' #eall_0.ddl_color = make.design.data(eall_0.processed_color)
#' 
#' # set models for Phi
#' Phi.dot = list(formula=~1, link="logit")
#' Phi.time = list(formula=~time, link="logit")
#' Phi.site = list(formula=~site, link="logit")
#' Phi.size = list(formula=~size, link='logit')
#' #Phi.stage = list(formula=~cap_stage, link="logit")
#' 
#' #Phi.capsize = list(formula=~cap_size, link="logit")
#' #Phi.color = list(formula=~cap_color, link="logit")
#' #Phi.size.plus.stage = list(formula=~cap_size+cap_stage, link='logit')
#' #Phi.size.plus.color = list(formula=~cap_size+cap_color, link='logit')
#' 
#' # set models for p
#' p.dot = list(formula=~1, link="logit")
#' p.time = list(formula=~time, link="logit")
#' p.dist = list(formula=~dist, link="logit")
#' p.site = list(formula=~site, link="logit")
#' p.size = list(formula=~size, link='logit')
#' p.size.plus.dist = list(formula=~size+dist, link='logit')
#' 
#' #p.size = list(formula=~cap_size, link="logit")
#' #p.stage = list(formula=~cap_stage, link="logit")
#' #p.size.plus.stage = list(formula=~cap_size+cap_stage, link="logit")
#' #p.time.plus.dist = list(formula=~time+dist, link="logit")
#' 
#' # run some models
#' # using mean-dist-for-NA and mean-size-for-NA (mean within each year) dataset
#' eall_mean.Phi.dot.p.dot = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot), prefix = 'eall_mean.Phi.dot.p.dot') #constant survival and recap, no covariates
#' eall_mean.Phi.dot.p.time = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.dot, p=p.time), prefix = 'eall_mean.Phi.dot.p.time') #constant survival, recapture time-dependent
#' eall_mean.Phi.time.p.dot = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.time, p=p.dot), prefix = 'eall_mean.Phi.time.p.dot') #time-varying survival, recapture constant
#' eall_mean.Phi.dot.p.dist = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist), prefix = 'eall_mean.Phi.dot.p.dist') #constant survival, recapture distance-dependent
#' eall_mean.Phi.size.p.size = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.size, p=p.size), prefix = 'eall_mean.Phi.size.p.size') #capture-size-dependent survival and recapture
#' eall_mean.Phi.size.p.dist = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.size, p=p.dist), prefix = 'eall_mean.Phi.size.p.dist') #size-dependent survival, distance-dependent recapture
#' eall_mean.Phi.size.p.size.plus.dist = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.size, p=p.size.plus.dist), prefix = 'eall_mean.Phi.size.p.size.plus.dist') #size-dependent survival, size-and-distance-dependent recapture
#' eall_mean.Phi.dot.p.size.plus.dist = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.dot, p=p.size.plus.dist), prefix = 'eall_mean.Phi.dot.p.size.plus.dist') #constant survival, size-and-distance-dependent recapture
#' eall_mean.Phi.site.p.size.plus.dist = mark(eall_mean.processed_site, eall_mean.ddl_site, model.parameters=list(Phi=Phi.site, p=p.size.plus.dist), prefix = 'eall_mean.Phi.site.p.size.plus.dist') #site-dependent survival, size-and-distance-dependent recapture
#' eall_mean.Phi.site.p.dot = mark(eall_mean.processed_site, eall_mean.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dot), prefix = 'eall_mean.Phi.site.p.dot') #site-dependent survival, constant recapture
#' #eall_mean.Phi.stage.p.size.plus.dist = mark(eall_mean.processed_stage, eall_mean.ddl_stage, model.parameters=list(Phi=Phi.stage, p=p.size.plus.dist), prefix = 'eall_mean.Phi.stage.p.size.plus.dist') #capture stage-dependent survival, size-and-distance-dependent recapture
#' #eall_mean.Phi.stage.p.dot = mark(eall_mean.processed_stage, eall_mean.ddl_stage, model.parameters=list(Phi=Phi.stage, p=p.dot), prefix = 'eall_mean.Phi.stage.p.dot') #capture stage-dependent survival, constant recapture
#' 
#' # using 0-dist-for-NA dataset
#' eall_0.Phi.dot.p.dot = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot), prefix = 'eall_0.Phi.dot.p.dot') #constant survival and recap, no covariates
#' eall_0.Phi.dot.p.time = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.dot, p=p.time), prefix = 'eall_0.Phi.dot.p.time') #constant survival, recapture time-dependent
#' eall_0.Phi.time.p.dot = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.time, p=p.dot), prefix = 'eall_0.Phi.time.p.dot') #time-varying survival, recapture constant
#' eall_0.Phi.dot.p.dist = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist), prefix = 'eall_0.Phi.dot.p.dist') #constant survival, recapture distance-dependent
#' eall_0.Phi.size.p.size = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.size, p=p.size), prefix = 'eall_0.Phi.size.p.size') #capture-size-dependent survival and recapture
#' eall_0.Phi.size.p.dist = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.size, p=p.dist), prefix = 'eall_0.Phi.size.p.dist') #size-dependent survival, distance-dependent recapture
#' eall_0.Phi.size.p.size.plus.dist = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.size, p=p.size.plus.dist), prefix = 'eall_0.Phi.size.p.size.plus.dist') #size-dependent survival, size-and-distance-dependent recapture
#' eall_0.Phi.dot.p.size.plus.dist = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.dot, p=p.size.plus.dist), prefix = 'eall_0.Phi.dot.p.size.plus.dist') #constant survival, size-and-distance-dependent recapture
#' eall_0.Phi.site.p.size.plus.dist = mark(eall_0.processed_site, eall_0.ddl_site, model.parameters=list(Phi=Phi.site, p=p.size.plus.dist), prefix = 'eall_0.Phi.site.p.size.plus.dist') #site-dependent survival, size-and-distance-dependent recapture
#' eall_0.Phi.site.p.dot = mark(eall_0.processed_site, eall_0.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dot), prefix = 'eall_0.Phi.site.p.dot') #site-dependent survival, constant recapture
#' #eall_0.Phi.stage.p.size.plus.dist = mark(eall_0.processed_stage, eall_0.ddl_stage, model.parameters=list(Phi=Phi.stage, p=p.size.plus.dist), prefix = 'eall_0.Phi.stage.p.size.plus.dist') #capture stage-dependent survival, size-and-distance-dependent recapture
#' #eall_0.Phi.stage.p.dot = mark(eall_0.processed_stage, eall_0.ddl_stage, model.parameters=list(Phi=Phi.stage, p=p.dot), prefix = 'eall_0.Phi.stage.p.dot') #capture stage-dependent survival, constant recapture
#' 
#' 
#' # Compare model AICc
#' model_comp_mean = data.frame(model = c('eall_mean.Phi.dot.p.dot','eall_mean.Phi.dot.p.time','eall_mean.Phi.time.p.dot',
#'                                        'eall_mean.Phi.dot.p.dist','eall_mean.Phi.size.p.size','eall_mean.Phi.size.p.dist',
#'                                        'eall_mean.Phi.size.p.size.plus.dist','eall_mean.Phi.dot.p.size.plus.dist',
#'                                        'eall_mean.Phi.site.p.dot','eall_mean.Phi.site.p.size.plus.dist'),
#'                              #'eall_mean.Phi.stage.p.dot','eall_mean.Phi.stage.p.size.plus.dist'),
#'                              AICc = c(eall_mean.Phi.dot.p.dot$results$AICc, eall_mean.Phi.dot.p.time$results$AICc, eall_mean.Phi.time.p.dot$results$AICc,
#'                                       eall_mean.Phi.dot.p.dist$results$AICc, eall_mean.Phi.size.p.size$results$AICc, eall_mean.Phi.size.p.dist$results$AICc,
#'                                       eall_mean.Phi.size.p.size.plus.dist$results$AICc, eall_mean.Phi.dot.p.size.plus.dist$results$AICc,
#'                                       eall_mean.Phi.site.p.dot$results$AICc, eall_mean.Phi.site.p.dot$results$AICc),
#'                              stringsAsFactors = FALSE)
#' #eall_mean.Phi.stage.p.dot$results$AICc, eall_mean.Phi.stage.p.size.plus.dist$results$AICc))
#' # , 
#' #                         Phi_intercept = c(eall_mean.Phi.dot.p.dot$results$beta$estimate[1], eall_mean.Phi.dot.p.time$results$beta$estimate[1],
#' #                                           eall_mean.Phi.time.p.dot$results$beta$estimate[1], eall_mean.Phi.dot.p.dist$results$beta$estimate[1],
#' #                                           eall_mean.Phi.size.p.size$results$beta$estimate[1], eall_mean.Phi.size.p.dist$results$beta$estimate[1]))
#' # ,
#' #                         p_intercept =  c(eall_mean.Phi.dot.p.dot$results$beta$estimate[2], eall_mean.Phi.dot.p.time$results$beta$estimate[2],
#' #                                          eall_mean.Phi.time.p.dot$results$beta$estimate[1], eall_mean.Phi.dot.p.dist$results$beta$estimate[1],
#' #                                          eall_mean.Phi.size.p.size$results$beta$estimate[1], eall_mean.Phi.size.p.dist$results$beta$estimate[1]))
#' 
#' model_comp_0s = data.frame(model = c('eall_0.Phi.dot.p.dot','eall_0.Phi.dot.p.time','eall_0.Phi.time.p.dot',
#'                                      'eall_0.Phi.dot.p.dist','eall_0.Phi.size.p.size','eall_0.Phi.size.p.dist',
#'                                      'eall_0.Phi.size.p.size.plus.dist','eall_0.Phi.dot.p.size.plus.dist',
#'                                      'eall_0.Phi.site.p.dot','eall_0.Phi.site.p.size.plus.dist'),
#'                            #'eall_0.Phi.stage.p.dot','eall_0.Phi.stage.p.size.plus.dist'),
#'                            AICc = c(eall_0.Phi.dot.p.dot$results$AICc, eall_0.Phi.dot.p.time$results$AICc, eall_0.Phi.time.p.dot$results$AICc,
#'                                     eall_0.Phi.dot.p.dist$results$AICc, eall_0.Phi.size.p.size$results$AICc, eall_0.Phi.size.p.dist$results$AICc,
#'                                     eall_0.Phi.size.p.size.plus.dist$results$AICc, eall_0.Phi.dot.p.size.plus.dist$results$AICc,
#'                                     eall_0.Phi.site.p.dot$results$AICc, eall_0.Phi.site.p.dot$results$AICc))
#' #eall_0.Phi.stage.p.dot$results$AICc, eall_0.Phi.stage.p.size.plus.dist$results$AICc))
#' 
#' # Arrange by increasing AICc
#' model_comp_mean <- model_comp_mean %>%
#'   mutate(dAICc = min(AICc) - AICc) %>%
#'   arrange(-dAICc)
#' 
#' model_comp_0s <- model_comp_0s %>%
#'   mutate(dAIC = min(AIC) - AIC) %>%
#'   arrange(-dAIC)
#' 
#' # organize output to get ready to make plots
#' ###### Constant survival and recapture prob, no covariates (eall_mean.Phi.dot.p.dot)
#' eall_mean.constant = as.data.frame(eall_mean.Phi.dot.p.dot$results$beta) %>% 
#'   mutate(param = c("Phi","p")) %>% #add a parameters column to make it easier to plot in ggplot
#'   mutate(upper = logit_recip(ucl), lower = logit_recip(lcl), est = logit_recip(estimate)) #do the reciprocal transform on upper and lower confidence limits (need to check what those are - 95? SE? and that transforming them just straight up is the right way to go)
#' 
#' eall_0.constant = as.data.frame(eall_0.Phi.dot.p.dot$results$beta) %>% 
#'   mutate(param = c("Phi","p")) %>% #add a parameters column to make it easier to plot in ggplot
#'   mutate(upper = logit_recip(ucl), lower = logit_recip(lcl), est = logit_recip(estimate)) #do the reciprocal transform on upper and lower confidence limits (need to check what those are - 95? SE? and that transforming them just straight up is the right way to go)
#' 
#' 
#' 
#' 
#' 
#' 
#' #################### OLD VERSION OF SCRIPT (Pre-Oct 2019) ####################
#' ### In Oct, tidied to take the encounter history creation out into Clownfish_encounters and distance-to-anem calculations into AnemDistFromDiveTrack
#' 
#' # Mark-recap analysis to estimate survival, recapture probability, and population size based on tags and genetic IDs
#' # Somewhat cleaned-up version of code from ClownfishSurvivalEstimates.R
#' 
#' #################### Set-up: ####################
#' source(here::here('Code', 'Constants_database_common_functions.R'))
#' # #Load relevant libraries
#' # library(RCurl) #allows running R scripts from GitHub
#' # library(RMySQL) #might need to load this to connect to the database?
#' library(dplyr)
#' # library(tidyr)
#' library(RMark)
#' # library(lubridate)
#' library(geosphere)
#' # #library(dbplyr)
#' # library(ggplot2)
#' # library(here)
#' # #library(rethinking)
#' 
#' #Load data input 
#' #load(file=here::here('Data','anemProcessed_tableIDs.RData')) #file with anems by anem_table_id with lat/lon info appended (use with gen+tag runs), generated by AnemLocations.R
#' #load(file=here::here("Data",'AnemAllInfowLatLon2.RData')) #file with anems, after 2018 anems matched, for anems with anem_ids (so use with tag-only runs), generated by AnemLocations.R
#' 
#' load(file=here::here('Data','anem.Processed_full.RData'))  # don't need to load this any more, anems_Processed is in Constants_database_common_functions
#' 
#' #load(file=here("Data", "AnemLatLonObsbyAnem.RData")) #file with lat-lon info appended - this has mean lat/lon for each anem across obs - not sure I need this here? Maybe use the mean lat/lon later?
#' 
#' #load(file=here("Data", "fish_Tagged.RData")) #file with distances appended
#' 
#' # Load files from new, simple growth analysis (Growth_analysis.R)
#' load(file = here::here("Data/Script_outputs", "recap_pairs_year.RData"))
#' load(file = here::here("Data/Script_outputs", "growth_info_estimate.RData"))
#' 
#' # Find mean Linf and K from those runs
#' Linf_mean <- mean(growth_info_estimate$Linf_est)  # 10.71 (old: 10.58)
#' k_mean <- mean(growth_info_estimate$k_est)  # 0.864 (old: 0.928)
#' 
#' # # From grow_1pair1month$results in Methods_notes.pdf
#' # Linf_Faber <- 10.62
#' # K_Faber <- 0.906
#' # s2error_Faber <- 0.676
#' # t0 <- 0  #just making this up for now... (time or age at which size is 0)
#' 
#' # Calculates second length, t_i is in terms of years - from Fabers section of Hampton et al
#' Fabers_model <- function(Linf, L_release, K, t_i, e_i) {
#'   growth <- (Linf - L_release)*(1 - exp(-(K*t_i))) + e_i
#'   length_out <- L_release + growth
#'   return(length_out)
#' }
#' 
#' # Growth-increment VBL
#' growthIncrementVBL <- function(Linf, L_release, K) {
#'   length_out = L_release + (Linf - L_release)*(1 - exp(-K))
#'   return(length_out)
#' }
#' #################### Functions: ####################
#' # # Functions and constants from my GitHub function/constant collection
#' # script <- getURL("https://raw.githubusercontent.com/pinskylab/Clownfish_data_analysis/master/Code/Common_constants_and_functions.R?token=AH_ZQJT5uCEwjgDGYOneY0W6Zdjol5axks5alHmBwA%3D%3D", ssl.verifypeer = FALSE)
#' # eval(parse(text = script))
#' # 
#' # # Functions from Michelle's GitHub helpers script, updated version in field repository
#' # script <- getURL("https://raw.githubusercontent.com/pinskylab/field/master/scripts/field_helpers.R", ssl.verifypeer = FALSE)
#' # eval(parse(text = script))
#' # # 
#' # # Finds the real parameter estimate from the logit estimate
#' # logit_recip <- function(logitval) {
#' #   recip = (exp(logitval))/(1 + exp(logitval))
#' #   return(recip)
#' # }
#' 
#' # Creates the summarized encounter history by fish_id: output is a data frame with 2 columns - fish_id (either capXX or tagXXXXX or sampleXXX) and summarized encounter history (i.e. 0010)
#' CreateEncounterSummary <- function(start.year, end.year, tagged.fish) { #start.year is first year a fish could have been tagged (either PIT or genetically), end.year is last year of sampling, tagged.fish is dataframe with all fish with cap_ids or tag_ids
#'   
#'   sample.years <- seq(start.year, end.year, 1) #make a vector of years the fish could have been seen
#'   encounters <- list(); #initialize an empty list to store the various encounter data frames by year
#'  
#'   for (i in 1:length(sample.years)) { #pull out encounter vector by tag for each year, store each as a data frame of tag ids and binary encounters in the encounter list
#'     year.name <- sample.years[i] #get year 
#'     var.name <- paste("sighted", as.character(sample.years[i]), sep=".") #create dynamic column names for encounters - sighted.[samplingyear]
#'     
#'     encounters[[i]] <- tagged.fish %>% #create data frames for each year that have a vector of fish_ids and a vector of encountered (1) or didn't (0) in 
#'       group_by(fish_indiv) %>% 
#'       summarise(!!var.name := ifelse(sum(year == sample.years[i])>0, 1, 0)) #encounter history
#'   }
#'   
#'   encounters.out <- as.data.frame(encounters[[1]]) #seed summary data frame with list of tag ids and encounter in 1st year
#'   colnames(encounters.out)<- c("fish_indiv","encounter.hist") #rename the columns so the encounter.hist column can be pasted to iteratively in next step
#'   
#'   for (i in 2:length(sample.years)) {
#'     encounters.out$encounter.hist <- paste(encounters.out$encounter.hist, encounters[[i]][[2]], sep="") #paste on the other encounter 1/0s so get overall encounter histories as strings
#'   }
#'   
#'   return(encounters.out)
#' }
#' 
#' #Go through table, match anem_table_id to entry in anem.Processed2, take lat/lon from there (should do this better, by using average lat/lons - partial work on this below) - from AnemDistFromDiveTrack.R but with anem_table_id instead of anem_id
#' # addLatLons <- function(anem_vec, anemdf, latorlon) {
#' #   out = rep(NA, length(anem_vec))
#' #   
#' #   for(i in 1:length(anem_vec)) {
#' #     testid = anem_vec[i]
#' #     
#' #     if(!is.na(testid)) { #if the anem_table_id isn't NA
#' #       matches = anemdf %>% filter(anem_table_id == testid) #see if there are matches - so are there no anem_table_ids in 2012? or does anem.Processed2 just not pull them correctly...
#' #       if(length(matches$dive_table_id) != 0) {
#' #         if(latorlon == "lat") {
#' #           out[i] = matches$lat[1]
#' #         } else if (latorlon == "lon") {
#' #           out[i] = matches$lon[1]
#' #         }
#' #       } else if (length(matches$dive_table_id) == 0) {
#' #         print(paste("No id matches for anem_table_id",testid, sep=" "))
#' #       }
#' #     }
#' #   }
#' #   return(out)
#' # }
#' 
#' addLatLons <- function(encounters_df, anemdf, latorlon) {
#'   nentries = length(encounters_df$capture_anem_id) #number of anems to go through, variable to use as short-hand length
#'   out = rep(NA, nentries) #output vector
#' 
#'     for(i in 1:nentries) {
#'     test_table_id = encounters_df$capture_anem_table_id[i] #anem_table_id of the anem in question
#'     test_id = encounters_df$capture_anem_id[i] #anem_id of the anem in question (if it has one)
#'     # 
#'     # # see if can get a lat/lon from the anem_id
#'     # if(is.na(test_id)) { #if the anem_id is NA
#'     #   if(is.na(test_table_id)) { #and the anem_table_id is NA
#'     #     print('Both anem_id',)
#'     #   }
#'     # }
#'     # 
#' 
#'     if(!is.na(test_id)) { #if the anem_id isn't NA
#'       matches = anemdf %>% filter(anem_id == test_id) #see if there are any matches
#'       if(length(matches$dive_table_id != 0)) { #if there are matches
#'         if(latorlon == "lat") {
#'           out[i] = matches$lat[1] #put the lat coord in if lat
#'         } else if (latorlon == "lon") {
#'           out[i] = matches$lon[1] #put the lon coord in if lon
#'         }
#'       }
#'     }
#'     
#'    # check to see if got a coordinate from the anem_id, if not try the anem_table_id
#'     if(is.na(out[i])) { #if the coordinate is still NA
#'       if(!is.na(test_table_id)) { #and there is an anem_table_id
#'         matches = anemdf %>% filter(anem_table_id == test_table_id) #see if there are any matches
#'         if(length(matches$dive_table_id != 0)) { #if there are matches
#'           if(latorlon == "lat") {
#'             out[i] = matches$lat[1] #put the lat coord in if lat
#'           } else if (latorlon == "lon") {
#'             out[i] = matches$lon[1] #put the lon coord in if lon
#'           }
#'         }
#'       }
#'     }
#'     
#'     # check if either of those got a coordinate
#'     if(is.na(out[i])) {
#'       print(paste("No coordinates for anem_id",test_id,'or anem_table_id',test_table_id, sep=" "))
#'     }
#'   }
#'   return(out)
#' }
#'     
#' encounters_all$dist_2012 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2012, site_visits)
#' 
#' # Find the min distance from the anem where the fish was caught or last caught to a recorded GPS point in each year - from AnemDistFromDiveTrack.R
#' findDist <- function(sitevec, latvec, lonvec, capyearvec, gpxdf, divedf, sampleyear, site_visits){
#'   out <- rep(NA, length(sitevec))
#'   nsteps <- length(sitevec)
#'   
#'   for(i in 1:nsteps) {
#'     testsite = sitevec[i]
#'     testlat = as.numeric(latvec[i])
#'     testlon = as.numeric(lonvec[i])
#'     testcapyear = as.numeric(capyearvec[i])
#'     
#'     if(sampleyear >= testcapyear) {
#'       if(!is.na(testlat)) {
#'         visited <- site_visits$sampled[which(site_visits$site == testsite & site_visits$year == sampleyear)] # pull out the sampled value from this year and site
#'         if(visited == 1) { # if the sample site was visited that year, pull the dives from that site
#'           relevantdives <- divedf %>%
#'             filter(site == testsite, as.numeric(year) == sampleyear)
#'         } else { #if not, pull dives from sites 1 to the north and south
#'           site_N <- site_visits$site[which(site_visits$site == testsite & site_visits$year == sampleyear) - 1] #site to the north
#'           site_S <- site_visits$site[which(site_visits$site == testsite & site_visits$year == sampleyear) + 1] #site to the south
#'           
#'           relevantdives <- divedf %>%
#'             filter(site %in% c(site_N, site_S), as.numeric(year) == sampleyear)
#'         }
#'        
#'         #filter out gps points from those dives
#'         # relevantgps <- gpxdf %>%
#'         #   filter(gps_date %in% relevantdives$date & gps_year == sampleyear)
#'         relevantgps <- gpxdf %>%
#'           filter(gps_date %in% relevantdives$dive_date)
#'         
#'         distvec <- rep(NA, length(relevantgps$lat)) #set a place to store distances
#'         #distvec2 <- rep(NA, length(relevantgps$lat))
#'         
#'         for(j in 1:length(distvec)) { #go through the list of coordinates
#'           #distvec2[j] = distGeo(c(testlon, testlat), c(relevantgps$lon[j], relevantgps$lat[j]))
#'           distvec[j] = distHaversine(c(testlon,testlat), c(as.numeric(relevantgps$lon[j]), as.numeric(relevantgps$lat[j])))
#'           #print(paste(j, "of", length(distvec)))
#'         }
#'         out[i] = min(distvec)
#'         #print(min(distvec2))
#'       }
#'     } else {
#'       out[i] = NA
#'     }
#'     print(paste(i, "out of", length(sitevec), "and", out[i]))
#'   } 
#'   return(out)
#' }
#' 
#' # library(geosphere)
#' # mat <- distm(list1[,c('longitude','latitude')], list2[,c('longitude','latitude')], fun=distVincentyEllipsoid)
#' # 
#' # # Create a distance column in meters from a data.frame that has both points 
#' # loc.df$dist <- distGeo(loc.df[,c('lon1', 'lat1')], loc.df[,c('lon2', 'lat2')])
#' 
#' #################### Running things: ####################
#' 
#' # #### I think I can comment out the stuff between this and the next comment like this - check! (line 263)
#' # ##### Pulling and setting up data
#' # leyte <- read_db("Leyte")
#' # 
#' # # Pull all clownfish observations from the caught-clownfish table (don't need the ones only observed b/c won't be able to link them, right?)
#' # allfish_fish <- leyte %>%
#' #   tbl('clownfish') %>%
#' #   select(fish_table_id, anem_table_id, fish_spp, sample_id, gen_id, anem_table_id, recap, tag_id, color, size, notes) %>%
#' #   collect() %>%
#' #   filter(fish_spp == 'APCL')
#' # 
#' # # and their corresponding anemones
#' # allfish_anems <- leyte %>%
#' #   tbl('anemones') %>%
#' #   select(anem_table_id, dive_table_id, anem_obs, anem_id, old_anem_id) %>%
#' #   collect() %>%
#' #   filter(anem_table_id %in% allfish_fish$anem_table_id)
#' #   
#' # # and the corresponding dives
#' # allfish_dives <- leyte %>%  
#' #   tbl("diveinfo") %>%
#' #   select(dive_table_id, dive_type, date, site, gps) %>%
#' #   collect() %>%
#' #   filter(dive_table_id %in% allfish_anems$dive_table_id) %>%
#' #   mutate(year = as.integer(substring(date, 1, 4)))
#' #   
#' # # then join them together
#' # allfish <- left_join(allfish_fish, allfish_anems, by="anem_table_id")
#' # allfish <- left_join(allfish, allfish_dives, by="dive_table_id")
#' # 
#' # # Make size numeric (rather than a chr) so can do means and such
#' # allfish$size <- as.numeric(allfish$size) #make size numeric (rather than a chr) so can do means and such
#' # 
#' # # Pull out dive info (for assigning distances to anems each year)
#' # dive.Info <- leyte %>%
#' #   tbl("diveinfo") %>%
#' #   select(dive_table_id, dive_type, date, site, gps) %>%
#' #   collect() %>%
#' #   mutate(year = as.integer(substring(date,1,4))) 
#' # 
#' # # Pull out GPS info (for assigning distances to anems each year)
#' # gps.Info <- leyte %>%
#' #   tbl("GPX") %>%
#' #   select(lat, lon, time, unit) %>%
#' #   collect(n = Inf) %>%
#' #   mutate(obs_time = force_tz(ymd_hms(time), tzone = "UTC")) %>% #tell it that it is in UTC time zone
#' #   mutate(month = month(obs_time), #and separate out useful components of the time (this and line above largely from Michelle's assign_db_gpx function)
#' #          day = day(obs_time), 
#' #          hour = hour(obs_time), 
#' #          min = minute(obs_time), 
#' #          sec = second(obs_time), 
#' #          year = year(obs_time)) %>%
#' #   separate(time, into = c("date", "time"), sep = " ") #pull out date separately as a chr string too
#' # 
#' # ##### I think I can comment out the stuff above this line... check!
#' 
#' ##### Join together recaptures of the same fish, whether through genetic recapture or tag recapture
#' 
#' # Right now (1/7/18) - gen_id is the same for instances where it is really certain the fish is the same
#' # If everything checks out but the fish has moved sites, the gen_id is different but the words "genetic recapture" are in the notes
#' # Michelle thinks there is some tag loss
#' 
#' # Create several data sets:
#' # 1) all recaptures, genetic + tag, but only certain ones
#' # 2) all recaptures, genetic + tag, but including uncertain (site-switching) ones
#' # 3) just tag recaptures
#' # 4) all recaptures, genetic + tag, only certain ones, size estimated for unrecaught fish
#' # 5) all recaptures, genetic + tag, including uncertain (site-switching) ones, size estimated for unrecaught fish
#' # 6) just tag recaptures, size estimated for unrecaught fish
#' 
#' ##### Pull out fish that are marked in some way (either tag or genetic) - 3997 fish
#' # allfish_mark <- allfish %>%
#' #   filter(!is.na(gen_id) | (tag_id != 'NA' & !is.na(tag_id))) %>% #pull out fish "tagged" in any way, either PIT or via genetic sample
#' #   mutate(fish_id = case_when(tag_id != 'NA' ~ paste('tag', tag_id, sep=''), #if has a tag_id, use it in the fish_id
#' #                              tag_id == 'NA' & !is.na(gen_id) ~ paste('gen', gen_id, sep=''))) # %>% #if it doesn't have a tag_id, use the gen_id in the fish_id - but some of these fish might later have a tag_id
#' 
#' # allfish_mark <- allfish_caught %>%  # edited to use allfish_caught (which is generated by Constants_database_common_functions)
#' #   filter(!is.na(gen_id) | (tag_id != 'NA' & !is.na(tag_id))) %>%  # pull out fish "tagged" in any way, either PIT or via genetic sample
#' #   mutate(fish_id = case_when(gen_id != 'NA' ~ paste('gen', gen_id, sep=''),  # if fish has a gen_id, use that as the fish_id
#' #                              is.na(gen_id) ~ paste('tag', tag_id, sep='')))  # f it doesn't have a gen_id, use the tag_id
#' 
#' # Pull out all fish marked in some way
#' marked_fish <- allfish_caught %>%
#'   filter(!is.na(fish_indiv)) 
#' saveRDS(marked_fish, file = here::here("Data/Script_outputs", "marked_fish.RData"))
#' 
#' # Exploring recaptures - gen_id 1394 gets genetically id-ed 3 times (2015, May 2016, June 2016), getting a new PIT tag each time, then gets recaught with the third tag again in 2018... nuts!
#' 
#' # mark_set <- as.data.frame(table((allfish_mark %>% group_by(gen_id))$gen_id))
#' # 
#' # #### Set one: all recaptures, genetic + tag, only certain (non-site-switching) genetic recaps
#' # allfish_set1 <- allfish_mark
#' # 
#' # # Link up genetic ids with tagged ones, if posible - this time, going to switch to the gen_id rather than the tag, since seems like fish sometimes have multiple tags... should look into how common that is...
#' # for(i in 1:length(allfish_set1$fish_id)) {
#' #   if(substring(allfish_set1$fish_id[i],1,3) == 'tag') {  # if it has a fish_id based on tag_id rather than gen_id...
#' #     tag_id_val <- allfish_set1$tag_id[i]  # pull out the tag_id
#' #     
#' #     matches <- allfish_set1 %>%
#' #       filter(tag_id == tag_id_val)  # filter out any other cases that match that tag_id
#' #     
#' #     for(j in 1:length(matches$fish_id)) {
#' #      if(!is.na(matches$gen_id[j])) {  # if any of the captures with that tag_id have a gen_id
#' #        allfish_set1$fish_id[i] <- matches$fish_id[j]  # update the fish_id to the fish_id based on the match with a gen_id
#' #      }
#' #     }
#' #   }
#' # }
#' # 
#' # marked_fish <- allfish_set1
#' # saveRDS(marked_fish, file = here::here("Data/Script_outputs", "marked_fish.RData"))
#' # This seems to be working pretty well at catching all tags that have the same gen_id at some point, still 721 tag observations (617 individual tags) that are not associated with gen_ids (looks like sequencing failed on first capture, then didn't take samples after that b/c was tagged)
#' # For fish that still have a tag_id-based fish_id, see if other captures of that tag
#'   
#' # for each gen_id, figure out if there are any associated tag_ids
#' # for those cases, put the fish_id as the gen_id
#' # need to go through each of those tag_ids and find cases when fish got caught just based on the tag
#'   
#' # # Do some spot-checking
#' # allfish_set1 %>% filter(gen_id == 1394) %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # allfish_set1 %>% filter(fish_id == 'gen1394') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # 
#' # set1 <- as.data.frame(table(allfish_set1$fish_id))
#' # 
#' # allfish_set1 %>% filter(substring(fish_id,1,3) == 'tag') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)  # still about 720 with tags but no gen_id - looks like cases where samples were taken but genotyping failed?
#' # 
#' # # Spot-checking, looks like in these cases, the gen_id as fish_id percolates through multiple tags so long they are linked to the gen_id at some point
#' # allfish_set1 %>% filter(fish_id == 'gen1015') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # allfish_set1 %>% filter(tag_id == '986112100170625') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # allfish_set1 %>% filter(tag_id == '985153000406699') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # allfish_set1 %>% filter(fish_id == 'gen1019') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # allfish_set1 %>% filter(tag_id == '982000411818704') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # allfish_set1 %>% filter(tag_id == '985153000404653') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # allfish_set1 %>% filter(fish_id == 'gen1020') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # allfish_set1 %>% filter(fish_id == 'gen1065') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # allfish_set1 %>% filter(fish_id == 'gen1071') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # allfish_set1 %>% filter(fish_id == 'gen1499') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # allfish_set1 %>% filter(fish_id == 'gen1496') %>% select(-fish_spp, -dive_table_id, -notes, -gps, -fish_table_id, -anem_table_id)
#' # 
#' # # For cases that still have a tag_id-based fish_id, double-check that there aren't other captures of those fish with that tag that have a gen_id
#' # tag_set1 <- allfish_set1 %>% filter(substring(fish_id,1,3) == 'tag')  # still about 720 with tags but no gen_id - looks like cases where samples were taken but genotyping failed?
#' # tag_set1_table <- as.data.frame(table(tag_set1$fish_id))  # 617 individual fish (or tags...)
#' # tag_set1 %>% filter(fish_id == "tag982000411818572")
#' # tag_set1 %>% filter(!is.na(gen_id))  # none that have a gen_id
#' # 
#' # tag_ids_set1 <- allfish_set1 %>%
#' #   filter(tag_id %in% tag_set1$tag_id)  # pull out any fish with the tags that are in the tag_set1 tag_ids
#' # tag_ids_set1 %>% filter(!is.na(gen_id))  # do any of them have a gen_id at any point? - looks like no
#' #   
#'   
#' 
#' ##### Prep data for MARK - make encounter histories  - haven't done tag-only yet so commented out for now
#' # Pull out all the fish_ids and their encounter history from 2012-2018, including new genetic data, add column to see which fish are IDed by gen at any point (gen) or just tag (tag) - set 1 above
#' # encounters_all <- CreateEncounterSummary(2012, 2018, allfish_set1) %>%  # use function to get encounter history - gives 3222 fish encountered (was 3049 before 2016-2018 genetic data added in)
#' #   mutate(id_type = case_when(substring(fish_id,1,1) == "t" ~ "tag", #1122 (gen at some point)
#' #                              substring(fish_id,1,1) == "g" ~ "gen")) #1927 (tag only)
#' 
#' encounters_all <- CreateEncounterSummary(2012, 2018, marked_fish)
#' 
#' #encounters_tag <- CreateEncounterSummary(2015, 2018, allfish_tag) #use function to get encounter history - gives 1927 fish encountered
#' 
#' # Make tables of encounter histories into data frames for some plotting
#' #genetic+tag recaps, only certain genetic matches (no site-switching)
#' encounters_all_table <- as.data.frame(table(encounters_all$encounter.hist)) %>%
#'   dplyr::rename(encounter.hist = Var1) %>%
#'   mutate(times_caught = (as.numeric(substring(encounter.hist,1,1)) + 
#'             as.numeric(substring(encounter.hist,2,2)) + 
#'             as.numeric(substring(encounter.hist,3,3)) + 
#'             as.numeric(substring(encounter.hist,4,4)) +
#'             as.numeric(substring(encounter.hist,5,5)) +
#'             as.numeric(substring(encounter.hist,6,6)) +
#'             as.numeric(substring(encounter.hist,7,7))))
#' # 
#' # #just tag recaps
#' # encounters_tag_table <- as.data.frame(table(encounters_tag$encounter.hist)) %>%
#' #   dplyr::rename(encounter.hist = Var1) %>%
#' #   mutate(times_caught = (as.numeric(substring(encounter.hist,1,1)) + 
#' #                            as.numeric(substring(encounter.hist,2,2)) + 
#' #                            as.numeric(substring(encounter.hist,3,3)) + 
#' #                            as.numeric(substring(encounter.hist,4,4))))
#' 
#' # # Check if recaptured sites are different - they are in some cases but seems like fish actually moved... going to go with first site where they were caught for site
#' # # Find and add in site - a bit of a convoluted method right now but figured I should check for fish caught at different sites, while I was editing this to see if it is much of a problem...
#' # site_info <- allfish_mark %>%
#' #   group_by(fish_id) %>%
#' #   distinct(fish_id, year, .keep_all = TRUE) %>% #just pull out one observation per year
#' #   mutate(nobs = n()) %>%
#' #   mutate(site1 = site[1]) %>%
#' #   mutate(site2 = case_when(nobs == 1 ~ site[1], #or NA? can't decide which is better...
#' #                            nobs > 1 ~ site[2])) %>%
#' #   mutate(site3 = case_when(nobs <= 2 ~ site[1],
#' #                            nobs > 2 ~ site[3])) %>%
#' #   mutate(site4 = case_when(nobs <= 3 ~ site[1],
#' #                            nobs > 3 ~ site[4])) %>%
#' #   mutate(site5 = case_when(nobs <= 4 ~ site[1], #think 5 is the max number of obs at this point...
#' #                            nobs > 4 ~ site[5]))
#' # 
#' # # Check if any of the sites of recapture are different - went through these 10/4/18 and they seem like fish that actually moved (either the tagged fish we know about or MagN/MagS/Cab fish)
#' # for(i in 1:length(site_info$fish_id)) {
#' #   if(site_info$site1[i] != site_info$site2[i] |
#' #      site_info$site1[i] != site_info$site3[i] |
#' #      site_info$site1[i] != site_info$site4[i] |
#' #      site_info$site1[i] != site_info$site5[i]) {
#' #     print(site_info$fish_id[i])
#' #   } 
#' # }
#' 
#' # Find fish traits: site, tail color, size, life stage (all at first time captured), join with encounter histories 
#' #first, for genetic and tagged fish combined
#' trait_info <- marked_fish %>% 
#'   group_by(fish_indiv) %>%
#'   arrange(year) %>%
#'   summarize(site = site[1],
#'             first_capture_year = min(year), #year this fish was first captured
#'             capture_size = size[1], #earlier did min(size) and didn't arrange by year, could go either way
#'             capture_color = color[1],
#'             capture_stage = sex[1])
#'             # capture_stage = case_when(capture_color == "YP" ~ "female", #matches breeding F color + rough size - should re-evaluate these capture stage cutoffs at some point
#'             #                           capture_size >= min_breeding_F_size & capture_color == "Y" ~ "female", #matches breeding F color but undifferentiated, rough size
#'             #                           capture_color == "O" ~ "male", #matches breeding M color
#'             #                           capture_color != "YP" & capture_color != "Y" & capture_color != "YR" & capture_color != "WR" & capture_size >= min_breeding_M_size ~ "male", #this is the most fishy...
#'             #                           capture_size < min_breeding_M_size & capture_color != "YP" & capture_color != "O" ~ "juvenile", #small, not breeding colors
#'             #                           is.na(capture_color) & capture_size <= min_breeding_M_size ~ 'juvenile',
#'             #                           is.na(capture_color) & capture_size > min_breeding_M_size & capture_size <= female_size_cutoff ~ 'unknown',
#'             #                           is.na(capture_color) & capture_size > female_size_cutoff ~ 'female',
#'             #                           capture_color == "YR" & capture_size >= breeding_F_YR_cutoff ~ 'female',
#'             #                           capture_color == "YR" & capture_size < breeding_F_YR_cutoff ~ 'juvenile',
#'             #                           is.na(capture_color) & is.na(capture_size) ~ 'unknown',
#'             #                           is.na(capture_size) & capture_color == "Y" ~ 'unknown'
#'             #                           ))
#' # #and now just for tagged fish...
#' # trait_info_tag <- allfish_tag %>% 
#' #   group_by(fish_id) %>%
#' #   arrange(year) %>%
#' #   summarize(site = site[1],
#' #             first_capture_year = min(year),
#' #             capture_size = size[1],
#' #             capture_color = color[1],
#' #             capture_stage = case_when(capture_color == "YP" ~ "female", #matches breeding F color + rough size - should re-evaluate these capture stage cutoffs at some point
#' #                             capture_size >= min_breeding_F_size & capture_color == "Y" ~ "female", #matches breeding F color but undifferentiated, rough size
#' #                             capture_color == "O" ~ "male", #matches breeding M color
#' #                             capture_color != "YP" & capture_color != "Y" & capture_color != "YR" & capture_color != "WR" & capture_size >= min_breeding_M_size ~ "male", #this is the most fishy...
#' #                             capture_size < min_breeding_M_size & capture_color != "YP" & capture_color != "O" ~ "juvenile", #small, not breeding colors
#' #                             is.na(capture_color) & capture_size <= min_breeding_M_size ~ 'juvenile',
#' #                             is.na(capture_color) & capture_size > min_breeding_M_size & capture_size <= female_size_cutoff ~ 'unknown',
#' #                             is.na(capture_color) & capture_size > female_size_cutoff ~ 'female',
#' #                             capture_color == "YR" & capture_size >= breeding_F_YR_cutoff ~ 'female',
#' #                             capture_color == "YR" & capture_size < breeding_F_YR_cutoff ~ 'juvenile',
#' #                             is.na(capture_color) & is.na(capture_size) ~ 'unknown',
#' #                             is.na(capture_size) & capture_color == "Y" ~ 'unknown'))
#' 
#' # Join trait info with encounter histories
#' encounters_all <- left_join(encounters_all, trait_info, by="fish_indiv") #gen+tag recaps
#' # encounters_tag <- left_join(encounters_tag, trait_info_tag, by="fish_id") #tag-only recaps
#' 
#' ##### Find distance from anemone in previous years
#' encounters_all <- encounters_all %>%
#'   mutate(capture_anem_id = NA,
#'          capture_anem_table_id = NA,
#'          capture_anem_lat = NA,
#'          capture_anem_lon = NA,
#'          dist_2012 = NA,
#'          dist_2013 = NA,
#'          dist_2014 = NA,
#'          dist_2015 = NA,
#'          dist_2016 = NA,
#'          dist_2017 = NA,
#'          dist_2018 = NA)
#' 
#' # encounters_all <- encounters_all %>%
#' #   mutate(capture_anem_id = rep(NA, length(fish_indiv)), #anem_id of anem where fish was first captured
#' #          capture_anem_table_id = rep(NA, length(fish_indiv)), #anem_table_id of anem where fish was first captured
#' #          capture_anem_lat = rep(NA, length(fish_indiv)), #lat coordinates of first-capture anem
#' #          capture_anem_lon = rep(NA, length(fish_indiv)), #lon coordinates of first-capture anem
#' #          dist_2012 = rep(NA, length(fish_id)), 
#' #          dist_2013 = rep(NA, length(fish_id)),
#' #          dist_2014 = rep(NA, length(fish_id)),
#' #          dist_2015 = rep(NA, length(fish_id)),
#' #          dist_2016 = rep(NA, length(fish_id)),
#' #          dist_2017 = rep(NA, length(fish_id)),
#' #          dist_2018 = rep(NA, length(fish_id))) 
#'            
#' # Go through list of tagged fish, find distance from anem where first caught for all years after capture (originally, updated anem each year to where it had been caught (code in AmenDistFromDiveTrack) but doesn't that bias fish that have been recaught to have lower distances?)
#' for(i in 1:length(encounters_all$fish_indiv)) {
#'   
#'   # #Check for multiple anems recorded in capture year, print message if find some - this just records mulitple observations, not actually different anem_ids... amend to fix that...
#'   # if(length((allfish_mark %>% filter(fish_id == encounters_all$fish_id[i], year == encounters_all$first_capture_year[i]))$anem_id) > 1) { 
#'   #   mult_anems_test <- (allfish_mark %>% filter(fish_id == encounters_all$fish_id[i], year == encounters_all$first_capture_year[i]))$anem_id
#'   #   for(l in 1:length(mult_anems_test))
#'   # 
#'   #   print(paste("Multiple matching anems for fish",encounters_all$fish_id[i],"in year",encounters_all$first_capture_year[i]))
#'   # }
#'   
#'   # Record capture anem (id and anem_table_id) in first capture year and find lat/lon coordinates
#'   encounters_all$capture_anem_id[i] = (marked_fish %>% filter(fish_indiv == encounters_all$fish_indiv[i], year == encounters_all$first_capture_year[i]))$anem_id[1] #record the anem where fish was first captured
#'   encounters_all$capture_anem_table_id[i] = (marked_fish %>% filter(fish_indiv == encounters_all$fish_indiv[i], year == encounters_all$first_capture_year[i]))$anem_table_id[1]
#' }
#' 
#' # Now fill in the lat/lon coordinates for the capture anem
#' #encounters_all$capture_anem_lat <- addLatLons(encounters_all, anems_Processed, "lat")
#' encounters_all$capture_anem_lat <- addLatLons(encounters_all, anem.Processed_full, "lat") #now 5 anems without coordinates (which seems much more reasonable) and same ones for both lat and lon: anem_table_id 363, 364, 365, 366, 8094, none with anem_id
#' encounters_all$capture_anem_lon <- addLatLons(encounters_all, anem.Processed_full, "lon")
#' # talked to Michelle about these ones - probably from when the battery died at one point (2012) or the gps unit was submerged (2016)
#' 
#' # replace the NA in site_visits with 0s
#' site_visits$sampled[is.na(site_visits$sampled)] <- 0  
#' 
#' ## NEED TO RUN THESE OVERNIGHT - TOO SLOW TO DO DURING THE DAY!
#' # Calculate the distances from tracks in each year to the capture anem (NA for years prior to first capture year) - at some point, should check into why they're not all 0 or NA in 2012...
#' encounters_all$dist_2012 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2012, site_visits)
#' encounters_all$dist_2013 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2013, site_visits)
#' encounters_all$dist_2014 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2014, site_visits)
#' encounters_all$dist_2015 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2015, site_visits)
#' encounters_all$dist_2016 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2016, site_visits)
#' encounters_all$dist_2017 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2017, site_visits)
#' encounters_all$dist_2018 <- findDist(encounters_all$site, encounters_all$capture_anem_lat, encounters_all$capture_anem_lon, encounters_all$first_capture_year, gps_Info, dives_db, 2018, site_visits)
#' 
#' #test <- findDist(encounters_all$site[1:10], encounters_all$capture_anem_lat[1:10], encounters_all$capture_anem_lon[1:10], encounters_all$first_capture_year[1:10], gps_Info, dives_db, 2012, site_visits)
#' #save(encounters_all, file=here('Data','encounters_all.RData'))
#' saveRDS(encounters_all, file=here('Data','encounters_all_with_dist.RData'))
#' 
#' # save just distances
#' min_survey_dist_to_anems <- encounters_all %>% select(first_capture_year, capture_anem_id, capture_anem_table_id, capture_anem_lat, capture_anem_lon,
#'                                                      dist_2012, dist_2013, dist_2014, dist_2015, dist_2016, dist_2017, dist_2018)
#' 
#' ########## LOAD encounters_all AND START HERE!!!
#' load(file=here('Data', 'encounters_all.RData'))
#' ### NEED TO PUT IN SIZES AT EACH YEAR FOR EACH FISH (HOW TO DEAL WITH FISH CAUGHT TWICE IN ONE YEAR?)
#' ### THEN NEED TO FILL IN SIZES FOR FISH NOT CAUGHT USING VBL GROWTH CURVE
#' 
#' # Rename to fit mark input expectations
#' encounters_all <- encounters_all %>% dplyr::rename(ch = encounter.hist)
#' 
#' ### Add in sizes
#' 
#' # First, find sizes for captured fish
#' size_2012 <- marked_fish %>% filter(fish_indiv %in% encounters_all$fish_indiv) %>%
#'   filter(year == 2012) %>%
#'   group_by(fish_indiv) %>%
#'   summarize(size_2012 = mean(size, rm.na = TRUE))
#' 
#' size_2013 <- marked_fish %>% filter(fish_indiv %in% encounters_all$fish_indiv) %>%
#'   filter(year == 2013) %>%
#'   filter(size != 'NA') %>%
#'   group_by(fish_indiv) %>%
#'   summarize(size_2013 = mean(size, rm.na = TRUE))
#' 
#' size_2014 <- marked_fish %>% filter(fish_indiv %in% encounters_all$fish_indiv) %>%
#'   filter(year == 2014) %>%
#'   filter(size != 'NA') %>%
#'   group_by(fish_indiv) %>%
#'   summarize(size_2014 = mean(size, rm.na = TRUE))
#' 
#' size_2015 <- marked_fish %>% filter(fish_indiv %in% encounters_all$fish_indiv) %>%
#'   filter(year == 2015) %>%
#'   filter(size != 'NA') %>%
#'   group_by(fish_indiv) %>%
#'   summarize(size_2015 = mean(size, rm.na = TRUE))
#' 
#' size_2016 <- marked_fish %>% filter(fish_indiv %in% encounters_all$fish_indiv) %>%
#'   filter(year == 2016) %>%
#'   filter(size != 'NA') %>%
#'   group_by(fish_indiv) %>%
#'   summarize(size_2016 = mean(size, rm.na = TRUE))
#' 
#' size_2017 <- marked_fish %>% filter(fish_indiv %in% encounters_all$fish_indiv) %>%
#'   filter(year == 2017) %>%
#'   filter(size != 'NA') %>%
#'   group_by(fish_indiv) %>%
#'   summarize(size_2017 = mean(size, rm.na = TRUE))
#' 
#' size_2018 <- marked_fish %>% filter(fish_indiv %in% encounters_all$fish_indiv) %>%
#'   filter(year == 2018) %>%
#'   filter(size != 'NA') %>%
#'   group_by(fish_indiv) %>%
#'   summarize(size_2018 = mean(size, rm.na = TRUE))
#' 
#' encounters_all <- left_join(encounters_all, size_2012, by='fish_indiv')
#' encounters_all <- left_join(encounters_all, size_2013, by='fish_indiv')
#' encounters_all <- left_join(encounters_all, size_2014, by='fish_indiv')
#' encounters_all <- left_join(encounters_all, size_2015, by='fish_indiv')
#' encounters_all <- left_join(encounters_all, size_2016, by='fish_indiv')
#' encounters_all <- left_join(encounters_all, size_2017, by='fish_indiv')
#' encounters_all <- left_join(encounters_all, size_2018, by='fish_indiv')
#' 
#' # Fill in estimated sizes for unrecaptured (or unmeasured?) fish
#' mean_size <- mean(c(size_2012$size_2012, size_2013$size_2013, size_2014$size_2014, size_2015$size_2015, size_2016$size_2016, size_2017$size_2017, size_2018$size_2018))
#' 
#' encounters_all_sizes <- encounters_all  # just so don't mess things up and have to re-load encounters_all...
#' encounters_all_sizes_0 <- encounters_all 
#' 
#' 
#' # projectSize <- function(df, Linf_Faber, K_Faber, t_i, s2error_Faber) {
#' #   
#' # }
#' 
#' # Fill in mean size or 0 for all fish not caught yet (which MARK should ignore, right?), otherwise project based on growth curve
#' # 2012
#' for(i in 1:length(encounters_all_sizes$size_2012)) {
#'   size_test <- encounters_all_sizes$size_2012[i]
#'   
#'   if(is.na(size_test)) {
#'     encounters_all_sizes$size_2012[i] <- mean_size
#'     encounters_all_sizes_0$size_2012[i] <- 0
#'   }
#' }
#' 
#' # 2013
#' for(i in 1:length(encounters_all_sizes$size_2013)) {
#'   size_test <- encounters_all_sizes$size_2013[i]
#'   first_cap <- encounters_all_sizes$first_capture_year[i]
#'   prev_size <- encounters_all_sizes$size_2012[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2013) {
#'       encounters_all_sizes$size_2013[i] <- mean_size
#'       encounters_all_sizes_0$size_2013[i] <- 0
#'     } else {
#'       #encounters_all_sizes$size_2013[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'       #encounters_all_sizes_0$size_2013[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'       encounters_all_sizes$size_2013[i] <- growthIncrementVBL(Linf_mean, prev_size, k_mean)
#'       encounters_all_sizes_0$size_2013[i] <- growthIncrementVBL(Linf_mean, prev_size, k_mean)
#'       }
#'   }
#' }
#' 
#' # 2014
#' for(i in 1:length(encounters_all_sizes$size_2014)) {
#'   size_test <- encounters_all_sizes$size_2014[i]
#'   first_cap <- encounters_all_sizes$first_capture_year[i]
#'   prev_size <- encounters_all_sizes$size_2013[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2014) {
#'       encounters_all_sizes$size_2014[i] <- mean_size
#'       encounters_all_sizes_0$size_2014[i] <- 0
#'     } else {
#'       #encounters_all_sizes$size_2014[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'       encounters_all_sizes$size_2014[i] <- growthIncrementVBL(Linf_mean, prev_size, k_mean)
#'     }
#'   }
#' }
#' 
#' # 2015
#' for(i in 1:length(encounters_all_sizes$size_2015)) {
#'   size_test <- encounters_all_sizes$size_2015[i]
#'   first_cap <- encounters_all_sizes$first_capture_year[i]
#'   prev_size <- encounters_all_sizes$size_2014[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2015) {
#'       encounters_all_sizes$size_2015[i] <- mean_size
#'       encounters_all_sizes_0$size_2015[i] <- 0
#'     } else {
#'       #encounters_all_sizes$size_2015[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'       encounters_all_sizes$size_2015[i] <- growthIncrementVBL(Linf_mean, prev_size, k_mean)
#'     }
#'   }
#' }
#' 
#' # 2016
#' for(i in 1:length(encounters_all_sizes$size_2016)) {
#'   size_test <- encounters_all_sizes$size_2016[i]
#'   first_cap <- encounters_all_sizes$first_capture_year[i]
#'   prev_size <- encounters_all_sizes$size_2015[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2016) {
#'       encounters_all_sizes$size_2016[i] <- mean_size
#'       encounters_all_sizes_0$size_2016[i] <- 0
#'     } else {
#'       #encounters_all_sizes$size_2016[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'       encounters_all_sizes$size_2016[i] <- growthIncrementVBL(Linf_mean, prev_size, k_mean)
#'     }
#'   }
#' }
#' 
#' # 2017
#' for(i in 1:length(encounters_all_sizes$size_2017)) {
#'   size_test <- encounters_all_sizes$size_2017[i]
#'   first_cap <- encounters_all_sizes$first_capture_year[i]
#'   prev_size <- encounters_all_sizes$size_2016[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2017) {
#'       encounters_all_sizes$size_2017[i] <- mean_size
#'       encounters_all_sizes_0$size_2017[i] <- 0
#'     } else {
#'       #encounters_all_sizes$size_2017[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'       encounters_all_sizes$size_2017[i] <- growthIncrementVBL(Linf_mean, prev_size, k_mean)
#'     }
#'   }
#' }
#' 
#' # 2018
#' for(i in 1:length(encounters_all_sizes$size_2018)) {
#'   size_test <- encounters_all_sizes$size_2018[i]
#'   first_cap <- encounters_all_sizes$first_capture_year[i]
#'   prev_size <- encounters_all_sizes$size_2017[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2018) {
#'       encounters_all_sizes$size_2018[i] <- mean_size
#'       encounters_all_sizes_0$size_2018[i] <- 0
#'     } else {
#'       encounters_all_sizes$size_2018[i] <- growthIncrementVBL(Linf_mean, prev_size, k_mean)
#'     }
#'   }
#' }
#' ###### should really figure out how to project size and fill in missing pre-first capture sizes in a function...
#' # Now do it again, filling in 0s instead of mean dist for fish not caught yet (which MARK should ignore, right?), otherwise project based on growth curve
#' # 2012
#' for(i in 1:length(encounters_all_sizes_0$size_2012)) {
#'   size_test <- encounters_all_sizes_0$size_2012[i]
#'   
#'   if(is.na(size_test)) {
#'     encounters_all_sizes_0$size_2012[i] <- 0
#'   }
#' }
#' 
#' # 2013
#' for(i in 1:length(encounters_all_sizes_0$size_2013)) {
#'   size_test <- encounters_all_sizes_0$size_2013[i]
#'   first_cap <- encounters_all_sizes_0$first_capture_year[i]
#'   prev_size <- encounters_all_sizes_0$size_2012[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2013) {
#'       encounters_all_sizes_0$size_2013[i] <- 0
#'     } else {
#'       encounters_all_sizes_0$size_2013[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'     }
#'   }
#' }
#' 
#' # 2014
#' for(i in 1:length(encounters_all_sizes_0$size_2014)) {
#'   size_test <- encounters_all_sizes_0$size_2014[i]
#'   first_cap <- encounters_all_sizes_0$first_capture_year[i]
#'   prev_size <- encounters_all_sizes_0$size_2013[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2014) {
#'       encounters_all_sizes_0$size_2014[i] <- 0
#'     } else {
#'       encounters_all_sizes_0$size_2014[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'     }
#'   }
#' }
#' 
#' # 2015
#' for(i in 1:length(encounters_all_sizes_0$size_2015)) {
#'   size_test <- encounters_all_sizes_0$size_2015[i]
#'   first_cap <- encounters_all_sizes_0$first_capture_year[i]
#'   prev_size <- encounters_all_sizes_0$size_2014[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2015) {
#'       encounters_all_sizes_0$size_2015[i] <- 0
#'     } else {
#'       encounters_all_sizes_0$size_2015[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'     }
#'   }
#' }
#' 
#' # 2016
#' for(i in 1:length(encounters_all_sizes_0$size_2016)) {
#'   size_test <- encounters_all_sizes_0$size_2016[i]
#'   first_cap <- encounters_all_sizes_0$first_capture_year[i]
#'   prev_size <- encounters_all_sizes_0$size_2015[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2016) {
#'       encounters_all_sizes_0$size_2016[i] <- 0
#'     } else {
#'       encounters_all_sizes_0$size_2016[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'     }
#'   }
#' }
#' 
#' # 2017
#' for(i in 1:length(encounters_all_sizes_0$size_2017)) {
#'   size_test <- encounters_all_sizes_0$size_2017[i]
#'   first_cap <- encounters_all_sizes_0$first_capture_year[i]
#'   prev_size <- encounters_all_sizes_0$size_2016[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2017) {
#'       encounters_all_sizes_0$size_2017[i] <- 0
#'     } else {
#'       encounters_all_sizes_0$size_2017[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'     }
#'   }
#' }
#' 
#' # 2018
#' for(i in 1:length(encounters_all_sizes_0$size_2018)) {
#'   size_test <- encounters_all_sizes_0$size_2018[i]
#'   first_cap <- encounters_all_sizes_0$first_capture_year[i]
#'   prev_size <- encounters_all_sizes_0$size_2017[i]
#'   
#'   if(is.na(size_test)) {
#'     if(first_cap > 2018) {
#'       encounters_all_sizes_0$size_2018[i] <- 0
#'     } else {
#'       encounters_all_sizes_0$size_2018[i] <- Fabers_model(Linf_Faber, prev_size, K_Faber, 1, s2error_Faber)
#'     }
#'   }
#' }
#' 
#' # save(encounters_all_sizes, file=here::here('Data','encounters_all_sizes.RData'))
#' 
#' # # Add in distances to anem (from fish.Tagged, loaded above)
#' # encounters_all <- left_join(encounters_all, (fish.Tagged %>% select(tag_id, year_tagged, dist_2016, dist_2017, dist_2018)), by="tag_id") 
#' 
#' # can't have NAs in the covariate columns so changing to mean of distances in each year (should also try replacing with 0, since it shouldn't matter)
#' mean2012 <- mean(encounters_all$dist_2012, na.rm=TRUE)  # do I need this?
#' mean2013 <- mean(encounters_all$dist_2013, na.rm=TRUE)
#' mean2014 <- mean(encounters_all$dist_2014, na.rm=TRUE)
#' mean2015 <- mean(encounters_all$dist_2015, na.rm=TRUE)
#' mean2016 <- mean(encounters_all$dist_2016, na.rm=TRUE)
#' mean2017 <- mean(encounters_all$dist_2017, na.rm=TRUE)
#' mean2018 <- mean(encounters_all$dist_2018, na.rm=TRUE)
#' mean_dist_allyears <- mean(c(encounters_all$dist_2012, encounters_all$dist_2013, encounters_all$dist_2014,
#'                              encounters_all$dist_2015, encounters_all$dist_2016, encounters_all$dist_2017,
#'                              encounters_all$dist_2018), na.rm=TRUE)
#' mean_dist_no2012 <- mean(c(encounters_all$dist_2013, encounters_all$dist_2014,
#'                            encounters_all$dist_2015, encounters_all$dist_2016, encounters_all$dist_2017,
#'                            encounters_all$dist_2018), na.rm=TRUE)
#' 
#' # replace the NAs (for years before the fish was tagged) with the mean distance overall (previously was doing mean for that year)
#' encounters_all_means <- encounters_all_sizes %>%
#'   dplyr::rename(dist2012 = dist_2012, dist2013 = dist_2013, dist2014 = dist_2014, dist2015 = dist_2015,
#'          dist2016 = dist_2016, dist2017 = dist_2017, dist2018 = dist_2018) %>% #first, rename distance columns so easier for MARK to find
#'   dplyr::rename(cap_size = capture_size, cap_stage = capture_stage, cap_color = capture_color) %>% #turns out MARK has a character cap of 10 for covariates?
#'   dplyr::rename(size2012 = size_2012, size2013 = size_2013, size2014 = size_2014, size2015 = size_2015,
#'                 size2016 = size_2016, size2017 = size_2017, size2018 = size_2018) %>%
#'   mutate(dist2013 = replace(dist2013, is.na(dist2013), mean_dist_allyears),
#'          dist2014 = replace(dist2014, is.na(dist2014), mean_dist_allyears),
#'          dist2015 = replace(dist2015, is.na(dist2015), mean_dist_allyears),
#'          dist2016 = replace(dist2016, is.na(dist2016), mean_dist_allyears),
#'          dist2017 = replace(dist2017, is.na(dist2017), mean_dist_allyears),
#'          dist2018 = replace(dist2018, is.na(dist2018), mean_dist_allyears),
#'          dist2012 = replace(dist2012, is.na(dist2012), mean_dist_allyears)) 
#' 
#' # replace the NAs (for years before the fish was caught) with 0 - could also just calculate distance to that anem, based on the first tag number in the year the fish was tagged in the original script that finds the distance
#' encounters_all_0s <- encounters_all_sizes_0 %>%
#'   dplyr::rename(dist2012 = dist_2012, dist2013 = dist_2013, dist2014 = dist_2014, dist2015 = dist_2015,
#'          dist2016 = dist_2016, dist2017 = dist_2017, dist2018 = dist_2018) %>% #first, rename distance columns so easier for MARK to find
#'   dplyr::rename(cap_size = capture_size, cap_stage = capture_stage, cap_color = capture_color) %>% #turns out MARK has a character cap of 10 for covariates?
#'   dplyr::rename(size2012 = size_2012, size2013 = size_2013, size2014 = size_2014, size2015 = size_2015,
#'                 size2016 = size_2016, size2017 = size_2017, size2018 = size_2018) %>%
#'   mutate(dist2013 = replace(dist2013, is.na(dist2013), 0),
#'          dist2014 = replace(dist2014, is.na(dist2014), 0),
#'          dist2015 = replace(dist2015, is.na(dist2015), 0),
#'          dist2016 = replace(dist2016, is.na(dist2016), 0),
#'          dist2017 = replace(dist2017, is.na(dist2017), 0),
#'          dist2018 = replace(dist2018, is.na(dist2018), 0),
#'          dist2012 = replace(dist2012, is.na(dist2012), 0)) 
#' 
#' #encounters_all_means <- encounters_all_means %>% dplyr::rename(ch = encounter.hist)
#' 
#' #################### Running things: MARK models ####################
#' #################### ALL ENCOUNTERS (GENETIC + TAG) ####################
#' #### Models with distance (with 0s put in for distance pre-fish first caught), capture size, capture tail color, capture stage
#' # Data
#' # Using 0 for distances for years before fish are caught
#' eall_0 <- encounters_all_0s %>%
#'   select(ch, site, cap_size, cap_color, cap_stage, dist2012, dist2013, dist2014, dist2015, dist2016, dist2017, dist2018,
#'          size2012, size2013, size2014, size2015, size2016, size2017, size2018) %>%
#'   mutate(site = as.factor(site), cap_color = as.factor(cap_color), cap_stage = as.factor(cap_stage))
#' 
#' eall_0 <- eall_0[complete.cases(eall_0),] #need complete cases, sometimes NA in size
#' 
#' 
#' # Using mean dist in that year for years before fish are caught
#' eall_mean <- encounters_all_means %>%
#'   select(ch, site, cap_size, cap_color, cap_stage, dist2012, dist2013, dist2014, dist2015, dist2016, dist2017, dist2018,
#'          size2012, size2013, size2014, size2015, size2016, size2017, size2018) %>%
#'   mutate(site = as.factor(site), cap_color = as.factor(cap_color), cap_stage = as.factor(cap_stage))
#' 
#' eall_mean <- eall_mean[complete.cases(eall_mean),]
#' 
#' # Process data and make ddl
#' # first for data with mean filled in for NAs in both fish size and distance
#' eall_mean.processed = process.data(eall_mean, model='CJS', begin.time=2012)
#' eall_mean.processed_site = process.data(eall_mean, model="CJS", begin.time=2012, groups="site")
#' #eall_mean.processed_stage = process.data(eall_mean, model='CJS', begin.time=2012, groups='cap_stage')
#' #eall_mean.processed_color = process.data(eall_mean, model="CJS", begin.time=2012, groups="cap_color")
#' 
#' eall_mean.ddl = make.design.data(eall_mean.processed)
#' eall_mean.ddl_site = make.design.data(eall_mean.processed_site)
#' #eall_mean.ddl_stage = make.design.data(eall_mean.processed_stage)
#' #eall_mean.ddl_color = make.design.data(eall_mean.processed_color)
#' 
#' # then for data with 0s filled in for distance and size NAs
#' eall_0.processed = process.data(eall_0, model='CJS', begin.time=2012)
#' eall_0.processed_site = process.data(eall_0, model="CJS", begin.time=2012, groups="site")
#' #eall_0.processed_stage = process.data(eall_0, model='CJS', begin.time=2012, groups='cap_stage')
#' #eall_0.processed_color = process.data(eall_0, model="CJS", begin.time=2012, groups="cap_color")
#' 
#' eall_0.ddl = make.design.data(eall_0.processed)
#' eall_0.ddl_site = make.design.data(eall_0.processed_site)
#' #eall_0.ddl_stage = make.design.data(eall_0.processed_stage)
#' #eall_0.ddl_color = make.design.data(eall_0.processed_color)
#' 
#' # set models for Phi
#' Phi.dot = list(formula=~1, link="logit")
#' Phi.time = list(formula=~time, link="logit")
#' Phi.site = list(formula=~site, link="logit")
#' Phi.size = list(formula=~size, link='logit')
#' #Phi.stage = list(formula=~cap_stage, link="logit")
#' 
#' #Phi.capsize = list(formula=~cap_size, link="logit")
#' #Phi.color = list(formula=~cap_color, link="logit")
#' #Phi.size.plus.stage = list(formula=~cap_size+cap_stage, link='logit')
#' #Phi.size.plus.color = list(formula=~cap_size+cap_color, link='logit')
#' 
#' # set models for p
#' p.dot = list(formula=~1, link="logit")
#' p.time = list(formula=~time, link="logit")
#' p.dist = list(formula=~dist, link="logit")
#' p.site = list(formula=~site, link="logit")
#' p.size = list(formula=~size, link='logit')
#' p.size.plus.dist = list(formula=~size+dist, link='logit')
#' 
#' #p.size = list(formula=~cap_size, link="logit")
#' #p.stage = list(formula=~cap_stage, link="logit")
#' #p.size.plus.stage = list(formula=~cap_size+cap_stage, link="logit")
#' #p.time.plus.dist = list(formula=~time+dist, link="logit")
#' 
#' # run some models
#' # using mean-dist-for-NA and mean-size-for-NA (mean within each year) dataset
#' eall_mean.Phi.dot.p.dot = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot), prefix = 'eall_mean.Phi.dot.p.dot') #constant survival and recap, no covariates
#' eall_mean.Phi.dot.p.time = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.dot, p=p.time), prefix = 'eall_mean.Phi.dot.p.time') #constant survival, recapture time-dependent
#' eall_mean.Phi.time.p.dot = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.time, p=p.dot), prefix = 'eall_mean.Phi.time.p.dot') #time-varying survival, recapture constant
#' eall_mean.Phi.dot.p.dist = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist), prefix = 'eall_mean.Phi.dot.p.dist') #constant survival, recapture distance-dependent
#' eall_mean.Phi.size.p.size = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.size, p=p.size), prefix = 'eall_mean.Phi.size.p.size') #capture-size-dependent survival and recapture
#' eall_mean.Phi.size.p.dist = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.size, p=p.dist), prefix = 'eall_mean.Phi.size.p.dist') #size-dependent survival, distance-dependent recapture
#' eall_mean.Phi.size.p.size.plus.dist = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.size, p=p.size.plus.dist), prefix = 'eall_mean.Phi.size.p.size.plus.dist') #size-dependent survival, size-and-distance-dependent recapture
#' eall_mean.Phi.dot.p.size.plus.dist = mark(eall_mean.processed, eall_mean.ddl, model.parameters=list(Phi=Phi.dot, p=p.size.plus.dist), prefix = 'eall_mean.Phi.dot.p.size.plus.dist') #constant survival, size-and-distance-dependent recapture
#' eall_mean.Phi.site.p.size.plus.dist = mark(eall_mean.processed_site, eall_mean.ddl_site, model.parameters=list(Phi=Phi.site, p=p.size.plus.dist), prefix = 'eall_mean.Phi.site.p.size.plus.dist') #site-dependent survival, size-and-distance-dependent recapture
#' eall_mean.Phi.site.p.dot = mark(eall_mean.processed_site, eall_mean.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dot), prefix = 'eall_mean.Phi.site.p.dot') #site-dependent survival, constant recapture
#' #eall_mean.Phi.stage.p.size.plus.dist = mark(eall_mean.processed_stage, eall_mean.ddl_stage, model.parameters=list(Phi=Phi.stage, p=p.size.plus.dist), prefix = 'eall_mean.Phi.stage.p.size.plus.dist') #capture stage-dependent survival, size-and-distance-dependent recapture
#' #eall_mean.Phi.stage.p.dot = mark(eall_mean.processed_stage, eall_mean.ddl_stage, model.parameters=list(Phi=Phi.stage, p=p.dot), prefix = 'eall_mean.Phi.stage.p.dot') #capture stage-dependent survival, constant recapture
#' 
#' # using 0-dist-for-NA dataset
#' eall_0.Phi.dot.p.dot = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot), prefix = 'eall_0.Phi.dot.p.dot') #constant survival and recap, no covariates
#' eall_0.Phi.dot.p.time = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.dot, p=p.time), prefix = 'eall_0.Phi.dot.p.time') #constant survival, recapture time-dependent
#' eall_0.Phi.time.p.dot = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.time, p=p.dot), prefix = 'eall_0.Phi.time.p.dot') #time-varying survival, recapture constant
#' eall_0.Phi.dot.p.dist = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist), prefix = 'eall_0.Phi.dot.p.dist') #constant survival, recapture distance-dependent
#' eall_0.Phi.size.p.size = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.size, p=p.size), prefix = 'eall_0.Phi.size.p.size') #capture-size-dependent survival and recapture
#' eall_0.Phi.size.p.dist = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.size, p=p.dist), prefix = 'eall_0.Phi.size.p.dist') #size-dependent survival, distance-dependent recapture
#' eall_0.Phi.size.p.size.plus.dist = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.size, p=p.size.plus.dist), prefix = 'eall_0.Phi.size.p.size.plus.dist') #size-dependent survival, size-and-distance-dependent recapture
#' eall_0.Phi.dot.p.size.plus.dist = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.dot, p=p.size.plus.dist), prefix = 'eall_0.Phi.dot.p.size.plus.dist') #constant survival, size-and-distance-dependent recapture
#' eall_0.Phi.site.p.size.plus.dist = mark(eall_0.processed_site, eall_0.ddl_site, model.parameters=list(Phi=Phi.site, p=p.size.plus.dist), prefix = 'eall_0.Phi.site.p.size.plus.dist') #site-dependent survival, size-and-distance-dependent recapture
#' eall_0.Phi.site.p.dot = mark(eall_0.processed_site, eall_0.ddl_site, model.parameters=list(Phi=Phi.site, p=p.dot), prefix = 'eall_0.Phi.site.p.dot') #site-dependent survival, constant recapture
#' #eall_0.Phi.stage.p.size.plus.dist = mark(eall_0.processed_stage, eall_0.ddl_stage, model.parameters=list(Phi=Phi.stage, p=p.size.plus.dist), prefix = 'eall_0.Phi.stage.p.size.plus.dist') #capture stage-dependent survival, size-and-distance-dependent recapture
#' #eall_0.Phi.stage.p.dot = mark(eall_0.processed_stage, eall_0.ddl_stage, model.parameters=list(Phi=Phi.stage, p=p.dot), prefix = 'eall_0.Phi.stage.p.dot') #capture stage-dependent survival, constant recapture
#' 
#' 
#' # Compare model AICc
#' model_comp_mean = data.frame(model = c('eall_mean.Phi.dot.p.dot','eall_mean.Phi.dot.p.time','eall_mean.Phi.time.p.dot',
#'                                   'eall_mean.Phi.dot.p.dist','eall_mean.Phi.size.p.size','eall_mean.Phi.size.p.dist',
#'                                   'eall_mean.Phi.size.p.size.plus.dist','eall_mean.Phi.dot.p.size.plus.dist',
#'                                   'eall_mean.Phi.site.p.dot','eall_mean.Phi.site.p.size.plus.dist'),
#'                                   #'eall_mean.Phi.stage.p.dot','eall_mean.Phi.stage.p.size.plus.dist'),
#'                         AICc = c(eall_mean.Phi.dot.p.dot$results$AICc, eall_mean.Phi.dot.p.time$results$AICc, eall_mean.Phi.time.p.dot$results$AICc,
#'                                 eall_mean.Phi.dot.p.dist$results$AICc, eall_mean.Phi.size.p.size$results$AICc, eall_mean.Phi.size.p.dist$results$AICc,
#'                                 eall_mean.Phi.size.p.size.plus.dist$results$AICc, eall_mean.Phi.dot.p.size.plus.dist$results$AICc,
#'                                 eall_mean.Phi.site.p.dot$results$AICc, eall_mean.Phi.site.p.dot$results$AICc),
#'                         stringsAsFactors = FALSE)
#'                                 #eall_mean.Phi.stage.p.dot$results$AICc, eall_mean.Phi.stage.p.size.plus.dist$results$AICc))
#' # , 
#' #                         Phi_intercept = c(eall_mean.Phi.dot.p.dot$results$beta$estimate[1], eall_mean.Phi.dot.p.time$results$beta$estimate[1],
#' #                                           eall_mean.Phi.time.p.dot$results$beta$estimate[1], eall_mean.Phi.dot.p.dist$results$beta$estimate[1],
#' #                                           eall_mean.Phi.size.p.size$results$beta$estimate[1], eall_mean.Phi.size.p.dist$results$beta$estimate[1]))
#' # ,
#' #                         p_intercept =  c(eall_mean.Phi.dot.p.dot$results$beta$estimate[2], eall_mean.Phi.dot.p.time$results$beta$estimate[2],
#' #                                          eall_mean.Phi.time.p.dot$results$beta$estimate[1], eall_mean.Phi.dot.p.dist$results$beta$estimate[1],
#' #                                          eall_mean.Phi.size.p.size$results$beta$estimate[1], eall_mean.Phi.size.p.dist$results$beta$estimate[1]))
#' 
#' model_comp_0s = data.frame(model = c('eall_0.Phi.dot.p.dot','eall_0.Phi.dot.p.time','eall_0.Phi.time.p.dot',
#'                                        'eall_0.Phi.dot.p.dist','eall_0.Phi.size.p.size','eall_0.Phi.size.p.dist',
#'                                        'eall_0.Phi.size.p.size.plus.dist','eall_0.Phi.dot.p.size.plus.dist',
#'                                        'eall_0.Phi.site.p.dot','eall_0.Phi.site.p.size.plus.dist'),
#'                                        #'eall_0.Phi.stage.p.dot','eall_0.Phi.stage.p.size.plus.dist'),
#'                              AICc = c(eall_0.Phi.dot.p.dot$results$AICc, eall_0.Phi.dot.p.time$results$AICc, eall_0.Phi.time.p.dot$results$AICc,
#'                                      eall_0.Phi.dot.p.dist$results$AICc, eall_0.Phi.size.p.size$results$AICc, eall_0.Phi.size.p.dist$results$AICc,
#'                                      eall_0.Phi.size.p.size.plus.dist$results$AICc, eall_0.Phi.dot.p.size.plus.dist$results$AICc,
#'                                      eall_0.Phi.site.p.dot$results$AICc, eall_0.Phi.site.p.dot$results$AICc))
#'                                      #eall_0.Phi.stage.p.dot$results$AICc, eall_0.Phi.stage.p.size.plus.dist$results$AICc))
#' 
#' # Arrange by increasing AICc
#' model_comp_mean <- model_comp_mean %>%
#'   mutate(dAICc = min(AICc) - AICc) %>%
#'   arrange(-dAICc)
#' 
#' model_comp_0s <- model_comp_0s %>%
#'   mutate(dAIC = min(AIC) - AIC) %>%
#'   arrange(-dAIC)
#' 
#' # organize output to get ready to make plots
#' ###### Constant survival and recapture prob, no covariates (eall_mean.Phi.dot.p.dot)
#' eall_mean.constant = as.data.frame(eall_mean.Phi.dot.p.dot$results$beta) %>% 
#'   mutate(param = c("Phi","p")) %>% #add a parameters column to make it easier to plot in ggplot
#'   mutate(upper = logit_recip(ucl), lower = logit_recip(lcl), est = logit_recip(estimate)) #do the reciprocal transform on upper and lower confidence limits (need to check what those are - 95? SE? and that transforming them just straight up is the right way to go)
#' 
#' eall_0.constant = as.data.frame(eall_0.Phi.dot.p.dot$results$beta) %>% 
#'   mutate(param = c("Phi","p")) %>% #add a parameters column to make it easier to plot in ggplot
#'   mutate(upper = logit_recip(ucl), lower = logit_recip(lcl), est = logit_recip(estimate)) #do the reciprocal transform on upper and lower confidence limits (need to check what those are - 95? SE? and that transforming them just straight up is the right way to go)
#' 
#' 
#' ###### Survival by size, recapture both by size and distance (eall_mean.Phi.size.p.size.plus.dist)
#' minsize = 1
#' maxsize = 15
#' size.values <- minsize+(0:30)*(maxsize-minsize)/30
#' 
#' mindist = 0
#' # maxdist1 = max(eall$dist2016,eall$dist2017,eall$dist2018) #seems kind of high and probably due to an error
#' maxdist2 = 200 #encompasses highest values in 2016, 2018 
#' # dist.values1 = mindist+(0:30)*(maxdist1-mindist)/30
#' dist.values2 = mindist+(0:30)*(maxdist2-mindist)/30
#' 
#' # Means
#' eall_mean.Phi.size.p.size.plus.dist.results = as.data.frame(eall_mean.Phi.size.p.size.plus.dist$results$beta)
#' 
#' Phibysize_Phi.size.p.size.plus.dist_means <- data.frame(size = size.values) %>%
#'   mutate(Phi_logit = eall_mean.Phi.size.p.size.plus.dist.results$estimate[1] + eall_mean.Phi.size.p.size.plus.dist.results$estimate[2]*size,
#'          Phi_lcl_logit = eall_mean.Phi.size.p.size.plus.dist.results$lcl[1] + eall_mean.Phi.size.p.size.plus.dist.results$lcl[2]*size,
#'          Phi_ucl_logit = eall_mean.Phi.size.p.size.plus.dist.results$ucl[1] + eall_mean.Phi.size.p.size.plus.dist.results$ucl[2]*size,
#'          Phi = logit_recip(Phi_logit),
#'          Phi_lcl = logit_recip(Phi_lcl_logit),
#'          Phi_ucl = logit_recip(Phi_ucl_logit))
#' 
#' pbysize_Phi.size.p.size.plus.dist_means <- data.frame(size = size.values) %>%
#'   mutate(p_logit = eall_mean.Phi.size.p.size.plus.dist.results$estimate[3] + eall_mean.Phi.size.p.size.plus.dist.results$estimate[4]*size,
#'          p_lcl_logit = eall_mean.Phi.size.p.size.plus.dist.results$lcl[3] + eall_mean.Phi.size.p.size.plus.dist.results$lcl[4]*size,
#'          p_ucl_logit = eall_mean.Phi.size.p.size.plus.dist.results$ucl[3] + eall_mean.Phi.size.p.size.plus.dist.results$ucl[4]*size,
#'          p = logit_recip(p_logit),
#'          p_lcl = logit_recip(p_lcl_logit),
#'          p_ucl = logit_recip(p_ucl_logit))
#' 
#' pbydist_Phi.size.p.size.plus.dist_means <- data.frame(dist = dist.values2) %>%
#'   mutate(p_logit = eall_mean.Phi.size.p.size.plus.dist.results$estimate[3] + eall_mean.Phi.size.p.size.plus.dist.results$estimate[5]*dist,
#'          p_lcl_logit = eall_mean.Phi.size.p.size.plus.dist.results$lcl[3] + eall_mean.Phi.size.p.size.plus.dist.results$lcl[5]*dist,
#'          p_ucl_logit = eall_mean.Phi.size.p.size.plus.dist.results$ucl[3] + eall_mean.Phi.size.p.size.plus.dist.results$ucl[5]*dist,
#'          p = logit_recip(p_logit),
#'          p_lcl = logit_recip(p_lcl_logit),
#'          p_ucl = logit_recip(p_ucl_logit))
#' 
#' # Zeros
#' eall_0.Phi.size.p.size.plus.dist.results = as.data.frame(eall_0.Phi.size.p.size.plus.dist$results$beta)
#' 
#' Phibysize_Phi.size.p.size.plus.dist_0s <- data.frame(size = size.values) %>%
#'   mutate(Phi_logit = eall_0.Phi.size.p.size.plus.dist.results$estimate[1] + eall_0.Phi.size.p.size.plus.dist.results$estimate[2]*size,
#'          Phi_lcl_logit = eall_0.Phi.size.p.size.plus.dist.results$lcl[1] + eall_0.Phi.size.p.size.plus.dist.results$lcl[2]*size,
#'          Phi_ucl_logit = eall_0.Phi.size.p.size.plus.dist.results$ucl[1] + eall_0.Phi.size.p.size.plus.dist.results$ucl[2]*size,
#'          Phi = logit_recip(Phi_logit),
#'          Phi_lcl = logit_recip(Phi_lcl_logit),
#'          Phi_ucl = logit_recip(Phi_ucl_logit))
#' 
#' pbysize_Phi.size.p.size.plus.dist_0s <- data.frame(size = size.values) %>%
#'   mutate(p_logit = eall_0.Phi.size.p.size.plus.dist.results$estimate[3] + eall_0.Phi.size.p.size.plus.dist.results$estimate[4]*size,
#'          p_lcl_logit = eall_0.Phi.size.p.size.plus.dist.results$lcl[3] + eall_0.Phi.size.p.size.plus.dist.results$lcl[4]*size,
#'          p_ucl_logit = eall_0.Phi.size.p.size.plus.dist.results$ucl[3] + eall_0.Phi.size.p.size.plus.dist.results$ucl[4]*size,
#'          p = logit_recip(p_logit),
#'          p_lcl = logit_recip(p_lcl_logit),
#'          p_ucl = logit_recip(p_ucl_logit))
#' 
#' pbydist_Phi.size.p.size.plus.dist_0s <- data.frame(dist = dist.values2) %>%
#'   mutate(p_logit = eall_0.Phi.size.p.size.plus.dist.results$estimate[3] + eall_0.Phi.size.p.size.plus.dist.results$estimate[5]*dist,
#'          p_lcl_logit = eall_0.Phi.size.p.size.plus.dist.results$lcl[3] + eall_0.Phi.size.p.size.plus.dist.results$lcl[5]*dist,
#'          p_ucl_logit = eall_0.Phi.size.p.size.plus.dist.results$ucl[3] + eall_0.Phi.size.p.size.plus.dist.results$ucl[5]*dist,
#'          p = logit_recip(p_logit),
#'          p_lcl = logit_recip(p_lcl_logit),
#'          p_ucl = logit_recip(p_ucl_logit))
#' 
#' ####### MOVE CODE UP FROM HERE
#' ###### Survival and recapture both vary by size (eall_mean.Phi.size.p.size)
#' # Means
#' eall_mean.Phi.size.p.size.results = as.data.frame(eall_mean.Phi.size.p.size$results$beta) 
#' Phibysize_means = data.frame(size = size.values) %>%
#'   mutate(Phi_logit = (eall_mean.Phi.size.p.size.results$estimate[1]) + (eall_mean.Phi.size.p.size.results$estimate[2]*size)) %>%
#'   mutate(Phi_lcl_logit = eall_mean.Phi.size.p.size.results$lcl[1] + (eall_mean.Phi.size.p.size.results$lcl[2]*size)) %>%
#'   mutate(Phi_ucl_logit = eall_mean.Phi.size.p.size.results$ucl[1] + (eall_mean.Phi.size.p.size.results$ucl[2]*size)) %>%
#'   mutate(Phi = logit_recip(Phi_logit)) %>%
#'   mutate(Phi_lcl = logit_recip(Phi_lcl_logit)) %>%
#'   mutate(Phi_ucl = logit_recip(Phi_ucl_logit)) %>%
#'   mutate(model = 'Phi.size.p.size',
#'          method = 'means')
#' 
#' pbysize_means = data.frame(size = size.values) %>%
#'   mutate(p_logit = (eall_mean.Phi.size.p.size.results$estimate[3]) + (eall_mean.Phi.size.p.size.results$estimate[4]*size)) %>%
#'   mutate(p_lcl_logit = eall_mean.Phi.size.p.size.results$lcl[3] + (eall_mean.Phi.size.p.size.results$lcl[4]*size)) %>%
#'   mutate(p_ucl_logit = eall_mean.Phi.size.p.size.results$ucl[3] + (eall_mean.Phi.size.p.size.results$ucl[4]*size)) %>%
#'   mutate(p = logit_recip(p_logit)) %>%
#'   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#'   mutate(p_ucl = logit_recip(p_ucl_logit)) %>%
#'   mutate(model = 'Phi.size.p.size',
#'          method = 'means')
#' 
#' # 0s
#' eall_0.Phi.size.p.size.results = as.data.frame(eall_0.Phi.size.p.size$results$beta) 
#' Phibysize_0 = data.frame(size = size.values) %>%
#'   mutate(Phi_logit = (eall_0.Phi.size.p.size.results$estimate[1]) + (eall_0.Phi.size.p.size.results$estimate[2]*size)) %>%
#'   mutate(Phi_lcl_logit = eall_0.Phi.size.p.size.results$lcl[1] + (eall_0.Phi.size.p.size.results$lcl[2]*size)) %>%
#'   mutate(Phi_ucl_logit = eall_0.Phi.size.p.size.results$ucl[1] + (eall_0.Phi.size.p.size.results$ucl[2]*size)) %>%
#'   mutate(Phi = logit_recip(Phi_logit)) %>%
#'   mutate(Phi_lcl = logit_recip(Phi_lcl_logit)) %>%
#'   mutate(Phi_ucl = logit_recip(Phi_ucl_logit)) %>%
#'   mutate(model = 'Phi.size.p.size',
#'          method = 'zeros')
#' 
#' pbysize_0 = data.frame(size = size.values) %>%
#'   mutate(p_logit = (eall_0.Phi.size.p.size.results$estimate[3]) + (eall_0.Phi.size.p.size.results$estimate[4]*size)) %>%
#'   mutate(p_lcl_logit = eall_0.Phi.size.p.size.results$lcl[3] + (eall_0.Phi.size.p.size.results$lcl[4]*size)) %>%
#'   mutate(p_ucl_logit = eall_0.Phi.size.p.size.results$ucl[3] + (eall_0.Phi.size.p.size.results$ucl[4]*size)) %>%
#'   mutate(p = logit_recip(p_logit)) %>%
#'   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#'   mutate(p_ucl = logit_recip(p_ucl_logit)) %>%
#'   mutate(model = 'Phi.size.p.size',
#'          method = 'zeros')
#' 
#' Phibysize <- rbind(Phibysize_0, Phibysize_means)
#' pbysize <- rbind(pbysize_0, pbysize_means)
#' 
#' ###### Survival by size, recapture by distance (eall_mean.Phi.size.p.dist)
#' mindist = 1
#' maxsize = 15
#' size.values <- minsize+(0:30)*(maxsize-minsize)/30
#' mindist = 0
#' # maxdist1 = max(eall$dist2016,eall$dist2017,eall$dist2018) #seems kind of high and probably due to an error
#' maxdist2 = 200 #encompasses highest values in 2016, 2018 
#' # dist.values1 = mindist+(0:30)*(maxdist1-mindist)/30
#' dist.values2 = mindist+(0:30)*(maxdist2-mindist)/30
#' 
#' eall_mean.Phi.size.p.dist.results = as.data.frame(eall_mean.Phi.size.p.dist$results$beta) 
#' 
#' Phibysize_Phisizepdist = data.frame(size = size.values) %>%
#'   mutate(Phi_logit = (eall_mean.Phi.size.p.dist.results$estimate[1]) + (eall_mean.Phi.size.p.dist.results$estimate[2]*size)) %>%
#'   mutate(Phi_lcl_logit = eall_mean.Phi.size.p.dist.results$lcl[1] + (eall_mean.Phi.size.p.dist.results$lcl[2]*size)) %>%
#'   mutate(Phi_ucl_logit = eall_mean.Phi.size.p.dist.results$ucl[1] + (eall_mean.Phi.size.p.dist.results$ucl[2]*size)) %>%
#'   mutate(Phi = logit_recip(Phi_logit)) %>%
#'   mutate(Phi_lcl = logit_recip(Phi_lcl_logit)) %>%
#'   mutate(Phi_ucl = logit_recip(Phi_ucl_logit)) 
#' 
#' pbydist_Phisizepdist = data.frame(dist = dist.values2) %>%
#'   mutate(p_logit = (eall_mean.Phi.size.p.dist.results$estimate[3]) + (eall_mean.Phi.size.p.dist.results$estimate[4]*dist)) %>%
#'   mutate(p_lcl_logit = eall_mean.Phi.size.p.dist.results$lcl[3] + (eall_mean.Phi.size.p.dist.results$lcl[4]*dist)) %>%
#'   mutate(p_ucl_logit = eall_mean.Phi.size.p.dist.results$ucl[3] + (eall_mean.Phi.size.p.dist.results$ucl[4]*dist)) %>%
#'   mutate(p = logit_recip(p_logit)) %>%
#'   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#'   mutate(p_ucl = logit_recip(p_ucl_logit)) 
#' 
#' ###### Survival constant, recapture varies by time (eall_mean.Phi.dot.p.time)
#' eall_mean.Phi.dot.p.time = as.data.frame(eall_mean.Phi.dot.p.time$results$beta) #%>% 
#'   mutate(param = c("Phi","p")) %>% #add a parameters column to make it easier to plot in ggplot
#'   mutate(upper = logit_recip(ucl), lower = logit_recip(lcl), est = logit_recip(estimate)) #do the reciprocal transform on upper and lower confidence limits (need to check what those are - 95? SE? and that transforming them just straight up is the right way to go)
#' 
#' ###### Survival constant, p varies by distance (eall_mean.Phi.dot.p.dist)
#' mindist = min(eall$dist2016,eall$dist2017,eall$dist2018) #this is 0, obv...
#' maxdist1 = max(eall$dist2016,eall$dist2017,eall$dist2018) #seems kind of high - not an error, though, it's the fish that moved from Wangag to Hicgop South!
#' maxdist2 = 200 #encompasses highest values in 2016, 2018 
#' #tail(sort(eall$dist2016),10)
#' #tail(sort(eall$dist2017),10) #really 2017 that has the high distances...
#' #tail(sort(eall$dist2018),10)
#' dist.values1 = mindist+(0:30)*(maxdist1-mindist)/30
#' dist.values2 = mindist+(0:30)*(maxdist2-mindist)/30
#' 
#' # Not getting this to work - it doesn't change with distance... maybe/probably I'm specifying the index incorrectly?
#' #pbydist1 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values),indices=c(7))
#' #pbydist2 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values2),indices=c(1))
#' #pbydist3 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values3),indices=c(1))
#' 
#' # Set up dataframe 
#' eall.Phi.dot.p.dist.results = as.data.frame(eall.Phi.dot.p.dist$results$beta) 
#' 
#' # Large range of distances (up to 5000m)
#' pbydist_longrange= data.frame(dist = dist.values1) %>%
#'   mutate(p_logit = (eall.Phi.dot.p.dist.results$estimate[2]) + (eall.Phi.dot.p.dist.results$estimate[3]*dist)) %>%
#'   mutate(p_lcl_logit = eall.Phi.dot.p.dist.results$lcl[2] + (eall.Phi.dot.p.dist.results$lcl[3]*dist)) %>%
#'   mutate(p_ucl_logit = eall.Phi.dot.p.dist.results$ucl[2] + (eall.Phi.dot.p.dist.results$ucl[3]*dist)) %>%
#'   mutate(p = logit_recip(p_logit)) %>%
#'   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#'   mutate(p_ucl = logit_recip(p_ucl_logit)) 
#' 
#' # Short range of distances (up to 200m)
#' pbydist_shortrange = data.frame(dist = dist.values2) %>%
#'   mutate(p_logit = (eall.Phi.dot.p.dist.results$estimate[2]) + (eall.Phi.dot.p.dist.results$estimate[3]*dist)) %>%
#'   mutate(p_lcl_logit = eall.Phi.dot.p.dist.results$lcl[2] + (eall.Phi.dot.p.dist.results$lcl[3]*dist)) %>%
#'   mutate(p_ucl_logit = eall.Phi.dot.p.dist.results$ucl[2] + (eall.Phi.dot.p.dist.results$ucl[3]*dist)) %>%
#'   mutate(p = logit_recip(p_logit)) %>%
#'   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#'   mutate(p_ucl = logit_recip(p_ucl_logit)) 
#' 
#' ###### Survival varies by tagging size, p constant (eall.Phi.size.p.dot)
#' minsize = min(eall$tag_size) #right now, this is 1.6 (b/c of a typo.... Michelle fixed - need to re-call data from db and run distance calcs again)
#' maxsize = max(eall$tag_size)
#' size.values = minsize+(0:30)*(maxsize-minsize)/30
#' Phibysize = covariate.predictions(eall.Phi.size.p.dot,data=data.frame(tag_size=size.values),indices=c(1)) #should do this the way I did distance too, make sure get the same thing...
#' 
#' ###### Survival varies by tail color, p constant (eall.Phi.color.p.dot)
#' eall.color = as.data.frame(eall.Phi.color.p.dot$results$real) %>%
#'   mutate(param = c(rep("Phi",6),"p")) %>%
#'   mutate(color = c("BW","O","W","Y","YP","YR","all")) %>%
#'   mutate(row = seq(1,7,1)) #makes plotting easier
#' 
#' ###### Survival varies by tagging size, p varies by distance (eall.Phi.size.p.dist)
#' minsize = min(eall$tag_size) #right now, this is 2.6 (prob still a typo? Michelle fixed the 1.6 one but are there more that are too small?)
#' maxsize = max(eall$tag_size)
#' size.values = minsize+(0:30)*(maxsize-minsize)/30
#' Phibysize_Phisizepdist = covariate.predictions(eall.Phi.size.p.dist,data=data.frame(tag_size=size.values),indices=c(1)) #should do this the way I did distance too, make sure get the same thing...
#' 
#' mindist = min(eall$dist2016,eall$dist2017,eall$dist2018) #this is 0, obv...
#' maxdist1 = max(eall$dist2016,eall$dist2017,eall$dist2018) #seems kind of high and probably due to an error
#' maxdist2 = 200 #encompasses highest values in 2016, 2018 
#' dist.values1 = mindist+(0:30)*(maxdist1-mindist)/30
#' dist.values2 = mindist+(0:30)*(maxdist2-mindist)/30
#' 
#' eall.Phi.size.p.dist.results = as.data.frame(eall.Phi.size.p.dist$results$beta) 
#' 
#' # Size
#' pbydistPhisize_size= data.frame(tag_size = size.values) %>%
#'   mutate(Phi_logit = (eall.Phi.size.p.dist.results$estimate[1]) + (eall.Phi.size.p.dist.results$estimate[2]*tag_size)) %>%
#'   mutate(Phi_lcl_logit = eall.Phi.size.p.dist.results$lcl[1] + (eall.Phi.size.p.dist.results$lcl[2]*tag_size)) %>%
#'   mutate(Phi_ucl_logit = eall.Phi.size.p.dist.results$ucl[1] + (eall.Phi.size.p.dist.results$ucl[2]*tag_size)) %>%
#'   mutate(Phi = logit_recip(Phi_logit)) %>%
#'   mutate(Phi_lcl = logit_recip(Phi_lcl_logit)) %>%
#'   mutate(Phi_ucl = logit_recip(Phi_ucl_logit)) 
#' 
#' # Large range of distances (up to 5000m)
#' pbydistPhisize_longrange= data.frame(dist = dist.values1) %>%
#'   mutate(p_logit = (eall.Phi.size.p.dist.results$estimate[3]) + (eall.Phi.size.p.dist.results$estimate[4]*dist)) %>%
#'   mutate(p_lcl_logit = eall.Phi.size.p.dist.results$lcl[3] + (eall.Phi.size.p.dist.results$lcl[4]*dist)) %>%
#'   mutate(p_ucl_logit = eall.Phi.size.p.dist.results$ucl[3] + (eall.Phi.size.p.dist.results$ucl[4]*dist)) %>%
#'   mutate(p = logit_recip(p_logit)) %>%
#'   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#'   mutate(p_ucl = logit_recip(p_ucl_logit)) 
#' 
#' # Short range of distances (up to 200m)
#' pbydistPhisize_shortrange = data.frame(dist = dist.values2) %>%
#'   mutate(p_logit = (eall.Phi.size.p.dist.results$estimate[3]) + (eall.Phi.size.p.dist.results$estimate[4]*dist)) %>%
#'   mutate(p_lcl_logit = eall.Phi.size.p.dist.results$lcl[3] + (eall.Phi.size.p.dist.results$lcl[4]*dist)) %>%
#'   mutate(p_ucl_logit = eall.Phi.size.p.dist.results$ucl[3] + (eall.Phi.size.p.dist.results$ucl[4]*dist)) %>%
#'   mutate(p = logit_recip(p_logit)) %>%
#'   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#'   mutate(p_ucl = logit_recip(p_ucl_logit)) 
#' 
#' ###### Survival varies by site, p constant (eall.Phi.site.p.dot)
#' eall.Phisite = as.data.frame(eall.Phi.site.p.dot$results$real) %>%
#'   mutate(param = c(rep("Phi",15),"p")) %>%
#'   mutate(site = c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
#'                  "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag","all")) %>%
#'   mutate(row = seq(1,16,1))
#' 
#' ###### Survival varies by site, p by distance (eall.Phi.site.p.dist) - doesn't look like this one parsed distance correctly?
#' eall.Phisitepdist_Phi = as.data.frame(eall.Phi.site.p.dist$results$real[1:15,]) %>%
#'   mutate(param = c(rep("Phi",15))) %>%
#'   mutate(site = c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
#'                   "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag")) %>%
#'   mutate(row = seq(1,15,1))
#' 
#' ###### Survival varies by site, p varies by site (eall.Phi.site.p.site)
#' eall.site = as.data.frame(eall.Phi.site.p.site$results$beta) %>%
#'   mutate(param = c(rep("Phi",15),rep("p",15))) %>% #add a parameter column to make it easier to plot in ggplot
#'   mutate(site = c(rep(c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
#'                         "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag"),2))) #add site column (Cabatoan is intercept)
#' 
#' eall.site.real = as.data.frame(eall.Phi.site.p.site$results$real) %>%
#'   mutate(param = c(rep("Phi",15),rep("p",15))) %>% #add a parameter column to make it easier to plot in ggplot
#'   mutate(site = c(rep(c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
#'                         "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag"),2))) %>% #add site column (Cabatoan is intercept)
#'   mutate(row = seq(1,30,1)) #add row numbers to make plotting easier
#' 
#' #################### Plots: ####################
#' ##### PLOTS WITH NEW GEN DATA
#' 
#' ###### Survival varies by size, p varies by size and distance - means instead of NA (for pre-capture records) (eall_mean.Phi.size.p.size.plus.dist)
#' ## Recap plot - effect of size (p)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_mean_Phisize_psizeanddist_newGen_p_size.pdf"))
#' ggplot(data = pbysize_Phi.size.p.size.plus.dist_means, aes(size, p)) +
#'   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("size of fish (cm)") + ylab("p estimate") + ggtitle("p size effect (Phi size, p size, dist)") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ## Recap plot - effect of dist (p) - this looks weird... not the shape I would have expected
#' pdf(file = here("Plots/PhiandpEstimates", "eall_mean_Phisize_psizeanddist_newGen_p_dist.pdf"))
#' ggplot(data = pbydist_Phi.size.p.size.plus.dist_means, aes(dist, p)) +
#'   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("dist to anem (m)") + ylab("p estimate") + ggtitle("p dist effect (Phi size, p size, dist)") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ## Survival plot - effect of size (Phi)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_mean_Phisize_psizeanddist_newGen_Phi_size.pdf"))
#' ggplot(data = Phibysize_Phi.size.p.size.plus.dist_means, aes(size, Phi)) +
#'   geom_ribbon(aes(ymin=Phi_lcl,ymax=Phi_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("size of fish (cm)") + ylab("Phi estimate") + ggtitle("Phi size effect (Phi size, p size, dist)") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ###### Constant survival and recapture prob, no covariates (eall_mean.Phi.dot.p.dot)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phiandp_constant_newGen.pdf"))
#' ggplot(data = eall_mean.constant, aes(param, est, color=param)) +
#'   geom_pointrange(aes(ymin=lower, ymax=upper), size=1) +
#'   xlab("parameter") + ylab("estimate") + ggtitle("Constant p and Phi (eall.Phi.dot.p.dot), 2016-2018 gen") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ###### Survival and p varies by size (eall_mean.Phi.size.p.size) 
#' ## Recap plot (p)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_mean_Phisize_psize_newGen_p.pdf"))
#' ggplot(data = pbysize, aes(size, p)) +
#'   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("size of fish (cm)") + ylab("p estimate") + ggtitle("Phi and p by size (eall.Phi.size.p.size)), 2016-2018 gen") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ## Survival plot (Phi)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_mean_Phisize_psize_newGen_Phi.pdf"))
#' ggplot(data = Phibysize, aes(size, Phi)) +
#'   geom_ribbon(aes(ymin=Phi_lcl,ymax=Phi_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("size of fish (cm)") + ylab("Phi estimate") + ggtitle("Phi and p by size (eall.Phi.size.p.size)), 2016-2018 gen") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ###### Survival by size, p by distance (eall_mean.Phi.size.p.dist) 
#' ## Recap plot (p)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_mean_Phisize_pdist_newGen_p.pdf"))
#' ggplot(data = pbydist_Phisizepdist, aes(dist, p)) +
#'   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("dist to anem (m)") + ylab("p estimate") + ggtitle("Phi size, p dist (eall.Phi.size.p.dist)), 2016-2018 gen") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ## Survival plot (Phi)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_mean_Phisize_pdist_newGen_Phi.pdf"))
#' ggplot(data = Phibysize_Phisizepdist, aes(size, Phi)) +
#'   geom_ribbon(aes(ymin=Phi_lcl,ymax=Phi_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("size of fish (cm)") + ylab("Phi estimate") + ggtitle("Phi size, p dist (eall.Phi.size.p.dist)), 2016-2018 gen") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' #################################### HAVEN'T RE-MADE PLOTS BELOW WITH NEW GENETIC DATA INCLUDED #####################################
#' ##### Encounter history histograms
#' # Histogram of encounter histories including genetic and tag recaptures (genetic only through 2015)
#' pdf(file = here::here('Plots/PhiandpEstimates', 'HistofEncounterHistories_TimesCaught.pdf'))
#' ggplot(data = encounters_all_table, aes(x=encounter.hist, y=Freq, fill=factor(times_caught))) +
#'   geom_bar(stat='identity') +
#'   ggtitle('Histogram of encounter histories (genetic + tag)') +
#'   xlab('encounter history') + ylab('# fish') +
#'   theme_bw() +
#'   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
#' dev.off()
#' 
#' # Histogram of encounter histories for tag-only fish (2015-2018)
#' pdf(file = here::here('Plots/PhiandpEstimates', 'HistofTagEncounterHistories_TimesCaught.pdf'))
#' ggplot(data = encounters_tag_table, aes(x=encounter.hist, y=Freq, fill=factor(times_caught))) +
#'   geom_bar(stat='identity') +
#'   ggtitle('Histogram of encounter histories (tag only)') +
#'   xlab('encounter history') + ylab('# fish') +
#'   theme_bw() +
#'   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#' dev.off()
#' 
#' # make some plots
#' ###### Constant survival and recapture prob, no covariates (eall.Phi.dot.p.dot)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phiandp_constant.pdf"))
#' ggplot(data = eall_mean.constant, aes(param, est, color=param)) +
#'   geom_pointrange(aes(ymin=lower, ymax=upper), size=1) +
#'   xlab("parameter") + ylab("estimate") + ggtitle("Constant p and Phi (eall.Phi.dot.p.dot)") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ###### Survival constant, p varies by distance (eall.Phi.dot.p.dist)
#' # all distances included, up to 5000m from anem (even though large ones are probably an error...)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phidot_pdist_alldists.pdf"))
#' ggplot(data = pbydist_longrange, aes(dist, p)) +
#'   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("Constant Phi, p by distance (eall.Phi.dot.p.dist))") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' # most distances included (all from 2016, 2018, not some of the higher ones from 2017)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phidot_pdist_shortdists.pdf"))
#' ggplot(data = pbydist_shortrange, aes(dist, p)) +
#'   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("Constant Phi, p by distance (eall.Phi.dot.p.dist))") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ###### Survival varies by tagging size, p constant (eall.Phi.size.p.dot)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phisize.pdf"))
#' ggplot(data = Phibysize$estimates, aes(covdata, estimate)) +
#'   geom_ribbon(aes(ymin=lcl,ymax=ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("size at tagging (cm)") + ylab("Phi estimate") + ggtitle("Phi varies by tagging size, p constant (eall.Phi.size.p.dot)") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' # same as above but formatted for MPE poster
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phisize_MPEposter.pdf"))
#' ggplot(data = Phibysize$estimates, aes(covdata, estimate)) +
#'   geom_ribbon(aes(ymin=lcl,ymax=ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black", size=3) +
#'   xlab("size at tagging (cm)") + ylab("survival estimate") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw() +
#'   theme(text = element_text(size=40)) +
#'   theme(axis.text.x = element_text(size=30)) + theme(axis.text.y = element_text(size=30)) 
#' dev.off()
#' 
#' 
#' ###### Survival varies by tail color, p constant (eall.Phi.color.p.dot)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phicolor.pdf"))
#' ggplot(data = eall.color, aes(row, estimate, color=color, shape=param)) +
#'   geom_pointrange(aes(ymin=lcl, ymax=ucl), size=1) +
#'   xlab("parameter") + ylab("estimate") + ggtitle("Phi by tail color, p constant (eall.Phi.color.p.dot)") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ###### Survival varies by tagging size, p varies by distance (eall.Phi.size.p.dist)
#' # all distances included, up to 5000m from anem (even though large ones are probably an error...)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phisize_pdist_alldists.pdf"))
#' ggplot(data = pbydistPhisize_longrange, aes(dist, p)) +
#'   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("p: Phi by size, p by distance (eall.Phi.size.p.dist))") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' # most distances included (all from 2016, 2018, not some of the higher ones from 2017)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phisize_pdist_shortdists.pdf"))
#' ggplot(data = pbydistPhisize_shortrange, aes(dist, p)) +
#'   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("p: Phi by size, p by distance (eall.Phi.size.p.dist))") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' # Phi by tag_size
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phisize_pdist_Phi.pdf"))
#' ggplot(data = pbydistPhisize_size, aes(tag_size, Phi)) +
#'   geom_ribbon(aes(ymin=Phi_lcl,ymax=Phi_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black") +
#'   xlab("size at tagging (cm)") + ylab("Phi estimate") + ggtitle("Phi: Phi by size, p by distance (eall.Phi.size.p.dist))") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' # Phi by tag_size, version for MPE poster
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phisize_pdist_Phi_MPEposter.pdf"))
#' ggplot(data = pbydistPhisize_size, aes(tag_size, Phi)) +
#'   geom_ribbon(aes(ymin=Phi_lcl,ymax=Phi_ucl),color="light blue",fill="light blue") +
#'   geom_line(color="black", size = 3) +
#'   xlab("size at tagging (cm)") + ylab("estimated survival") + ggtitle("Size effect on survival") +
#'   scale_y_continuous(limits = c(0, 1)) +
#'   theme(text = element_text(size=40)) +
#'   theme(axis.text.x = element_text(size=30)) + theme(axis.text.y = element_text(size=30)) +
#'   theme_bw()
#' dev.off()
#' 
#' ###### Survival varies by site, p constant (eall.Phi.site.p.dot)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phisite.pdf"))
#' ggplot(data = eall.Phisite, aes(row, estimate, color=site, shape=param)) +
#'   geom_pointrange(aes(ymin=lcl, ymax=ucl), size=1) +
#'   xlab("parameter") + ylab("Phi estimate") + ggtitle("Phi by site, p constant (eall.Phi.site.p.dot)") +
#'   #scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ###### Survival varies by site, p by distance (eall.Phi.site.p.dist) (doesn't look like this one did distance right - has a couple of params for it...)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phisitepdist_Phi.pdf"))
#' ggplot(data = eall.Phisitepdist_Phi, aes(row, estimate, color=site)) +
#'   geom_pointrange(aes(ymin=lcl, ymax=ucl), size=1) +
#'   xlab("parameter") + ylab("Phi estimate") + ggtitle("Phi: Phi by site, p by dist (eall.Phi.site.p.dist)") +
#'   #scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ###### Survival varies by site, p varies by site (eall.Phi.site.p.site)
#' pdf(file = here("Plots/PhiandpEstimates", "eall_Phiandp_site.pdf"))
#' ggplot(data = eall.site.real, aes(row, estimate, color=site, shape=param)) +
#'   geom_pointrange(aes(ymin=lcl, ymax=ucl), size=1) +
#'   xlab("parameter") + ylab("estimate") + ggtitle("Phi and p by site (eall.Phi.site.p.site)") +
#'   #scale_y_continuous(limits = c(0, 1)) +
#'   theme_bw()
#' dev.off()
#' 
#' ###### Comparing Phi and p across models
#' 
#' #################### Saving files and output ####################
#' saveRDS(encounters_all_means, file=here::here("Data/Script_outputs", "encounters_all_means.RData"))
#' saveRDS(eall_mean.Phi.size.p.size.plus.dist.results, file=here::here("Data/Script_outputs","eall_mean_Phi_size_p_size_plus_dist.RData"))
#' saveRDS(model_comp_mean, file=here::here("Data/Script_outputs","model_comp_mean.RData"))
#' 
#' save(min_survey_dist_to_anems, file=here::here("Data/Script_outputs", "min_survey_dist_to_anems.RData"))
#' 
#' save(encounters_all, file=here("Data", "encounters_all.RData"))
#' save(encounters_all_means, file=here("Data", "encounters_all_means.RData"))
#' save(encounters_all_0s, file=here("Data", "encounters_all_0s.RData"))
#' save(allfish, file=here("Data", "allfish.RData"))
#' 
#' load(file=here('Data','encounters_all_means.RData'))
#' load(file=here('Data','encounters_all_0s.RData'))
#' 
#' save(eall.Phi.size.p.dist.results, file=here("Data","eall_Phi_size_p_dist_results.RData"))
#' 
#' save(eall_mean.Phi.size.p.size.plus.dist, file=here('Data','eall_mean_Phi_size_p_size_plus_dist.RData'))
#' 
#' save(encounters_all_fulldata, file=here("Data", "encounters_all_fulldata.RData")) #for KC Michigan workshop
#' 
#' 
#' #### TO DOS
#' # Do runs with size for Phi, site for Phi?, dist for p
#' # Make plots to look at those
#' # Compare models (AIC? goodness of fit?) and summarize in some way
#' # Check the outlier distances, re-run with any errors fixed there and a fresh pull from the db
#' # Move on to metrics - write-up and placeholder script for those!
#' 
#' 
#' #################### OLD CODE - BOTH MARK AND NOT ####################
#' # eall_mean.Phi.size.plus.stage.p.dot = mark(eall_mean.processed_stage, eall_mean.ddl_stage, model.parameters=list(Phi=Phi.size.plus.stage, p=p.dot), prefix = 'eall_mean.Phi.size.plus.stage.p.dot') #size and stage dependent survival, constant recapture
#' # eall_mean.Phi.size.plus.stage.p.size.plus.stage = mark(eall_mean.processed_stage, eall_mean.ddl_stage, model.parameters=list(Phi=Phi.size.plus.stage, p=p.size.plus.stage), prefix = 'eall_mean_Phi.size.plus.stage.p.size.plus.stage') #size and stage dependent survival and recap
#' # eall_mean.Phi.size.plus.stage.p.size = mark(eall_mean.processed_stage, eall_mean.ddl_stage, model.parameters=list(Phi=Phi.size.plus.stage, p=p.size), prefix = 'eall_mean.Phi.size.plus.stage.p.size') #size and stage dependent survival, size-dependent recapture
#' # eall_mean.Phi.size.plus.stage.p.dist = mark(eall_mean.processed_stage, eall_mean.ddl_stage, model.parameters=list(Phi=Phi.size.plus.stage, p=p.dist), prefix = 'eall_mean.Phi.size.plus.stage.p.dist') #size and stage dependent survival, distance-dependent recapture
#' # eall_mean.Phi.size.plus.stage.p.stage = mark(eall_mean.processed_stage, eall_mean.ddl_stage, model.parameters=list(Phi=Phi.size.plus.stage, p=p.stage), prefix = 'eall_mean.Phi.size.plus.stage.p.stage') #size and stage dependent survival, stage-dependent recapture
#' # eall_mean.Phi.stage.p.dot = mark(eall_mean.processed_stage, eall_mean.ddl_stage, model.parameters=list(Phi=Phi.stage, p=p.dot), prefix = 'eall_mean.Phi.stage.p.dot') #stage-dependent survival, constant recapture
#' # eall_mean.Phi.stage.p.dist = mark(eall_mean.processed_stage, eall_mean.ddl_stage, model.parameters=list(Phi=Phi.stage, p=p.dist), prefix = 'eall_mean.Phi.stage.p.dist') #stage-dependent survival, distance-dependent recapture
#' # eall_mean.Phi.size.plus.color.p.dot = mark(eall_mean.processed_color, eall_mean.ddl_color, model.parameters=list(Phi=Phi.size.plus.color, p=p.dot), prefix = 'eall_mean.Phi.size.plus.color.p.dot') #size and color dependent survival, constant recapture
#' # eall_mean.Phi.size.plus.color.p.dist = mark(eall_mean.processed_color, eall_mean.ddl_color, model.parameters=list(Phi=Phi.size.plus.color, p=p.dist), prefix = 'eall_mean.Phi.size.plus.color.p.dist') #size and color dependent survival, distance-depedent recapture
#' # eall_mean.Phi.color.p.dot = mark(eall_mean.processed_color, eall_mean.ddl_color, model.parameters=list(Phi=Phi.color, p=p.dot), prefix = 'eall_mean.Phi.color.p.dot') #color dependent survival, constant recapture
#' # eall_mean.Phi.color.p.dist = mark(eall_mean.processed_color, eall_mean.ddl_color, model.parameters=list(Phi=Phi.color, p=p.dist), prefix = 'eall_mean.Phi.color.p.dist') #color-dependent survival, distance-dependent recapture
#' 
#' # eall_0.Phi.dot.p.dist = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist), prefix = 'eall_0.Phi.dot.p.dist') #constant survival, recapture distance-dependent
#' # eall_0.Phi.size.p.dist = mark(eall_0.processed, eall_0.ddl, model.parameters=list(Phi=Phi.size, p=p.dist), prefix = 'eall_0.Phi.size.p.dist') #size-dependent survival, distance-dependent recapture
#' # eall_0.Phi.size.plus.stage.p.dist = mark(eall_0.processed_stage, eall_0.ddl_stage, model.parameters=list(Phi=Phi.size.plus.stage, p=p.dist), prefix = 'eall_0.Phi.size.plus.stage.p.dist') #size and stage dependent survival, distance-dependent recapture
#' # eall_0.Phi.size.plus.color.p.dist = mark(eall_0.processed_color, eall_0.ddl_color, model.parameters=list(Phi=Phi.size.plus.color, p=p.dist), prefix = 'eall_0.Phi.size.plus.color.p.dist') #size and color dependent survival, distance-depedent recapture
#' # eall_0.Phi.color.p.dist = mark(eall_0.processed_color, eall_0.ddl_color, model.parameters=list(Phi=Phi.color, p=p.dist), prefix = 'eall_0.Phi.color.p.dist') #color-dependent survival, distance-dependent recapture
#' 
#' 
#' # ######## Try with a function
#' # fit.eall.models <- function() {
#' #   
#' #   # Process data
#' #   eall_mean.processed = process.data(eall_mean, model = "CJS", groups = "cap_stage")
#' #   eall_mean.ddl = make.design.data(eall_mean.processed)
#' #   
#' #   # set models for Phi
#' #   Phi.dot = list(formula=~1, link="logit")
#' #   Phi.time = list(formula=~time, link="logit")
#' #   Phi.size = list(formula=~cap_size, link="logit")
#' #   Phi.stage = list(formula=~cap_stage, link="logit")
#' #   Phi.size.plus.stage = list(formula=~cap_size+cap_stage, link='logit')
#' #   
#' #   # set models for p
#' #   p.dot = list(formula=~1, link="logit")
#' #   p.time = list(formula=~time, link="logit")
#' #   p.dist = list(formula=~dist, link="logit")
#' #   p.size = list(formula=~cap_size, link="logit")
#' #   p.stage = list(formula=~cap_stage, link="logit")
#' #   p.size.plus.stage = list(formula=~cap_size+cap_stage, link="logit")
#' #   
#' #   # Create model list
#' #   cml = create.model.list("CJS")
#' #   
#' #   # Run and return marklist of models
#' #   return(mark.wrapper(cml, data = eall_mean.processed, ddl = eall_mean.ddl))
#' # }
#' # 
#' # # Fit models
#' # eall_mean.models <- fit.eall.models()
#' 
#' # eall.Phi.dot.p.time = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.dot, p=p.time)) #has issues running, can't estimate p in 2018
#' # eall.Phi.dot.p.dist = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist)) #think this ran!!!
#' # eall.Phi.size.p.dot = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.size, p=p.dot))
#' # eall.Phi.size.p.site = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.size, p=p.site))
#' # eall.Phi.size.p.dist = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.size, p=p.dist))
#' # eall.Phi.site.p.site = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.site, p=p.site)) #this one has issues converging for some of the sites
#' # eall.Phi.color.p.dot = mark(eall.processed2, eall.ddl2, model.parameters=list(Phi=Phi.color, p=p.dot))
#' # eall.Phi.site.p.dot = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.site, p=p.dot)) #hmmm, something seems funny with this one...
#' # eall.Phi.site.p.dist = mark(eall.processed, eall.ddl, model.parameters=list(Phi=Phi.site, p=p.dist))
#' # 
#' # #edist.Phi.dot.p.dist.plus.time = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.dot, p=p.time.plus.dist))
#' # edist.Phi.time.p.dot = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.dot))
#' # edist.Phi.time.p.time = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.time))
#' # edist.Phi.time.p.dist = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.dist))
#' # edist.Phi.time.p.time.plus.dist = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.time.plus.dist))
#' 
#' 
#' # model_comp = data.frame(model = c('eall_mean.Phi.dot.p.dot','eall_mean.Phi.dot.p.time','eall_mean.Phi.time.p.dot',
#' #                                   'eall_mean.Phi.dot.p.dist','eall_mean.Phi.size.p.size','eall_mean.Phi.size.p.dist',
#' #                                   'eall_mean.Phi.size.plus.stage.p.dot','eall_mean.Phi.size.plus.stage.p.size.plus.stage',
#' #                                   'eall_mean.Phi.size.plus.stage.p.size','eall_mean.Phi.size.plus.stage.p.dist',
#' #                                   'eall_mean.Phi.size.plus.stage.p.stage','eall_mean.Phi.stage.p.dot', 'eall_mean.Phi.stage.p.dist',
#' #                                   'eall_mean.Phi.size.plus.color.p.dot','eall_mean.Phi.size.plus.color.p.dist',
#' #                                   'eall_mean.Phi.color.p.dot','eall_mean.Phi.color.p.dist'),
#' #                         AIC = c(eall_mean.Phi.dot.p.dot$results$AICc, eall_mean.Phi.dot.p.time$results$AICc, eall_mean.Phi.time.p.dot$results$AICc,
#' #                                 eall_mean.Phi.dot.p.dist$results$AICc, eall_mean.Phi.size.p.size$results$AICc, eall_mean.Phi.size.p.dist$results$AICc,
#' #                                 eall_mean.Phi.size.plus.stage.p.dot$results$AICc, eall_mean.Phi.size.plus.stage.p.size.plus.stage$results$AICc,
#' #                                 eall_mean.Phi.size.plus.stage.p.size$results$AICc, eall_mean.Phi.size.plus.stage.p.dist$results$AICc,
#' #                                 eall_mean.Phi.size.plus.stage.p.stage$results$AICc, eall_mean.Phi.stage.p.dot$results$AICc,
#' #                                 eall_mean.Phi.stage.p.dist$results$AICc, eall_mean.Phi.size.plus.color.p.dot$results$AICc,
#' #                                 eall_mean.Phi.size.plus.color.p.dist$results$AICc, eall_mean.Phi.color.p.dot$results$AICc, eall_mean.Phi.color.p.dist$results$AICc))
#' 
#' ### Old from before 2016-2018 genetic data was in, trying to remove fish that had no chance of being recaught
#' # # remove fish that are were clipped but not tagged in 2015 (since no chance for us to know if we have recaptured them yet) - fish_id is "genXXXXX", encounter hist is 0001000
#' # encounters_all_0s <- encounters_all_0s %>%
#' #   mutate(rem_2015 = if_else(ch == '0001000' & substring(fish_id,1,3) == 'gen', 1, 0)) #only get flagged for removal if only seen in 2015 and identified by gen_id (spot-checking, looks like ~ 100 juvenile fish, which makes sense (actually filters out 107 fish))
#' # encounters_all_0s <- encounters_all_0s %>%
#' #   filter(rem_2015 == 0)
#' # 
#' # encounters_all_means <- encounters_all_means %>%
#' #   mutate(rem_2015 = if_else(ch == '0001000' & substring(fish_id,1,3) == 'gen', 1, 0)) #only get flagged for removal if only seen in 2015 and identified by gen_id (spot-checking, looks like ~ 100 juvenile fish, which makes sense (actually filters out 107 fish))
#' # encounters_all_means <- encounters_all_means %>%
#' #   filter(rem_2015 == 0)
#' 
#' ######### Old code for linking up genetically-ided fish with tagged fish
#' # # Link up genetic ids with tag ones, if possible: go through allfish_mark, for any fish that has a gen_id-based fish_id, see if a fish with the same gen_id later has a tag_id and update the fish_id if so
#' # for(i in 1:length(allfish_mark$fish_id)) {
#' #   if(substring(allfish_mark$fish_id[i],1,3) == "gen") { #if it has a fish_id based on gen_id rather than tag_id...
#' #     gen_id_val <- allfish_mark$gen_id[i] #pull out the gen_id
#' #     
#' #     matches <- allfish_mark %>%
#' #       filter(gen_id == gen_id_val) #filter out other cases
#' #     
#' #     for(j in 1:length(matches$fish_id)){ #go through the fish entries that match the gen_id
#' #       if(!is.na(matches$tag_id[j])){ #if one of them has a tag_id
#' #         allfish_mark$fish_id[i] <- matches$fish_id[j] #update the fish_id to the fish_id based on the match with a tag_id
#' #       }
#' #     }
#' #   }
#' # }
#' 
#' ##### Incorporate genetic recaptures with tag recaptures - make identifier that uses cap_id if present, otherwise tag_id, or sample_id if neither 
#' # Possible cases: (new gen_id variable to indicate whether fish has a genotype, replaces cap_id, don't need to use sample_id b/c fish might have been clipped but not successfully sequenced)
#' # 1) clipped but not tagged, never caught again (sample_id Y, cap_id NA, tag_id NA), sampleXXXX - now gen_id Y, tag_id NA: genXXXXX  
#' # 2) clipped but not tagged, caught + clipped again (sample_id Y, cap_id Y, tag_id NA), capXXX - now gen_id Y, tag_id NA: genXXXX 
#' # 3) clipped but not tagged, tagged + clipped (sample_id Y, cap_id Y, tag_id NA in one Y in other), capXXXX - now gen_id Y, tag_id Y and NA: genXXX
#' # 4) clipped + tagged, never caught again (sample_id Y, tag_id Y, cap_id NA), tagXXXX - now gen_id Y, tag_id Y: genXXXX (should be tag? check if gen id carries through...)
#' # 5) tagged but not clipped, never caught again (sample_id NA, tag_id Y, cap_id NA), tagXXXX, now gen_id NA, tag_id Y: tagXXX
#' # 6) clipped + tagged 2016 or later, caught again via tag (sample_id NA (but soon!), tag_id Y, cap_id NA), tagXXXXX: now, gen_id NA, tag_id Y: tagXXXX
#' 
#' # # Fish with genetic and PIT markers - 3655 fish observations
#' # allfish_mark <- allfish %>%
#' #   filter(!is.na(gen_id) | !is.na(tag_id)) %>% #pull out fish "tagged" in any way, either PIT or via genetic sample
#' #   mutate(fish_id = case_when(!is.na(tag_id) ~ paste('tag', tag_id, sep=''), #if has a tag_id, use it in the fish_id
#' #                              is.na(tag_id) & !is.na(gen_id) ~ paste('gen', gen_id, sep=''))) # %>% #if it doesn't have a tag_id, use the gen_id in the fish_id - but some of these fish might later have a tag_id
#' # 
#' # # Fish with just PIT markers (to compare estimates) - 2376 fish observations
#' # allfish_tag <- allfish %>%
#' #   filter(!is.na(tag_id)) %>% #just pull out fish tagged with PIT-tags
#' #   dplyr::rename(fish_id = tag_id) #rename tag_id fish_id so works with CreateEncounterSummary function
#' # 
#' # # old version: this doesn't work b/c gen_id doesn't carry through all observations of the fish (like cap_id used to) - allfish_mark with the version is 3655 fish
#' # allfish_mark <- allfish %>% 
#' #   filter(!is.na(gen_id) | !is.na(tag_id)) %>% #pull out fish "tagged" in any way, either PIT or via genetic sample
#' #   mutate(fish_id = case_when(!is.na(gen_id) ~ paste('gen', gen_id, sep=''), #cases 1,2,3,4 above
#' #                              is.na(gen_id) & !is.na(tag_id) ~ paste('tag', tag_id, sep=''))) #cases 5,6 above
#' 
#' # # Link up genetic ids with tag ones, if possible: go through allfish_mark, for any fish that has a gen_id-based fish_id, see if a fish with the same gen_id later has a tag_id and update the fish_id if so
#' # for(i in 1:length(allfish_mark$fish_id)) {
#' #   if(substring(allfish_mark$fish_id[i],1,3) == "gen") { #if it has a fish_id based on gen_id rather than tag_id...
#' #     gen_id_val <- allfish_mark$gen_id[i] #pull out the gen_id
#' #     
#' #     matches <- allfish_mark %>%
#' #       filter(gen_id == gen_id_val) #filter out other cases
#' #     
#' #     for(j in 1:length(matches$fish_id)){ #go through the fish entries that match the gen_id
#' #       if(!is.na(matches$tag_id[j])){ #if one of them has a tag_id
#' #         allfish_mark$fish_id[i] <- matches$fish_id[j] #update the fish_id to the fish_id based on the match with a tag_id
#' #       }
#' #     }
#' #   }
#' # }
#' #   
#' # # Bit of investigating, testing to make sure the fish_id and for loop are working as expected...
#' # # look at a few cases to see if this is working as expected (not comprehensive, just looking at the dataframe)
#' # allfish_mark %>% filter(tag_id == 986112100165961) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -anem_table_id)
#' # allfish_mark %>% filter(tag_id == 985153000371766) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -anem_table_id) #hmm, this fish and the one below are the same, not linking up with my current fish_id
#' # allfish_mark %>% filter(gen_id == 1417) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -anem_table_id)
#' # allfish_mark %>% filter(tag_id == 985153000370613) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -anem_table_id)
#' # allfish_mark %>% filter(gen_id == 1392) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -anem_table_id)
#' # allfish_mark %>% filter(tag_id == 985153000355603) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -anem_table_id)
#' # allfish_mark %>% filter(gen_id == 1583) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -anem_table_id)
#' # 
#' # test_matches <- allfish_mark %>% filter(is.na(tag_id)) %>% #gen_ids 678, 772
#' #   mutate(id_type = substring(fish_id,1,3))
#' # table(test_matches$id_type) #86 that have a gen_id and get a tag_id-based fish_id added in - seems like the for loop might be working?
#' # 
#' # allfish_mark %>% filter(gen_id == 678) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -anem_table_id)
#' # allfish_mark %>% filter(gen_id == 772) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -anem_table_id)
#' 
#' ######## ORIGINAL WORK WITH SIZE AS A COVARIATE
#' # ##### Try running a MARK model with size-at-tagging as an individual constraint, site as grouping
#' # site_size_at_tagging <- function() {
#' #   # select relevant data
#' #   e2 <- encounters_all %>%
#' #     select(ch, tag_size, site) %>%
#' #     mutate(site = as.factor(site))
#' #   
#' #   e2_2 <- e2[complete.cases(e2),] #for using tag size, need no NAs...
#' # 
#' #   
#' #   # process data, make ddl
#' #   e2.processed = process.data(e2, begin.time=2015, groups="site")
#' #   e2.ddl = make.design.data(e2.processed)
#' #   
#' #   # process data, make ddl for complete cases
#' #   e2_2.processed = process.data(e2_2, begin.time=2015, groups="site")
#' #   e2_2.ddl = make.design.data(e2_2.processed)
#' #   
#' #   
#' #   # define models for Phi and p
#' #   Phi.dot = list(formula=~1, link="logit") #one survival for all sizes and sites and times
#' #   Phi.site = list(formula=~site, link="logit") #survival depends on site
#' #   Phi.size = list(formula=~tag_size, link="logit") #survival depends on size
#' #   Phi.size.x.site = list(formula=~tag_size*site, link="logit") #survival depends on site and site-dependent slope for size
#' #   Phi.time = list(formula=~time, link="logit") #survival varies with time, same for all sites
#' #   Phi.time.plus.site = list(formula=~time+site, link="logit")
#' #   
#' #   p.dot = list(formula=~1, link="logit") #one recap prob for all sizes, sites, and time
#' #   p.site = list(formula=~site, link="logit")
#' #   p.time = list(formula=~time, link="logit")
#' #   
#' #   # create model list
#' #   #e2.cml = create.model.list("CJS")
#' #   
#' #   #run and return models
#' #   #return(mark.wrapper(e2.cml, data=e2.processed, ddl=e2.ddl))
#' #   
#' #   #create models
#' #   e2.Phi.dot.p.dot = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot))
#' #   e2.Phi.site.p.dot = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.site, p=p.dot))
#' #   e2.Phi.size.p.dot = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.size, p=p.dot))
#' #   e2.Phi.size.p.site = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.size, p=p.site))
#' #   e2.Phi.size.x.site.p.site = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.size.x.site, p=p.site))
#' #   e2.Phi.site.p.site = mark(e2.processed, e2.ddl, model.parameters=list(Phi=Phi.site, p=p.site))
#' #   
#' #   e2_2.Phi.size.p.dot = mark(e2_2.processed, e2_2.ddl, model.parameters=list(Phi=Phi.size, p=p.dot))
#' #   e2_2.Phi.size.p.site = mark(e2_2.processed, e2_2.ddl, model.parameters=list(Phi=Phi.size, p=p.site))
#' #   e2_2.Phi.size.x.site.p.site = mark(e2_2.processed, e2_2.ddl, model.parameters=list(Phi=Phi.size.x.site, p=p.site))
#' #   
#' # }
#' # 
#' # e2.results = site_size_at_tagging()
#' # 
#' # 
#' # # Make some plots
#' # # e2.Phi.dot.p.dot (constant survival and recapture prob the whole time)
#' # e2.constant = as.data.frame(e2.Phi.dot.p.dot$results$beta) %>% 
#' #   mutate(param = c("Phi","p")) %>% #add a parameters column to make it easier to plot in ggplot
#' #   mutate(upper = logit_recip(ucl), lower = logit_recip(lcl), est = logit_recip(estimate)) #do the reciprocal transform on upper and lower confidence limits (need to check what those are - 95? SE? and that transforming them just straight up is the right way to go)
#' #   
#' # e2.site = as.data.frame(e2.Phi.site.p.site$results$beta) %>%
#' #   mutate(param = c(rep("Phi",15),rep("p",15))) %>% #add a parameter column to make it easier to plot in ggplot
#' #   mutate(site = c(rep(c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
#' #                         "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag"),2))) #add site column (Cabatoan is intercept)
#' # 
#' # e2.site.real = as.data.frame(e2.Phi.site.p.site$results$real) %>%
#' #   mutate(param = c(rep("Phi",15),rep("p",15))) %>% #add a parameter column to make it easier to plot in ggplot
#' #   mutate(site = c(rep(c("Cabatoan","Cardid Cemetery","Elementary School","Gabas","Haina","Hicgop South","Magbangon","Palanas",
#' #                         "Poroc Rose","Poroc San Flower","San Agustin","Sitio Baybayon","Tamakin Dacot","Visca","Wangag"),2))) %>% #add site column (Cabatoan is intercept)
#' #   mutate(row = seq(1,30,1)) #add row numbers to make plotting easier
#' # 
#' # #Think this is handled for now, should check again to be sure...
#' # #NEED TO ADD THE EFFECT OF EACH SITE TO THE INTERCEPT
#' # #### NEED TO SWITCH IT BACK FROM THE LOGIT SCALE FIRST!! AND THAT AFFECTS THE CONFIDENCE INTERVALS, RIGHT? CHECK THAT AGAIN....
#' # #### ALSO SHOULD PLOT +- SE TOO, SEE HOW THAT COMPARES/MATCHES UP TO LCL and UCL
#' # 
#' # # Constant Phi and p across time and site, no covariates (e2.Phi.dot.p.dot)
#' # pdf(file = here("Plots/PhiandpEstimates", "e2_Phiandp_constant.pdf"))
#' # ggplot(data = e2.constant, aes(param, est, color=param)) +
#' #   geom_pointrange(aes(ymin=lower, ymax=upper), size=1) +
#' #   xlab("parameter") + ylab("estimate") + ggtitle("e2.Phi.dot.p.dot (constant p and Phi)") +
#' #   scale_y_continuous(limits = c(0, 1)) +
#' #   theme_bw()
#' # dev.off()
#' #   
#' # # Constant Phi and p across time, variable across site, no covariates (e2.Phi.site.p.site)
#' # pdf(file = here("Plots/PhiandpEstimates", "e2_Phiandp_site.pdf"))
#' # ggplot(data = e2.site.real, aes(row, estimate, color=site, shape=param)) +
#' #   geom_pointrange(aes(ymin=lcl, ymax=ucl), size=1) +
#' #   xlab("parameter") + ylab("estimate") + ggtitle("e2.Phi.site.p.site (p and Phi by site)") +
#' #   #scale_y_continuous(limits = c(0, 1)) +
#' #   theme_bw()
#' # dev.off()
#' # 
#' # # Phi has tagging size as a covariate, not dependent on site or time, p constant (e2_2.Phi.size.p.dot)
#' # minsize = min(e2_2$tag_size) #right now, this is 1.6 (prob b/c of a typo.... Michelle is looking at it)
#' # maxsize = max(e2_2$tag_size)
#' # size.values = minsize+(0:30)*(maxsize-minsize)/30
#' # Phibysize = covariate.predictions(e2_2.Phi.size.p.dot,data=data.frame(tag_size=size.values),indices=c(1))
#' # 
#' # pdf(file = here("Plots/PhiandpEstimates", "e2_2Phisize.pdf"))
#' # ggplot(data = Phibysize$estimates, aes(covdata, estimate)) +
#' #   geom_ribbon(aes(ymin=lcl,ymax=ucl),color="light blue",fill="light blue") +
#' #   geom_line(color="black") +
#' #   xlab("size at tagging (cm)") + ylab("Phi estimate") + ggtitle("e2_2.Phi.size.p.dot") +
#' #   scale_y_continuous(limits = c(0, 1)) +
#' #   theme_bw()
#' # dev.off()
#' 
#' ###### ORIGINAL WORK WITH DISTANCE AS A TIME-VARYING INDIVIDUAL CONSTRAINT
#' # ##### Try running time-varying individual constraints - distance
#' # # data
#' # edist <- fish.Tagged %>%
#' #   select(ch, site, dist_2016, dist_2017, dist_2018) %>% #select just site and distance (not including 2015, first potential tagging year)
#' #   rename(dist2016 = dist_2016, dist2017 = dist_2017, dist2018 = dist_2018) %>% 
#' #   mutate(dist2016 = replace(dist2016, is.na(dist2016), 0)) %>% #replace NAs (before a fish was tagged) with 0s
#' #   mutate(dist2017 = replace(dist2017, is.na(dist2017), 0)) %>%
#' #   mutate(dist2018 = replace(dist2018, is.na(dist2018), 0))
#' # 
#' # # process data and make ddl
#' # edist.processed = process.data(edist, model="CJS", begin.time=2015)
#' # edist.ddl = make.design.data(edist.processed)
#' # 
#' # # set models for Phi
#' # Phi.dot = list(formula=~1, link="logit")
#' # Phi.time = list(formula=~time, link="logit")
#' # 
#' # # set models for p
#' # p.dot = list(formula=~1, link="logit")
#' # p.time = list(formula=~time, link="logit")
#' # p.dist = list(formula=~dist, link="logit")
#' # p.time.plus.dist = list(formula=~time+dist, link="logit")
#' # 
#' # # run some models
#' # edist.Phi.dot.p.dot = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot))
#' # edist.Phi.dot.p.time = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.dot, p=p.time)) #has issues running, can't estimate p in 2018
#' # edist.Phi.dot.p.dist = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist)) #think this ran!!!
#' # # edist.Phi.dot.p.dist.plus.time = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.dot, p=p.time.plus.dist))
#' # # edist.Phi.time.p.dot = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.dot))
#' # # edist.Phi.time.p.time = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.time))
#' # # edist.Phi.time.p.dist = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.dist))
#' # # edist.Phi.time.p.time.plus.dist = mark(edist.processed, edist.ddl, model.parameters=list(Phi=Phi.time, p=p.time.plus.dist))
#' # 
#' # # make some plots
#' # 
#' # # Phi constant, not dependent on site or time, p depends on distance to anem (edist.Phi.dot.p.dist)
#' # mindist = min(edist$dist2016,edist$dist2017,edist$dist2018) 
#' # maxdist = max(edist$dist2016,edist$dist2017,edist$dist2018) #seems kind of high and probably due to an error
#' # maxdist2 = 400 #encompasses all but the highest value (which is an order of magnitude higher - 4874 compared to next highest of 352, both in 2017)
#' # maxdist3 = 200 #encompasses highest values in 2016, 2018 
#' # #tail(sort(edist$dist2016),10)
#' # #tail(sort(edist$dist2017),10) #really 2017 that has the high distances...
#' # #tail(sort(edist$dist2018),10)
#' # dist.values = mindist+(0:30)*(maxdist-mindist)/30
#' # dist.values2 = mindist+(0:30)*(maxdist2-mindist)/30
#' # dist.values3 = mindist+(0:30)*(maxdist3-mindist)/30
#' # 
#' # # Not getting this to work - it doesn't change with distance... maybe/probably I'm specifying the index incorrectly?
#' # #pbydist1 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values),indices=c(7))
#' # #pbydist2 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values2),indices=c(1))
#' # #pbydist3 = covariate.predictions(edist.Phi.dot.p.dist,data=data.frame(dist=dist.values3),indices=c(1))
#' # 
#' # edist.Phi.dot.p.dist.results = as.data.frame(edist.Phi.dot.p.dist$results$beta) 
#' # 
#' # # Large range of distances (up to 5000m)
#' # pbydist_range1 = data.frame(dist = dist.values) %>%
#' #   mutate(p_logit = (edist.Phi.dot.p.dist.results$estimate[2]) + (edist.Phi.dot.p.dist.results$estimate[3]*dist)) %>%
#' #   mutate(p_lcl_logit = edist.Phi.dot.p.dist.results$lcl[2] + (edist.Phi.dot.p.dist.results$lcl[3]*dist)) %>%
#' #   mutate(p_ucl_logit = edist.Phi.dot.p.dist.results$ucl[2] + (edist.Phi.dot.p.dist.results$ucl[3]*dist)) %>%
#' #   mutate(p = logit_recip(p_logit)) %>%
#' #   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#' #   mutate(p_ucl = logit_recip(p_ucl_logit)) 
#' # 
#' # # Mid range of distances (up to 400m)
#' # pbydist_range2 = data.frame(dist = dist.values2) %>%
#' #   mutate(p_logit = (edist.Phi.dot.p.dist.results$estimate[2]) + (edist.Phi.dot.p.dist.results$estimate[3]*dist)) %>%
#' #   mutate(p_lcl_logit = edist.Phi.dot.p.dist.results$lcl[2] + (edist.Phi.dot.p.dist.results$lcl[3]*dist)) %>%
#' #   mutate(p_ucl_logit = edist.Phi.dot.p.dist.results$ucl[2] + (edist.Phi.dot.p.dist.results$ucl[3]*dist)) %>%
#' #   mutate(p = logit_recip(p_logit)) %>%
#' #   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#' #   mutate(p_ucl = logit_recip(p_ucl_logit)) 
#' #   
#' # # Short range of distances (up to 200m)
#' # pbydist_range3 = data.frame(dist = dist.values3) %>%
#' #   mutate(p_logit = (edist.Phi.dot.p.dist.results$estimate[2]) + (edist.Phi.dot.p.dist.results$estimate[3]*dist)) %>%
#' #   mutate(p_lcl_logit = edist.Phi.dot.p.dist.results$lcl[2] + (edist.Phi.dot.p.dist.results$lcl[3]*dist)) %>%
#' #   mutate(p_ucl_logit = edist.Phi.dot.p.dist.results$ucl[2] + (edist.Phi.dot.p.dist.results$ucl[3]*dist)) %>%
#' #   mutate(p = logit_recip(p_logit)) %>%
#' #   mutate(p_lcl = logit_recip(p_lcl_logit)) %>%
#' #   mutate(p_ucl = logit_recip(p_ucl_logit)) 
#' # 
#' # # Long dist values (up to 5000m from anem)
#' # pdf(file = here("Plots/PhiandpEstimates", "edist_Phidot_pdist_distvalues1.pdf"))
#' # ggplot(data = pbydist_range1, aes(dist, p)) +
#' #   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#' #   geom_line(color="black") +
#' #   xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("edist.Phi.dot.p.dist (constant Phi, p by dist)") +
#' #   scale_y_continuous(limits = c(0, 1)) +
#' #   theme_bw()
#' # dev.off()
#' # 
#' # # Mid range dist values (up to 400m from anem)
#' # pdf(file = here("Plots/PhiandpEstimates", "edist_Phidot_pdist_distvalues2.pdf"))
#' # ggplot(data = pbydist_range2, aes(dist, p)) +
#' #   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#' #   geom_line(color="black") +
#' #   xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("edist.Phi.dot.p.dist (constant Phi, p by dist)") +
#' #   scale_y_continuous(limits = c(0, 1)) +
#' #   theme_bw()
#' # dev.off()
#' # 
#' # # Short range dist values (up to 200m from anem)
#' # pdf(file = here("Plots/PhiandpEstimates", "edist_Phidot_pdist_distvalues3.pdf"))
#' # ggplot(data = pbydist_range3, aes(dist, p)) +
#' #   geom_ribbon(aes(ymin=p_lcl,ymax=p_ucl),color="light blue",fill="light blue") +
#' #   geom_line(color="black") +
#' #   xlab("distance from anem (m)") + ylab("p estimate") + ggtitle("edist.Phi.dot.p.dist (constant Phi, p by dist)") +
#' #   scale_y_continuous(limits = c(0, 1)) +
#' #   theme_bw()
#' # dev.off()
#' 
#' 
#' ##### Try running MARK (again), this time with time-varying individual constraints!
#' # Trying one without site as a group, just with distance to anem as the time-varying individual constraint
#' # data
#' e1 <- encounters_means %>% 
#'   select(ch, dist2016, dist2017, dist2018)
#' 
#' # process data and make ddl
#' fish.processed = process.data(e1, model="CJS", begin.time=2015)
#' fish.ddl = make.design.data(fish.processed)
#' 
#' # set models for Phi 
#' Phi.dot = list(formula=~1, link="logit")
#' Phi.time = list(formula=~time, link="logit")
#' 
#' # set models for p
#' p.dot = list(formula=~1, link="logit")
#' p.time = list(formula=~time, link="logit")
#' p.dist = list(formula=~dist, link="logit")
#' p.dist.plus.time = list(formula=~dist+time, link="logit")
#' 
#' # run models individually
#' e1.Phi.dot.p.dot = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot))
#' e1.Phi.dot.p.time = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.dot, p=p.time))
#' e1.Phi.dot.p.dist = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist))
#' e1.Phi.dot.p.dist.plus.time = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.dot, p=p.dist.plus.time))
#' e1.Phi.time.p.dot = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.time, p=p.dot))
#' e1.Phi.time.p.time = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.time, p=p.time))
#' e1.Phi.time.p.dist = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.time, p=p.dist))
#' e1.Phi.time.p.dist.plus.time = mark(fish.processed, fish.ddl, model.parameters=list(Phi=Phi.time, p=p.dist.plus.time))
#' 
#' # run models using clm and wrapper
#' e1.cml = create.model.list("CJS")
#' e1.model.wrap = mark.wrapper(e1.cml, data=fish.processed, ddl=fish.ddl)
#' 
#' e1.model.wrap[[1]]$design.matrix # this is one of the ones with distance... not convinced it worked right..
#' e1.model.wrap[[2]]$design.matrix #
#' e1.model.wrap[[5]]$design.matrix
#' 
#' e1.collect = collect.models(table=TRUE)
#' # make plots to compare the models
#' 
#' # compare using AIC?
#' 
#' # run models using site as a group
#' 
#' # run models with tagging size as covariate
#' 
#' # run models with actual size as a covariate (can fill in missing data for un-recaught fish using Michelle's growth? or just assume it stays the same if not ever caught again, mean between two sizes on either size if one caputure missed, not sure what would do for fish pre-capture)
#' 
#' # run models with distance where NAs were replaced with 0 and where actual distance to anem was calculated whether or not fish had been tagged yet, compare all three methods to make sure it doesn't matter
#' 
#' # make some sort of data frame with estimate, 95% confidence bounds, with different assumptions (time, distance, site, etc. matters) shown, so can use as want in persistence metrics and analyses
#' 
#' ##### Now, try running MARK! (original attempt)
#' # specifying site as group
#' # make process and design data
#' #with site
#' tagged.site.process = process.data(encounters_all, model="CJS", begin.time=2015, groups="site")
#' tagged.site.ddl = make.design.data(tagged.site.process)
#' 
#' #without site
#' tagged.process = process.data(encounters_all, model="CJS", begin.time=2015)
#' tagged.ddl = make.design.data(tagged.process)
#' 
#' # set different relationships
#' Phi.dot = list(formula=~1)
#' Phi.time = list(formula=~time)
#' Phi.site = list(formula=~site)
#' Phi.siteplustime = list(formua=~site+time)
#' 
#' p.dot = list(formula=~1)
#' p.time = list(formula=~time)
#' p.site = list(formula=~site)
#' p.siteplustime = list(formua=~site+time) # does this still pool sites at some level, assuming they come from the same dist and are related (like would in rethinking) or just break them out completely separately?
#' 
#' # run some models
#' tagged.phi.dot.p.dot = mark(tagged.process, tagged.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot))
#' tagged.site.phi.dot.p.dot = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot))
#' tagged.site.phi.time.p.dot = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.time, p=p.dot))
#' tagged.site.phi.dot.p.time = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.dot, p=p.time))
#' tagged.site.phi.time.p.time = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.time, p=p.time))
#' tagged.site.phi.site.p.site = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.site, p=p.site))
#' tagged.site.phi.time.p.siteplustime = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.site, p=p.siteplustime))
#' 
#' #these don't run... must have something wrong? 
#' tagged.site.phi.siteplustime.p.time = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.siteplustime, p=p.time))
#' tagged.site.phi.siteplustime.p.siteplustime = mark(tagged.site.process, tagged.site.ddl, model.parameters=list(Phi=Phi.siteplustime, p=p.siteplustime))
#' 
#' #tagged.site <- mark(data=encounters_all, model = "CJS")
#' 
#' # Using the output
#' ## Note from MARK book - the Beta estimates at the top relate to the link function - thethey're the values actually estimated before the actual paramter values are re-constituted; they're estimates on the logit scale
#' ## So for survival and recap, the Beta estimates are the logit values (ln(p/(1-p))), for linear model constraints, they are the coefficients (slopes) for each term in the linear model
#' ## See pp.226-229 in MARK book for details
#' #################################### Testing out doing mark-recap in rethinking ############# 
#' 
#' 
#' 
#' 
#' #################################### Old stuff below ############# 
#' 
#' ##### Previous work using own likelihood functions...
#' # Now, create encounter histories for fish at all sites combined
#' encountersAll <- CreateEncounterSummary(2015, 2017, allfish) #simpler way of getting encounter history like above using function
#' nvals_allfish <- as.data.frame(table(encountersAll$encounter.hist))
#' 
#' # Try gradient descent on this
#' params0 <- c(0.5, 0.5, 0.5, 0.5) #s1, s2, p2, p3
#' params0_2 <- c(0.5, 0.5) #s, p
#' highVal <- 0.999
#' lowVal <- 0.001
#' eps <- 0.001
#' del <- 0.001
#' p_vec <- seq(0.01, 0.99, 0.01) #recapture prob
#' s_vec <- seq(0.01, 0.99, 0.01) #survival prob
#' Val1 <- c("001","010","011","100","101","110","111")
#' 
#' 
#' all_GradDesc <- GradientDescent3SY4Params(nvals_allfish, params0, 50, 0.001, 0.001, 0.999, 0.001) #works! except param 2 and param 4 are exactly the same...
#' #param1 (s1) = 0.15, param2 (s2) = 0.36, param3 (p2) = 0.35, param4 (p3) = 0.36, like = -345.4786
#' all_GradDesc_1000tsteps <- GradientDescent3SY4Params(nvals_allfish, params0, 1000, 0.001, 0.001, 0.999, 0.001) #param 2 and 4 stll exactly the same... seems weird
#' all_GradDesc_1000tsteps_2p <- GradientDescent3SY2Params(nvals_allfish, params0_2, 1000, 0.001, 0.001, 0.999, 0.001) #param 2 and 4 stll exactly the same... seems weird
#' 
#' # Now try for different sites
#' # Cabatoan
#' encountersCabatoan <- CreateEncounterSummary(2015, 2017, (allfish %>% filter(site == site_vec[Cabatoan]))) #simpler way of getting encounter history like above using function
#' nvals_Cabatoan <- as.data.frame(table(encountersCabatoan$encounter.hist))
#' nvals_Cabatoan_clean <- c(0, 15, 0, 31, 0, 0, 0)
#' nvals_Cabat <- as.data.frame(Val1)
#' nvals_Cabat$Freq <- nvals_Cabatoan_clean
#' 
#' Cabat_ParamRangeScan <- ParamRangeCheck3SY2P_dfoutput(s_vec, p_vec, nvals_Cabat)
#' smaxCabat <- Cabat_ParamRangeScan$survival[which.max(Cabat_ParamRangeScan$like)]
#' pmaxCabat <- Cabat_ParamRangeScan$capture[which.max(Cabat_ParamRangeScan$like)]
#' Cabat_GradDesc_1000tsteps <- GradientDescent3SY4Params(nvals_Cabat, params0, 1000, 0.001, 0.001, 0.999, 0.001) #param 2 and 4 stll exactly the same... seems weird
#' 
#' 
#' N010 <- 2
#' N011 <- 3
#' N100 <- 4
#' N101 <- 5
#' N110 <- 6
#' N111 <- 7
#' # Try gradient descent on this
#' params0 <- c(0.5, 0.5, 0.5, 0.5) #s1, s2, p2, p3
#' highVal <- 0.999
#' lowVal <- 0.001
#' eps <- 0.001
#' del <- 0.001
#' all_GradDesc <- GradientDescent3SY4Params(nvals_allfish, params0, 50, 0.001, 0.001, 0.999, 0.001) #works! except param 2 and param 4 are exactly the same...
#' #param1 (s1) = 0.15, param2 (s2) = 0.36, param3 (p2) = 0.35, param4 (p3) = 0.36, like = -345.4786
#' all_GradDesc_1000tsteps <- GradientDescent3SY4Params(nvals_allfish, params0, 1000, 0.001, 0.001, 0.999, 0.001) #param 2 and 4 stll exactly the same... seems weird
#' 
#' 
#' # Caridad Cemetery
#' 
#' # Caridad Proper
#' 
#' # Elementary School
#' 
#' # Gabas
#' 
#' # Haina
#' 
#' # Hicgop South
#' 
#' # Magbangon
#' 
#' # Palanas
#' 
#' # Poroc Rose
#' 
#' # Poroc San Flower
#' 
#' 
#' # Pop size
#' Pal <- allfish %>% filter(site == "Palanas")
#' Pal_pop <- as.data.frame(table(Pal$year))
#' 
#' Wan <- allfish %>% filter(site == "Wangag")
#' Wan_pop <- as.data.frame(table(Wan$year))
#' 
#' SB <- allfish %>% filter(site == "Sitio Baybayon")
#' SB_pop <- as.data.frame(table(SB$year))
#' 
#' Vis <- allfish %>% filter(site == "Visca")
#' Vis_pop <- as.data.frame(table(Vis$year))
#' 
#' 
#' 
#' 
#' #################### Old code: ####################
#' # # Remove the NA tag id at the end (for 1111 fish)
#' # encounters_all <- encounters_all[1:(length(encounters_all$tag_id)-1),] 
#' 
#' # Data sorting, from when there was cap_id and not gen_id
#' # # old code from when there was cap_id and before gen_id         
#' # allfish_mark <- allfish %>%
#' #   filter(!is.na(cap_id) | !is.na(sample_id) | !is.na(tag_id)) %>% #pull out fish "tagged" in any way, either PIT or via fin-clip
#' #   mutate(fish_id = case_when(!is.na(sample_id) & is.na(cap_id) & is.na(tag_id) ~ paste("sample", sample_id, sep = ""), #case 1 above
#' #                              !is.na(sample_id) & !is.na(cap_id) ~ paste("cap", cap_id, sep = ""), #cases 2, 3 above
#' #                              is.na(cap_id) & !is.na(tag_id) ~ paste("tag", tag_id, sep = ""))) #cases 4, 5, 6 above
#' #                              
#' # # do a few tests to see if this is working as expected... (not comprehensive, just looking at the dataframe)
#' # allfish_mark %>% filter(cap_id == 46) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -fish_table_id, -anem_table_id)
#' # allfish_mark %>% filter(sample_id == "APCL12_090") %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -fish_table_id, -anem_table_id)
#' # allfish_mark %>% filter(tag_id == 985153000371766) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -fish_table_id, -anem_table_id) #this fish was caught in 4 years!
#' # allfish_mark %>% filter(tag_id == 985153000370613) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -fish_table_id, -anem_table_id) #seems unlikely this was YR at 12.2cm... was YP at 12.3 the next year...
#' # allfish_mark %>% filter(cap_id == 34) %>% select(-anem_obs, -anem_id, -old_anem_id, -dive_type, -fish_table_id, -anem_table_id)
#' # 
#' # dim(allfish %>% filter(!is.na(sample_id) | !is.na(tag_id) | !is.na(gen_id))) #4788 fish
#' # dim(allfish %>% filter(!is.na(gen_id) | !is.na(tag_id))) #3655 fish - seems like a lot are lost!
#' # dim(allfish %>% filter(!is.na(tag_id))) #2376 fish
#' # dim(allfish %>% filter(!is.na(sample_id))) #4370 fish
#' # dim(allfish %>% filter(!is.na(sample_id) & is.na(gen_id) & is.na(tag_id))) #1133 - this is a lot of fish - does this include the sample_ids just put in the sequencer now? In that case, though, seems too small...
#' # table((allfish %>% filter(!is.na(sample_id)))$year) #2231 for 2016 (814), 2017 (647), and 2018 (770)
#' # table((allfish %>% filter(!is.na(gen_id)))$year)
#' 
#' # # compare how many more fish this includes
#' # # just tags - 2376 fish
#' # allfish_tags <- allfish %>%
#' #   filter(!is.na(tag_id)) %>%
#' #   group_by(tag_id) %>%
#' #   mutate(nobs = n())
#' # # tags + cap_ids  - 4788 fish - lots more!
#' # allfish_tagsfins <- allfish %>%
#' #   filter(!is.na(cap_id) | !is.na(sample_id) | !is.na(tag_id)) %>% #pull out fish "tagged" in any way, either PIT or via fin-clip
#' #   mutate(fish_id = case_when(!is.na(sample_id) & is.na(cap_id) & is.na(tag_id) ~ paste("sample", sample_id, sep = ""), #case 1 above
#' #                              !is.na(sample_id) & !is.na(cap_id) ~ paste("cap", cap_id, sep = ""), #cases 2, 3 above
#' #                              is.na(cap_id) & !is.na(tag_id) ~ paste("tag", tag_id, sep = ""))) %>% #cases 4, 5, 6 above
#' #   group_by(fish_id) %>%
#' #   mutate(nobs = n())
#' 
#' # Old, from when there were cap_ids
#' # # Pull out all the tags and their encounter history from 2012-2018, add column to see which fish are IDed by tag, cap, or sample
#' # encounters_all <- CreateEncounterSummary(2012, 2018, allfish_mark) %>% #simpler way of getting encounter history like above using function
#' #   mutate(id_type = case_when(substring(fish_id,1,1) == "t" ~ "tag", #1851 fish
#' #                              substring(fish_id,1,1) == "c" ~ "cap", #135 fish
#' #                              substring(fish_id,1,1) == "s" ~ "sample")) #2173 fish
#' 
#' # ######### Old size and color by year code
#' # # Re-adding the tail color by year and size by year data for KC to use at Michigan workshop
#' # encounters_all_fulldata <- encounters_all
#' # ## Code that adds in tail color by year, problem with that is that we only have it for fish that are caught, which isn't allowed for time-varying individual covariates (have to have a value for each fish at each time)
#' # ##might be able to get around that by just doing the color at the last time it was caught but not sure that's worth it so for now, commenting out...
#' # # Add in tail color (could probably do these all together in one command... look into that later..)
#' # tail_color_2015 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#' #   filter(year == 2015) %>%
#' #   group_by(tag_id) %>%
#' #   summarize(tail_color_2015 = color[1])
#' # 
#' # tail_color_2016 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#' #   filter(year == 2016) %>%
#' #   group_by(tag_id) %>%
#' #   summarize(tail_color_2016 = color[1])
#' # 
#' # tail_color_2017 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#' #   filter(year == 2017) %>%
#' #   group_by(tag_id) %>%
#' #   summarize(tail_color_2017 = color[1])
#' # 
#' # tail_color_2018 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#' #   filter(year == 2018) %>%
#' #   group_by(tag_id) %>%
#' #   summarize(tail_color_2018 = color[1])
#' # 
#' # #encounters[[i]] <- tagged.fish %>% group_by(tag_id) %>% summarise(!!var.name := ifelse(sum(year == sample.years[i])>0, 1, 0)) #create data frames for each year that have a vector of tag ids and a vector of encountered (1) or didn't (0) in
#' # 
#' # encounters_all_fulldata <- left_join(encounters_all_fulldata, tail_color_2015, by="tag_id")
#' # encounters_all_fulldata <- left_join(encounters_all_fulldata, tail_color_2016, by="tag_id")
#' # encounters_all_fulldata <- left_join(encounters_all_fulldata, tail_color_2017, by="tag_id")
#' # encounters_all_fulldata <- left_join(encounters_all_fulldata, tail_color_2018, by="tag_id")
#' # 
#' # # 
#' # # ## Same issue with size in each year as with color above - don't have a value for the fish that weren't recaptured
#' # # #Might be able to add some in based on growth curves using the relationships Michelle has, could get into that at some point but for now ignoring
#' # # # Add in size in each of the years
#' # size_2015 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#' #   filter(year == 2015) %>%
#' #   group_by(tag_id) %>%
#' #   summarize(size_2015 = mean(size, rm.na = TRUE))
#' # 
#' # size_2016 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#' #   filter(year == 2016) %>%
#' #   group_by(tag_id) %>%
#' #   summarize(size_2016 = mean(size, rm.na = TRUE))
#' # 
#' # size_2017 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#' #   filter(year == 2017) %>%
#' #   group_by(tag_id) %>%
#' #   summarize(size_2017 = mean(size, rm.na = TRUE))
#' # 
#' # size_2018 <- allfish %>% filter(tag_id %in% encounters_all$tag_id) %>%
#' #   filter(year == 2018) %>%
#' #   group_by(tag_id) %>%
#' #   summarize(size_2018 = mean(size, rm.na = TRUE))
#' # 
#' # encounters_all_fulldata <- left_join(encounters_all_fulldata, size_2015, by="tag_id")
#' # encounters_all_fulldata <- left_join(encounters_all_fulldata, size_2016, by="tag_id")
#' # encounters_all_fulldata <- left_join(encounters_all_fulldata, size_2017, by="tag_id")
#' # encounters_all_fulldata <- left_join(encounters_all_fulldata, size_2018, by="tag_id")
#' 
#' 
#' # # Function to calculate likelihood for 3 years of sampling (check this!!) with fish tagged in year 2 included too
#' # Like3SampleYears4vars <- function(n111, n110, n101, n100, n011, n010, s1, s2, p2, p3) {
#' #   #like <- ((s1*p2*s2*p3)^n111)*((s1*p2*(1-s2*p3))^n110)*((s1*(1-p2)*s2*p3)^n101)*((1-s1*p2-s1*(1-p2)*s2*p3)^n100)
#' #   x2 <- 1 - s2 + s2*(1 - p3) #probability of not being seen after the second sampling session
#' #   x1 <- 1 - s1 + s1*(1-p2)*x2 #probability of not being seen after the first sampling session
#' #   
#' #   p111 <- s1*p2*s2*p3 #probability of encounter history 111
#' #   p110 <- s1*p2*x2 #probability of encounter history 110
#' #   p101 <- s1*(1 - p2)*s2*p3 #probability of encounter history 101
#' #   p100 <- x1 #probability of encounter history 100
#' #   p011 <- s2*p3 #probability of encounter history 011
#' #   p010 <- x2 #probability of encounter history 010
#' #   
#' #   like <- n111*log(p111) + n110*log(p110) + n101*log(p101) + n100*log(p100) + n011*log(p011) + n010*log(p010)
#' #   return(like)
#' # }
#' # 
#' # # Same as above but for only two variables (assuming that survival and recapture probs are same through time)
#' # Like3SampleYears2vars <- function(n111, n110, n101, n100, n011, n010, s, p) {
#' #   #like <- ((s1*p2*s2*p3)^n111)*((s1*p2*(1-s2*p3))^n110)*((s1*(1-p2)*s2*p3)^n101)*((1-s1*p2-s1*(1-p2)*s2*p3)^n100)
#' #   x2 <- 1 - s + s*(1 - p) #probability of not being seen after the second sampling session
#' #   x1 <- 1 - s + s*(1-p)*x2 #probability of not being seen after the first sampling session
#' #   
#' #   p111 <- s*p*s*p #probability of encounter history 111
#' #   p110 <- s*p*x2 #probability of encounter history 110
#' #   p101 <- s*(1 - p)*s*p #probability of encounter history 101
#' #   p100 <- x1 #probability of encounter history 100
#' #   p011 <- s*p #probability of encounter history 011
#' #   p010 <- x2 #probability of encounter history 010
#' #   
#' #   like <- n111*log(p111) + n110*log(p110) + n101*log(p101) + n100*log(p100) + n011*log(p011) + n010*log(p010)
#' #   return(like)
#' # }
#' # 
#' # # Just cycles through vectors of possible survival and recapture probabilities and calculates the likelihood
#' # ParamRangeCheck3SY2P_dfoutput <- function(s_vec, p_vec, nvals) { 
#' #   #set index values for various encouter histories (this is the way R orders them using table)
#' #   N010 <- 2
#' #   N011 <- 3
#' #   N100 <- 4
#' #   N101 <- 5
#' #   N110 <- 6
#' #   N111 <- 7
#' #   
#' #   #get numbers of each type of encounter history from the inputs
#' #   n111 <- nvals$Freq[N111]
#' #   n110 <- nvals$Freq[N110]
#' #   n101 <- nvals$Freq[N101]
#' #   n100 <- nvals$Freq[N100]
#' #   n011 <- nvals$Freq[N011]
#' #   n010 <- nvals$Freq[N010]
#' #   
#' #   survival <- data.frame(survival=rep(s_vec,length(p_vec)))
#' #   out <- survival
#' #   capture_list <- rep(p_vec[1], length(s_vec))
#' #   for (i in 2:length(p_vec)) {
#' #     cap_toadd <- rep(p_vec[i], length(s_vec))
#' #     capture_list <- c(capture_list, cap_toadd)
#' #   }
#' #   out$capture <- capture_list
#' #   like = NULL
#' #   
#' #   for (i in 1:dim(out)[1]) {
#' #     like1 <- Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, out$survival[i], out$capture[i])
#' #     like <- c(like, like1)
#' #   }
#' #   out$like=like
#' #   return(out)
#' # }
#' # 
#' # # Gradient descent for 3 sample years, estimating 4 parameters (doesn't really seem to work, oscillates between all params essentially at 0 or all at 1)
#' # GradientDescent3SY4Params <- function(nvals, params0, nsteps, eps, del, highVal, lowVal) {
#' #   param1 <- rep(NA,1,nsteps)
#' #   out <- as.data.frame(param1)
#' #   out$param2 <- rep(NA,1,nsteps)
#' #   out$param3 <- rep(NA,1,nsteps)
#' #   out$param4 <- rep(NA,1,nsteps)
#' #   out$like <- rep(NA,1,nsteps)
#' #   
#' #   params <- params0
#' #   
#' #   #set index values for various encouter histories (this is the way R orders them using table) (001 is 1)
#' #   N010 <- 2
#' #   N011 <- 3
#' #   N100 <- 4
#' #   N101 <- 5
#' #   N110 <- 6
#' #   N111 <- 7
#' #   
#' #   #get numbers of each type of encounter history from the inputs
#' #   n111 <- nvals$Freq[N111]
#' #   n110 <- nvals$Freq[N110]
#' #   n101 <- nvals$Freq[N101]
#' #   n100 <- nvals$Freq[N100]
#' #   n011 <- nvals$Freq[N011]
#' #   n010 <- nvals$Freq[N010]
#' #   
#' #   for (i in 1:nsteps) { #doesn't work with Palanas data b /c get NANs
#' #     #store parameter values and likelihood
#' #     out$param1[i] <- params[1]
#' #     out$param2[i] <- params[2]
#' #     out$param3[i] <- params[3]
#' #     out$param4[i] <- params[4]
#' #     out$like[i] <- Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4])
#' #     
#' #     #calculate the gradient and update the parameters
#' #     df <- rep(NA, length(params))
#' #     df[1] <- (Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1]+eps, params[2], params[3], params[4]) - Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4]))/eps
#' #     df[2] <- (Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2]+eps, params[3], params[4]) - Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4]))/eps
#' #     df[3] <- (Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3]+eps, params[4]) - Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4]))/eps
#' #     df[4] <- (Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4]+eps) - Like3SampleYears4vars(n111, n110, n101, n100, n011, n010, params[1], params[2], params[3], params[4]))/eps
#' #     
#' #     params <- params+del*df
#' #     
#' #     #if parameters have become <= 0 or >= 1, set to low val or high val
#' #     for (i in 1:length(params)) {
#' #       if (params[i] >= 1) {
#' #         params[i] <- highVal
#' #       }
#' #       if (params[i] <= 0) {
#' #         params[i] <- lowVal
#' #       }
#' #     }
#' #     #print(Like3SampleYears(n111, n110, n101, n100, params[1], params[2], params[3], params[4]))
#' #   }
#' #   return(out)
#' # }
#' # 
#' # # Gradient Descent 3 years, 2 params
#' # GradientDescent3SY2Params <- function(nvals, params0, nsteps, eps, del, highVal, lowVal) {
#' #   survival_p <- rep(NA,1,nsteps)
#' #   out <- as.data.frame(survival_p)
#' #   out$recapture_p <- rep(NA,1,nsteps)
#' #   #out$param3 <- rep(NA,1,nsteps)
#' #   #out$param4 <- rep(NA,1,nsteps)
#' #   out$like <- rep(NA,1,nsteps)
#' #   
#' #   params <- params0
#' #   
#' #   #set index values for various encouter histories (this is the way R orders them using table) (001 is 1)
#' #   N010 <- 2
#' #   N011 <- 3
#' #   N100 <- 4
#' #   N101 <- 5
#' #   N110 <- 6
#' #   N111 <- 7
#' #   
#' #   #get numbers of each type of encounter history from the inputs
#' #   n111 <- nvals$Freq[N111]
#' #   n110 <- nvals$Freq[N110]
#' #   n101 <- nvals$Freq[N101]
#' #   n100 <- nvals$Freq[N100]
#' #   n011 <- nvals$Freq[N011]
#' #   n010 <- nvals$Freq[N010]
#' #   
#' #   for (i in 1:nsteps) { #doesn't work with Palanas data b /c get NANs
#' #     #store parameter values and likelihood
#' #     out$survival_p[i] <- params[1]
#' #     out$recapture_p[i] <- params[2]
#' #     out$like[i] <- Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, params[1], params[2])
#' #     
#' #     #calculate the gradient and update the parameters 
#' #     df <- rep(NA, length(params))
#' #     df[1] <- (Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, params[1]+eps, params[2]) - Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, params[1], params[2]))/eps
#' #     df[2] <- (Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, params[1], params[2]+eps) - Like3SampleYears2vars(n111, n110, n101, n100, n011, n010, params[1], params[2]))/eps
#' #     
#' #     params <- params+del*df
#' #     
#' #     #if parameters have become <= 0 or >= 1, set to low val or high val
#' #     for (i in 1:length(params)) {
#' #       if (params[i] >= 1) {
#' #         params[i] <- highVal
#' #       }
#' #       if (params[i] <= 0) {
#' #         params[i] <- lowVal
#' #       }
#' #     }
#' #   }
#' #   return(out)
#' # }
#' # 
#' # 
#' # # This version over-biases recaps b/c doesn't include fish that just have a sample_id - were clipped (so genetic "tag" taken) but not recaptured
#' # allfish_mark <- allfish %>%
#' #   filter(!is.na(cap_id) | !is.na(tag_id)) %>% #filter for fish with a cap_id, sample_id, or tag_id
#' #   mutate(fish_id = ifelse(is.na(cap_id), paste("tag", tag_id, sep=""), paste("cap", cap_id, sep = "")))
#' 
#' # #################################### Old R MARK work ############# 
#' # # need to get data in format for R Mark
#' # # column ch - has cenounter history (like 1011001)
#' # # then columns with characteristics - should have tag_id, size (at first capture), site, could have distance to any anem it was caught it to that year's track for each year...
#' # #encounters[[i]] <- tagged.fish %>% group_by(tag_id) %>% summarise(!!var.name := ifelse(sum(year == sample.years[i])>0, 1, 0)) #create data frames for each year that have a vector of tag ids and a vector of encountered (1) or didn't (0) in 
#' # 
#' # 
#' # #Load practice data sets
#' # data(dipper)
#' # 
#' # #Example from the RMark pdf, just to test mark
#' # dipper.Phidot.pdot=mark(dipper,threads=1)
#' # 
#' # data(example.data)
#' 
#' # # Go through the list of tagged fish, if it could have been caught in a year, assign it the anem where it was caught or most recently caught
#' # findAnemByYear <- function(df, genortag, fish.AllInfo) { #df is data frame with fish encounter history (encounters_all or _tag), genortag should be either "gentag" or "tag", fish.AllInfo is df
#' #   # go down different path, depending whether data is from gen+tag combined or tag-only
#' #   if(genortag == "gentag"){
#' #     for(i in 1:length(df$fish_id)) {
#' #       id <- df$fish_id[i]
#' #       year_tag <- df$first_capture_year[i]
#' #   
#' #       #Check for multiple anems in capture year (should I take recapture dives out?)
#' #       if(length((allfish_mark %>% filter(fish_id == id, year == df$first_capture_year[i]))$anem_id) > 1) { print(paste("Multiple matching anems for fish",tag,"in year",year_tag))}
#' #       df$capture_anem[i] = (fish.AllInfo %>% filter(fish_id == id, year == year_tag))$anem_id[1] #record the anem where fish was first captured
#' #       
#' #       # Go through the years - hmm, this is taking the actual anem it was caught at in each of the years (which is known for fish caught and not known for fish not caught). What about just using the original anem caught at and calculating distance from that in each of the years?
#' #       if(year_tag == 2013){
#' #         df$anem_2013[i] = (fish.AllInfo %>% filter(fish_id == id, year == year_tag))$anem_id[1]
#' #         df$anem_2014[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2014))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2014))$anem_id, df$anem_2013[i])
#' #         df$anem_2015[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2015))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2015))$anem_id, df$anem_2014[i])
#' #         df$anem_2016[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2016))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2016))$anem_id, df$anem_2015[i])
#' #         df$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2017))$anem_id, df$anem_2016[i])
#' #         df$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2018))$anem_id, df$anem_2017[i])
#' #       } else if(year_tag == 2014) {
#' #         df$anem_2014[i] = (fish.AllInfo %>% filter(fish_id == id, year == year_tag))$anem_id[1]
#' #         df$anem_2015[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2015))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2015))$anem_id, df$anem_2014[i])
#' #         df$anem_2016[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2016))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2016))$anem_id, df$anem_2015[i])
#' #         df$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2017))$anem_id, df$anem_2016[i])
#' #         df$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2018))$anem_id, df$anem_2017[i])
#' #       } else if(year_tag == 2015) {
#' #         df$anem_2015[i] = (fish.AllInfo %>% filter(fish_id == id, year == year_tag))$anem_id[1]
#' #         df$anem_2016[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2016))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2016))$anem_id, df$anem_2015[i])
#' #         df$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2017))$anem_id, df$anem_2016[i])
#' #         df$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2018))$anem_id, df$anem_2017[i])
#' #       } else if(year_tag == 2016) {
#' #         df$anem_2016[i] = (fish.AllInfo %>% filter(fish_id == id, year == year_tag))$anem_id[1]
#' #         df$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2017))$anem_id, df$anem_2016[i])
#' #         df$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(fish_id == tag, year == 2018))$anem_id, df$anem_2017[i])
#' #       } else if (year_tag ==)
#' #       
#' #         
#' #         
#' #         
#' #         df$anem_2016[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2016))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2016))$anem_id, fish.Tagged$anem_2015[i])
#' #         df$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id, fish.Tagged$anem_2016[i])
#' #         df$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(fish_id == id, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
#' #         
#' #       }
#' #       if(year_tag == 2015){
#' #         if(length((fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id) > 1) { print(paste("Multiple matching anems for fish",tag,"in year",year_tag))}
#' #         fish.Tagged$anem_2015[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#' #         fish.Tagged$anem_2016[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2016))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2016))$anem_id, fish.Tagged$anem_2015[i])
#' #         fish.Tagged$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id, fish.Tagged$anem_2016[i])
#' #         fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
#' #       } else if(year_tag == 2016){
#' #         fish.Tagged$anem_2016[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#' #         fish.Tagged$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id, fish.Tagged$anem_2016[i])
#' #         fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
#' #       } else if(year_tag == 2017){
#' #         fish.Tagged$anem_2017[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#' #         fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
#' #       } else if(year_tag == 2018){
#' #         fish.Tagged$anem_2018[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#' #       }
#' #       
#' #     }
#' # 
#' #     
#' #   } elseif(genortag == "tag") {
#' #     
#' #   }
#' #   
#' # }
#' # for(i in 1:length(dist_info$tag_id)) {
#' #   
#' #   if(year_tag == 2015){
#' #     if(length((fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id) > 1) { print(paste("Multiple matching anems for fish",tag,"in year",year_tag))}
#' #     fish.Tagged$anem_2015[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#' #     fish.Tagged$anem_2016[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2016))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2016))$anem_id, fish.Tagged$anem_2015[i])
#' #     fish.Tagged$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id, fish.Tagged$anem_2016[i])
#' #     fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
#' #   } else if(year_tag == 2016){
#' #     fish.Tagged$anem_2016[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#' #     fish.Tagged$anem_2017[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2017))$anem_id, fish.Tagged$anem_2016[i])
#' #     fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
#' #   } else if(year_tag == 2017){
#' #     fish.Tagged$anem_2017[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#' #     fish.Tagged$anem_2018[i] = ifelse(length((fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id) != 0, (fish.AllInfo %>% filter(tag_id == tag, year == 2018))$anem_id, fish.Tagged$anem_2017[i])
#' #   } else if(year_tag == 2018){
#' #     fish.Tagged$anem_2018[i] = (fish.AllInfo %>% filter(tag_id == tag, year == year_tag))$anem_id[1]
#' #   }
#' # }
