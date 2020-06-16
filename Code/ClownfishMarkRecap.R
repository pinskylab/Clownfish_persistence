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
load(file = here::here("Data/Script_outputs", "encounters_size_0.RData"))  # encounter histories with sizes and missing ones filled in with projections or 0s
load(file = here::here("Data/Script_outputs", "encounters_size_means_by_year.RData"))  # encounter histories with sizes and missing ones filled in with projections or by-year means

# Size vector for plotting
n_size_steps <- 30
size_values <- min_size+(0:n_size_steps)*(max_size-min_size)/n_size_steps  # min and max size are in Constants_database_common_functions.R

# Distance vector for plotting
n_dist_steps <- 100
min_dist = 0
max_dist = 200  # actual max dist is very large, much larger than distance where recapture prob goes to essentially zero, so use this instead
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
encounters_size_means_by_year <- encounters_size_means_by_year %>% filter(!site %in% one_time_sites)

##### Create data frames for analysis, joining trait info and distance info
# Mean size (by year), mean distance (by year) filled in for fish pre-first capture
eall_meanYsize_meanYdist <- left_join(encounters_size_means_by_year, encounters_dist_mean_by_year %>% select(fish_indiv, dist2013, dist2014, dist2015, dist2016, dist2017, dist2018),
                                      by = "fish_indiv") %>%
  mutate(site = as.factor(site), cap_color = as.factor(cap_color), cap_stage = as.factor(cap_stage)) %>%
  dplyr::rename(ch = encounter.hist)
  
eall_meanYsize_meanYdist <- eall_meanYsize_meanYdist[complete.cases(eall_meanYsize_meanYdist),]

eall_meanYsize_meanYdist <- eall_meanYsize_meanYdist %>%  # MARK doesn't handle spaces in site names well so remove them
  mutate(site = case_when(site == "Cabatoan" ~ "Cabatoan",
                          site == "Caridad Cemetery" ~ "CaridadCemetery",
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

