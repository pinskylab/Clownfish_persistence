# Looking at observed population size (and size-structure?) across sampling years - do pops seem to be persistent?

#################### Set-up: ####################
source(here::here('Code', 'Constants_database_common_functions.R'))

##### Load libraries
library(ggplot2)
library(cowplot)
library(lme4)

##### Load other output or data
# Load new proportion habitat sampled file (using proportion anems sampled rather than converting to area)
load(file=here::here("Data/Script_outputs", "anems_visited_by_year.RData"))

# Find average probability of capturing a fish (prob_r) from KC code
prob_r_avg <- mean(prob_r)

##### Set parameters 
sample_months <- c(3,4,5,6,7,8) #set months to avoid winter 2015 anem surveys

# Number of sites
n_sites = 19

# Make a site vec that has "Tamakin Dacot" rather than "Tomakin Dako" and one with alphabetical order to match the order of random effects in model fitting
site_vec_NS_TD <- site_vec_NS 
site_vec_NS_TD[17] <- "Tomakin Dako"

#################### Running things: ####################
##### Number of females (just from clownfish table for now - could change this later)
## Find raw number of females in database, now using sex column instead of color/size 
females_df_F <- allfish_caught %>%
  #filter(dive_type %in% dive_list) %>%  # not filtering by dives b/c most of the fish in 2017 were getting removed? but there are duplicates of tagged fish and such?? not sure how to deal with this...  filter(sex == "F") %>%
  filter(sex == "F") %>%
  group_by(year,site) %>%
  summarize(nFemalesRaw = n())

# ## Find raw number of females in database, using color/size
# females_df <- allfish_caught %>%
#   #filter(dive_type %in% dive_list) %>%  # not filtering by dives b/c most of the fish in 2017 were getting removed? but there are duplicates of tagged fish and such?? not sure how to deal with this...
#   filter(color %in% c('YP','Y')) %>%
#   filter(size >= min_female_size) %>%
#   group_by(year,site) %>%
#   summarize(nFemalesRaw = n())

## Scale by proportion of habitat sampled and probability of capturing a fish
females_df_F <- left_join(females_df_F, anems_visited_by_year %>%
                          filter(method == "metal tags") %>%
                          select(site, year, prop_hab_sampled_tidied), by = c('year', 'site'))

females_df_F <- females_df_F %>%
  mutate(prob_r = prob_r_avg) %>%
  mutate(nFemalesScaled = nFemalesRaw/(prob_r*prop_hab_sampled_tidied)) 

# For sites that are sampled but have Inf for scaled females b/c proportion hab sampled was 0, make the females estimated the raw # instead of NA
females_df_F <- females_df_F %>%
  mutate(nFemalesEstimated = ifelse(is.infinite(nFemalesScaled), nFemalesRaw, nFemalesScaled)) %>%
  mutate(nF = round(nFemalesEstimated))

# Find mean abundance at each site
females_df_F_mean <- females_df_F %>%
  group_by(site) %>%
  summarize(mean_abundance = mean(nFemalesEstimated, na.rm = TRUE))

# Set up data to fit models - seems to converge better if use 1-7 for year instead of 2012-2018, fix Tomakin Dako so shows up correctly in plots
abundance_mod_data <- females_df_F %>%
  mutate(year_order = year-2011) %>%
  mutate(site = if_else(site == "Tamakin Dacot", "Tomakin Dako", site)) %>%
  select(year, year_order, nF, site)

##### Fit mixed linear model
ab_mod <- glmer(nF ~ year_order + (year_order|site), data=abundance_mod_data, family=poisson)

# Find mean number of fish at site (intercept) and multiplier of the mean with time (slope), make a trend of mean F through time so easier to plot
site_trends <- data.frame(site = site_vec_order$site_name, mean_nF = exp(coef(ab_mod)$site[,1]), mean_multiplier = exp(coef(ab_mod)$site[,2]), stringsAsFactors = FALSE)

# Find line for average site
site_trends_all <- data.frame(site = "all", year=1:7) %>%
  mutate(mean_nF = exp(coef(summary(ab_mod))[1,1])*(exp(coef(summary(ab_mod))[2,1]))^year)

# Find line for individual sites
site_trends_time <- data.frame(site = site_vec_order$site_name[1], year=1:7, stringsAsFactors = FALSE) %>%
  mutate(mean_nF = (site_trends$mean_multiplier[1])^year*site_trends$mean_nF[1])

for(i in 2:n_sites) {
  trend_df <- data.frame(site = site_vec_order$site_name[i], year=1:7, stringsAsFactors = FALSE) %>%
    mutate(mean_nF = (site_trends$mean_multiplier[i])^year*site_trends$mean_nF[i])
  
  site_trends_time <- rbind(site_trends_time, trend_df)
}

# Fix Tomakin Dako name here too
site_trends_time <- site_trends_time %>%
  mutate(site = if_else(site == "Tamakin Dacot", "Tomakin Dako", site))

#################### Plots: ####################

##### Trend line of average site, plus individual sites in grey - no raw data
pdf(file=here::here("Plots/FigureDrafts","Abundance_through_time.pdf"))
ggplot(data = site_trends_time, aes(x=year, y=mean_nF, group=site)) +
  geom_line(color="grey") +
  geom_line(data=site_trends_all, aes(x=year, y=mean_nF), color = "black", size=1.5) +
  xlab("year") + ylab("# females") +
  scale_x_continuous(breaks=c(2,4,6), labels=c("2013","2015","2017"))
dev.off()

##### WSN plot: trend line of average site, plus individual sites in grey - no raw data
pdf(file=here::here("Plots/WSN_2019","WSN_2019_abundance_through_time.pdf"))
ggplot(data = site_trends_time, aes(x=year, y=mean_nF, group=site)) +
  geom_line(color="grey") +
  geom_line(data=site_trends_all, aes(x=year, y=mean_nF), color = "black", size=1.5) +
  xlab("year") + ylab("# females") +
  scale_x_continuous(breaks=c(2,4,6), labels=c("2013","2015","2017")) +
  theme_bw() +
  theme(text=element_text(size=25)) 
dev.off()

##### Multi-paneled figured with individual sites (doesn't currently include uncertainty/error bars but could...)
# Make a plot for each site
plot_list <- list()
for(i in 1:length(site_vec_NS_TD)) {  # works for the first sixteen, then has an issue with one of the final ones...
  site_i = site_vec_NS_TD[i]
  plot_list[[i]] <- ggplot(data = abundance_mod_data %>% filter(site == site_i), aes(x=year_order, y=nF)) +
    geom_point() +
    geom_line(data=site_trends_time %>% filter(site==site_i), aes(x=year, y=mean_nF)) +
    ggtitle(site_i) +
    ylab("# females") +  xlab("year") + 
    scale_x_continuous(breaks=c(2,4,6), labels=c("2013","2015","2017")) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
}

sites_together <- ggplot(data = site_trends_time, aes(x=year, y=mean_nF, group=site)) +
  geom_line(color="grey") +
  geom_line(data=site_trends_all, aes(x=year, y=mean_nF), color = "black", size=1.5) +
  xlab("year") + ylab("# females") +
  scale_x_continuous(breaks=c(2,4,6), labels=c("2013","2015","2017")) +
  theme_bw()

# And combine
pdf(file = here::here("Plots/FigureDrafts", "Time_series_scaled_F_by_site_with_lines.pdf"), height = 8.5, width = 10)
plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], 
          plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[11]],
          plot_list[[12]], plot_list[[13]], plot_list[[14]], plot_list[[15]], plot_list[[16]], plot_list[[17]],
          plot_list[[18]], plot_list[[19]], sites_together,
          labels = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t"), 
          nrow = 5)
dev.off()


#################### Saving output: ####################
save(females_df_F, file=here("Data/Script_outputs", "females_df_F.RData"))
save(site_trends_time, file=here::here("Data/Script_outputs", "site_trends_time.RData"))
save(site_trends_all, file=here::here("Data/Script_outputs", "site_trends_all.RData"))

#################### Old code: ####################
# Load in site_areas
# load(file=here::here('Data','site_areas.RData'))

# Load file with proportion habitat sampled estimates (from TotalAnemsAtSite.R)
#load(file=here("Data",'anem_sampling_table.RData')) #file with anems, after 2018 anems matched

# Goal - what does the pop look like over time?:
# Metric 1: number of females - at site level and overall
# Metric 2: number of males - at site level and overall
# Metric 3: size structure through time - at site level and overall

# # Join in whether or not a site was sampled in a particular year
# females_df <- left_join(site_visits, females_df, by = c('year', 'site'))
# females_df_F <- left_join(site_visits, females_df_F, by = c("year", "site"))

# # Fix sites that are sampled but with 0 proportion habitat sampled (S. Magbangon in 2017)
# females_df <- females_df %>%
#   mutate(nFemalesEstimated = ifelse(is.infinite(nFemalesScaled), NA, nFemalesScaled))
# 
# females_df_F <- females_df_F %>%
#   mutate(nFemalesEstimated = ifelse(is.infinite(nFemalesScaled), NA, nFemalesScaled))


###### Previous version of analysis and figure where I fit a linear model to each site separately
# # For time series, only consider sites with more than one sampling year (so not Caridad Proper, Sitio Lonas, or Sitio Tugas)
# site_vec_timeseries <- c("Cabatoan", "Caridad Cemetery", "Elementary School", "Gabas", "Haina", "Hicgop South",
#               "N. Magbangon", "S. Magbangon", "Palanas", "Poroc Rose", "Poroc San Flower", "San Agustin", "Sitio Baybayon", 
#               "Tamakin Dacot", "Visca", "Wangag")

# # Add in "all sites combined" 
# # For sites considered, find total proportion of habitat each makes up
# anems_visited_by_year_totalprophab <- anems_visited_by_year %>%
#   filter(site %in% site_vec_timeseries, method =  ) %>%
#   mutate(prop_total_hab_in_site = n_total_anems/sum(n_total_anems))
# 
# # Add an 'all sites' factor
# females_df_allsites <- females_df %>% 
#   group_by(year) %>%
#   summarize(prop_hab_sampled_as = sum(sampled))
#             prop_hab_sampled = )

# # Fit linear models for each site and for overall - here for fish identified as F by size and color
# females_df_models <- data.frame(site = site_vec_timeseries, stringsAsFactors = FALSE) %>%
#   mutate(intercept = NA,
#          coeff =  NA,
#          intercept_se = NA,
#          coeff_se = NA)
# 
# for(i in 1:length(site_vec_timeseries)) {
#   site_val = site_vec_timeseries[i]
#   df <- females_df %>% 
#     filter(site == site_val) %>%
#     filter(!is.na(nFemalesEstimated))
#   
#   testm <- lm(nFemalesEstimated ~ year, data=df)
#   testm_S <- summary(testm)
#   
#   females_df_models$intercept[i] = testm$coefficients[1]
#   females_df_models$coeff[i] = testm$coefficients[2]
#   females_df_models$intercept_se[i] = testm_S$coefficients[1,2]
#   females_df_models$coeff_se[i] = testm_S$coefficients[2,2]
# }
# 
# # Fit linear models for each site and for overall - here for fish identified as F by sex column
# females_df_F_models <- data.frame(site = site_vec_timeseries, stringsAsFactors = FALSE) %>%
#   mutate(intercept = NA,
#          coeff =  NA,
#          intercept_se = NA,
#          coeff_se = NA)
# 
# for(i in 1:length(site_vec_timeseries)) {
#   site_val = site_vec_timeseries[i]
#   df <- females_df_F %>% 
#     filter(site == site_val) %>%
#     filter(!is.na(nFemalesEstimated))
#   
#   testm <- lm(nFemalesEstimated ~ year, data=df)
#   testm_S <- summary(testm)
#   
#   females_df_F_models$intercept[i] = testm$coefficients[1]
#   females_df_F_models$coeff[i] = testm$coefficients[2]
#   females_df_F_models$intercept_se[i] = testm_S$coefficients[1,2]
#   females_df_F_models$coeff_se[i] = testm_S$coefficients[2,2]
# }
# 
# # Change Tamakin Dacot name spelling
# females_df$site <- replace(females_df$site, females_df$site=="Tamakin Dacot", "Tomakin Dako")
# females_df_F$site <- replace(females_df_F$site, females_df_F$site=="Tamakin Dacot", "Tomakin Dako")
# females_df_models$site <- replace(females_df_models$site, females_df_models$site=="Tamakin Dacot", "Tomakin Dako")
# females_df_F_models$site <- replace(females_df_F_models$site, females_df_F_models$site=="Tamakin Dacot", "Tomakin Dako")
# 




##### Previous attempts at fitting mixed linear model
# abundance_model <- glmer(nF ~ year + (year|site), data = females_df_F, family = poisson)
# abundance_model_a <- glmer(nF ~ year + (year|site), data = females_df_F, family = poisson, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
# abundance_model2 <- glmer(nF ~ year + (year|site), data = females_df, family = poisson)
# abundance_model3 <- glmer(nF ~ year, data = females_df_F, family = poisson)
# 
# 
# # Try just with sites that have multiple years, just to see how that goes...
# females_df_F_subset <- data.frame(abundance_mod_data_complete %>% filter(site %in% site_vec_timeseries))
# 
# am1 <- glmer(nF ~ year_order + (year_order|site), data = females_df_F_subset, family = poisson)
# am2 <- glmer(nF ~ year_order + (year_order|site), data=abundance_mod_data_complete, family=poisson)

# # sites_together <- ggplot(data = abundance_mod_data_complete, aes(x=year_order, y=nF, color=site)) +
# #   geom_point() +
# #   geom_line(data=site_trends_time, aes(x=year, y=mean_nF, group = site), color="grey")
# #   #geom_abline(slope=exp(coef(summary(ab_mod))[2,1]), intercept=exp(coef(summary(ab_mod))[1,1]))
# 
# years_order <- 1:7
# 
# for(i in 1:n_sites) {
#   sites_together <- sites_together +
#     #geom_abline(slope=exp(coef(ab_mod)$site[i,2]), intercept=exp(coef(ab_mod)$site[i,1]), color = "grey")
#     geom_line(x=1:7, y=exp(coef(ab_mod)$site[i,1])*year_order*exp(coef(ab_mod)$site[i,2]))
# }
# 
# # Example for Hicgop South, where points are clearly down but line is up
# # Think that's because of a misunderstanding of the slope - the slope is saying, how much does mean nF at Hicgop South change with increased time
# # But how can you get a negative effect with a Poisson if exp() can only give positive?
# ggplot(data = abundance_mod_data_complete %>% filter(site == "Hicgop South"), aes(x=year_order, y=nF)) +
#   geom_point() +
#   geom_abline(slope=exp(coef(ab_mod)$site[7,2]), intercept=exp(coef(ab_mod)$site[7,1]), color = "grey")


# ##### All sites together
# all_sites <- ggplot(data = abundance_mod_data_complete, aes(x=year_order, y=nF, color=site)) +
#   geom_point() +
#   geom_abline(slope=bY_est_overall, intercept=a_est_overall)
# 
# for(i in 1:n_sites) {
#   all_sites <- all_sites +
#     geom_abline(slope=site_trends$bY_est[i], intercept=site_trends$a_est, color = "grey")
# }



##### Attempt at fitting mixed effects model using rethinking
# ####### Try using rethinking to fit

# # # This one doesn't fit very well...
# # abundance_mod_m1 <- map2stan(
# #   alist(
# #     nF ~ dpois(lambda),
# #     log(lambda) <- a + bY*year,
# #     a ~ dnorm(0,100),
# #     bY ~ dnorm(0,5)
# #   ), data=abundance_mod_data_complete)
# 
# # Adding in varying intercept and slope by site
# abundance_mod_m2 <- map2stan(
#   alist(
#     nF ~ dpois(lambda),
#     log(lambda) <- a_site[site] + bY_site[site]*year_order,
#     c(a_site, bY_site)[site] ~ dmvnorm2(c(a,bY), sigma_site, Rho),
#     a ~ dnorm(0,10),
#     bY ~ dnorm(0,10),
#     sigma_site ~ dcauchy(0,2),
#     Rho ~ dlkjcorr(2)
#   ), data=abundance_mod_data_complete, warmup=1000, iter=6000, chains=2)
# 
# # Extract samples
# post <- extract.samples(abundance_mod_m2)
# 
# # Trend by site - I'm not sure these parameters are what I think they are...
# site_trends <- data.frame(site_no = 1:19, a_est = NA, bY_est = NA)
# a_est_overall = exp(mean(post$a))
# bY_est_overall = exp(mean(post$bY))
# 
# for(i in 1:n_sites) {
#   site_trends$a_est[i] = exp(mean(post$a + post$a_site[i]))
#   site_trends$bY_est[i] = exp(mean(post$bY + post$bY_site[i]))
# }

#################### Plots: ####################

# # Using cowplot to make site plots individually so y-axis scale can vary by site
# 
# Cabatoan_F <- ggplot(data = females_df_F %>% filter(site == "Cabatoan"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Cabatoan"), aes(intercept = intercept, slope = coeff)) +
#   ggtitle('Cabatoan') +
#   ylab('# scaled females') +  xlab('year') + 
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# CaridadCemetery_F <- ggplot(data = females_df_F %>% filter(site == "Caridad Cemetery"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Caridad Cemetery"), aes(intercept = intercept, slope = coeff)) +
#   ylab('# scaled females') + ggtitle('Caridad Cemetery') +
#   xlab('year') +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# ElementarySchool_F <- ggplot(data = females_df_F %>% filter(site == "Elementary School"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Elementary School"), aes(intercept = intercept, slope = coeff)) +
#   ggtitle('Elementary School') +
#   xlab('year') + ylab('# scaled females') + 
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# Gabas_F <- ggplot(data = females_df_F %>% filter(site == "Gabas"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Gabas"), aes(intercept = intercept, slope = coeff)) +
#   xlab('year') + ggtitle('Gabas') +
#   ylab('# scaled females') +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# Haina_F <- ggplot(data = females_df_F %>% filter(site == "Haina"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Haina"), aes(intercept = intercept, slope = coeff)) +
#   xlab('year') + ggtitle("Haina") +
#   ylab('# scaled females') +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# HicgopSouth_F <- ggplot(data = females_df_F %>% filter(site == "Hicgop South"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Hicgop South"), aes(intercept = intercept, slope = coeff)) +
#   ggtitle('Hicgop South') +
#   xlab('year') + ylab('# scaled females') +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# NMagbangon_F <- ggplot(data = females_df_F %>% filter(site == "N. Magbangon"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "N. Magbangon"), aes(intercept = intercept, slope = coeff)) +
#   ggtitle('N. Magbangon') +
#   xlab('year') + ylab('# scaled females') + 
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# Palanas_F <- ggplot(data = females_df_F %>% filter(site == "Palanas"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Palanas"), aes(intercept = intercept, slope = coeff)) +
#   ylab('# scaled females') + ggtitle('Palanas') +
#   xlab('year') + 
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# PorocRose_F <- ggplot(data = females_df_F %>% filter(site == "Poroc Rose"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Poroc Rose"), aes(intercept = intercept, slope = coeff)) +
#   xlab('year') + ylab('# scaled females') + ggtitle('Poroc Rose') +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# PorocSanFlower_F <- ggplot(data = females_df_F %>% filter(site == "Poroc San Flower"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Poroc San Flower"), aes(intercept = intercept, slope = coeff)) +
#   ggtitle('Poroc San Flower') +
#   xlab('year') + ylab('# scaled females') + 
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# SMagbangon_F <- ggplot(data = females_df_F %>% filter(site == "S. Magbangon"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "S. Magbangon"), aes(intercept = intercept, slope = coeff)) +
#   ggtitle('S. Magbangon') +
#   xlab('year') + ylab('# scaled females') + 
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# SanAgustin_F <- ggplot(data = females_df_F %>% filter(site == "San Agustin"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "San Agustin"), aes(intercept = intercept, slope = coeff)) +
#   ggtitle('San Agustin') +
#   xlab('year') + ylab('# scaled females') + 
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# SitioBaybayon_F <- ggplot(data = females_df_F %>% filter(site == "Sitio Baybayon"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Sitio Baybayon"), aes(intercept = intercept, slope = coeff)) +
#   xlab('year') + ylab('# scaled females') + ggtitle('Sitio Baybayon') +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# TomakinDako_F <- ggplot(data = females_df_F %>% filter(site == "Tomakin Dako"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Tomakin Dako"), aes(intercept = intercept, slope = coeff)) +
#   xlab('year') + ggtitle('Tomakin Dako') +
#   ylab('# scaled females') + 
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# Visca_F <- ggplot(data = females_df_F %>% filter(site == "Visca"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Visca"), aes(intercept = intercept, slope = coeff)) +
#   xlab('year') + ggtitle("Visca") +
#   ylab('# scaled females') + 
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# Wangag_F <- ggplot(data = females_df_F %>% filter(site == "Wangag"), aes(x=year, y=nFemalesEstimated)) +
#   geom_point() +
#   geom_abline(data = females_df_F_models %>% filter(site == "Wangag"), aes(intercept = intercept, slope = coeff)) +
#   ggtitle('Wangag') +
#   xlab('year') + ylab('# scaled females') + 
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# # Histogram of slopes
# slopes_F_hist <- ggplot(data = females_df_F_models, aes(x=coeff)) +
#   geom_histogram(binwidth = 0.5) +
#   xlab("slope") + ggtitle("Histogram of slopes") +
#   theme_bw()
# 
# # ggplot(data = females_df_models, aes(x=coeff)) +
# #   geom_histogram(binwidth = 0.5) 
# 
# # And combine
# pdf(file = here::here("Plots/FigureDrafts", "Time_series_scaled_F_by_site_with_lines.pdf"), height = 8.5, width = 10)
# plot_grid(Palanas_F, Wangag_F, NMagbangon_F, SMagbangon_F, Cabatoan_F, CaridadCemetery_F, 
#           HicgopSouth_F, ElementarySchool_F, SanAgustin_F, PorocSanFlower_F, PorocRose_F,
#           Visca_F, Gabas_F, TomakinDako_F, Haina_F, SitioBaybayon_F, slopes_F_hist, 
#           labels = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q"), 
#           nrow = 4)
# dev.off()
# 

# # Scaled-up number of females at each site through time
# pdf(file = here('Plots/DataCharacteristics', 'Scaled_F_through_time_by_site.pdf'))
# ggplot(data = females_df, aes(x=year, y=nFemalesEstimated)) +
#   #geom_bar(stat='identity') +
#   geom_point() +
#   geom_abline(data = females_df_models, aes(intercept = intercept, slope = coeff)) +
#   #geom_ribbon(data = females_df_models, aes(ymin = ))
#   xlab('year') + ylab('# scaled females') + ggtitle('Females by site through time') +
#   facet_wrap(~site) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# dev.off()
# 
# # Scaled-up number of females at each site through time - attempts at figures
# pdf(file = here('Plots/FigureDrafts', 'Time_series_scaled_F_by_site_with_lines.pdf'))
# ggplot(data = females_df %>% filter(site %in% site_vec_timeseries), aes(x=year, y=nFemalesEstimated)) +
#   #geom_bar(stat='identity') +
#   geom_point() +
#   geom_abline(data = females_df_models, aes(intercept = intercept, slope = coeff)) +
#   #geom_ribbon(data = females_df_models, aes(ymin = ))
#   xlab('year') + ylab('# scaled females') + #ggtitle('Females by site through time') +
#   facet_wrap(~site) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# dev.off()
# 
# # Scaled-up number of females at each site through time - attempts at figures
# pdf(file = here('Plots/FigureDrafts', 'Time_series_scaled_F_by_site_no_lines.pdf'))
# ggplot(data = females_df, aes(x=year, y=nFemalesEstimated)) +
#   #geom_bar(stat='identity') +
#   geom_point() +
#   #geom_abline(data = females_df_models, aes(intercept = intercept, slope = coeff)) +
#   #geom_ribbon(data = females_df_models, aes(ymin = ))
#   xlab('year') + ylab('# scaled females') + ggtitle('Females by site through time') +
#   facet_wrap(~site) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# dev.off()


#################### Code from elsewhere that might be useful: ####################
# ##### Find number of breeding females
# # Find number of breeding females seen or captured
# breedingF <- allAPCL %>% 
#   group_by(site,year) %>%
#   filter(color %in% c("Y","YP")) %>% #just filter out YP fish (and Y for 2012, when there was no YP)
#   filter(size >= min_size_breedingF) %>% #but also make sure they're big enough to be a breeder (sometimes small fish get marked YP too)
#   summarise(NbreedingF = n())
# 
# # Find number of breeding females tagged
# breedingF_tagged <- allAPCL %>%
#   group_by(year,site) %>%
#   filter(color %in% c("Y","YP")) %>%
#   filter(!is.na(tag_id)) %>%
#   summarise(Ntagged_breedingF = n()) 
# 
# # Combine - since didn't tag all years, might not have any females for some site/years
# breedingF_info <- left_join(breedingF, breedingF_tagged, by = c("site" = "site", "year" = "year")) # no 2012 obs here - why?
# 
# # Find average prob_r from KC code
# prob_r_avg <- mean(prob_r)
# 
# 
# demog_info_recruits <- recruits_info %>%
#   group_by(site) %>%
#   summarize(mean_est_R = mean(totalR_est_metalTA, na.rm=TRUE),
#             high_est_R = max(totalR_est_metalTA, na.rm=TRUE),
#             low_est_R = min(totalR_est_metalTA, na.rm=TRUE),
#             mean_raw_R = mean(Nrecruits, na.rm=TRUE),
#             high_raw_R = max(Nrecruits, na.rm=TRUE),
#             low_raw_R = min(Nrecruits, na.rm=TRUE))
# 
# demog_info_eggs <- breedingF_info %>%
#   group_by(site) %>%
#   summarize(mean_est_eggs = mean(est_eggs_metalTA, na.rm=TRUE),
#             high_est_eggs = max(est_eggs_metalTA, na.rm=TRUE),
#             low_est_eggs = min(est_eggs_metalTA, na.rm=TRUE),
#             mean_raw_F = mean(NbreedingF_combo, na.rm=TRUE),
#             high_raw_F = max(NbreedingF_combo, na.rm=TRUE),
#             low_raw_F = min(NbreedingF_combo, na.rm=TRUE))
# 
