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
  filter(sex == "F") %>%
  group_by(year,site) %>%
  summarize(nFemalesRaw = n())

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

# Find average number of females at each site
average_females_by_site <- site_trends_time %>% 
  group_by(site) %>%
  summarize(avg_nF = mean(mean_nF))

##### Does population summed across patches stay stable?
# Sites sampled every year
all_years_2012_to_2018 <- c("Palanas", "Wangag", "San Agustin", "Visca")
most_years_2013_to_2018 <- c(all_years_2012_to_2018, "Elementary School", "Tamakin Dacot", "Haina", "Sitio Baybayon")

# Sum
females_df_F_summed_all <- females_df_F %>%
  filter(site %in% all_years_2012_to_2018) %>%
  group_by(year) %>%
  summarize(nF_metapop = sum(nF))

females_df_F_summed_most <- females_df_F %>%
  filter(site %in% most_years_2013_to_2018) %>%
  group_by(year) %>%
  summarize(nF_metapop = sum(nF)) 

females_df_F_summed_most <- rbind(females_df_F_summed_most, data.frame(year = 2012, nF_metapop = NA))

#################### Plots: ####################

##### Trend line of average site, plus individual sites in grey - no raw data
pdf(file=here::here("Plots/FigureDrafts","Abundance_through_time.pdf"))
ggplot(data = site_trends_time, aes(x=year, y=mean_nF, group=site)) +
  geom_line(color="grey") +
  geom_line(data=site_trends_all, aes(x=year, y=mean_nF), color = "black", size=1.5) +
  xlab("Year") + ylab("# Females") +
  scale_x_continuous(breaks=c(2,4,6), labels=c("2013","2015","2017"))
dev.off()

##### Summed across patches (for those sampled all or most years)
pdf(file=here::here("Plots/FigureDrafts", "Summed_abundance_through_time.pdf"))
ggplot(data =  females_df_F_summed_all, aes(x=year, y=nF_metapop)) +
  geom_line(color = "blue") +
  geom_line(data=females_df_F_summed_most, aes(x=year, y=nF_metapop), color = "orange") +
  xlab("year") + ylab("# summed females") +
  theme_bw()
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
    ylab("# Females") +  xlab("Year") + 
    scale_x_continuous(breaks=c(2,4,6), labels=c("2013","2015","2017")) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
}

sites_together <- ggplot(data = site_trends_time, aes(x=year, y=mean_nF, group=site)) +
  geom_line(color="grey") +
  geom_line(data=site_trends_all, aes(x=year, y=mean_nF), color = "black", size=1.5) +
  xlab("Year") + ylab("# Females") +
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

# For correct journal format
FigD6_plot <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], 
                        plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[11]],
                        plot_list[[12]], plot_list[[13]], plot_list[[14]], plot_list[[15]], plot_list[[16]], plot_list[[17]],
                        plot_list[[18]], plot_list[[19]], sites_together,
                        labels = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t"), 
                        nrow = 5)

ggplot2::ggsave(filename = here::here("For_submission/Final_submission", "FigD6.pdf"), plot = FigD6_plot, scale=1, width=8.5, height=10, units = "in",
                dpi=1000)

#################### Saving output: ####################
save(females_df_F, file=here("Data/Script_outputs", "females_df_F.RData"))
save(site_trends_time, file=here::here("Data/Script_outputs", "site_trends_time.RData"))
save(site_trends_all, file=here::here("Data/Script_outputs", "site_trends_all.RData"))

