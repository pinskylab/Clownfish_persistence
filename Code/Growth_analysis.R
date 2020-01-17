# Fitting a growth curve (pulling code from growth_newdata.R in Growth repository, figure out how to credit Michelle correctly)

#################### Set-up: ####################
source(here::here("Code", "Constants_database_common_functions.R"))

# These are already loaded in Constants_database_common_functions.R
# library(dplyr)
# library(tidyr)
# library(lubridate)
# library(stringr)

# library(fishmethods) - # this is the library Michelle used, could go back and do that if want to make models more complicated...
# library(readr)

library(ggplot2)

#library(grid)
#library(gridExtra)

#library(mgcv)  # used by Rees et al. 2014
#library(lme4)  # used by Hart et al. 2009

##### Load files from other scripts within this repository or source those scripts (below, commented out)

##### Set up parameters

# Pull out a set where fish were captured about a year apart
year_lower_limit = 345
year_upper_limit = 385

month_days = 30 

#################### Functions: ####################

# Function to find the first size, second size, and time elapsed for a pair of encounters 
findRecapPair <- function(recaps_df, marked_fish_df, encounter1, encounter2) {
  out <- marked_fish_df %>%
    filter(fish_id %in% recaps_df$fish_id) %>% 
    group_by(fish_id) %>%
    arrange(date) %>%
    summarize(L1 = size[encounter1],
              L2 = size[encounter2],
              tal = as.Date(date[encounter2]) - as.Date(date[encounter1]))
  return(out)
}

# Choose one recap pair per fish, run lm, save coefficients, pairs already selected for particular time limit
chooseRecapFitModel_TimeLimitPreFiltered <- function(n_runs, pairs_df, multiple_year_pairs) {  # pairs df is something like recap_pairs_year with all pairs caught within the time frame, multiple_year_pairs is like multiple_year_recap_fish filtered for n_caps > 1
  
  # Set up an output data frame
  out = data.frame(run = seq(from=1, to=n_runs, by=1),
                   intercept_est = NA,
                   slope_est = NA,
                   intercept_se = NA,
                   slope_se = NA)
  
  # Pull out the pairs that were only caught once
  pairs_1_set <- pairs_df %>% 
    filter(fish_id %in% ((multiple_year_pairs %>% filter(n_year_caps == 1))$fish_id))
  
  # In for each run, choose pairs, run model, save output
  for(i in 1:n_runs) {
    
    # Choose a pair for the fish caught multiple times
    set_selected <- multiple_year_pairs %>%
      filter(n_year_caps > 1) %>%
      mutate(pair_chosen = NA)
    
    # Create a data frame to collect the chosen pairs from the fish with 2 or more (this is a little weird, should figure out a better way...)
    pairs_2plus_set <- data.frame(fish_id = NA, L1 = NA, L2 = NA,
                                  tal_days = NA, tal = NA, growth = NA,
                                  growth_per_year = NA, L2_one_year_later = NA)
    
    for(j in 1:length(set_selected$fish_id)) {
      
      # select a pair
      set_selected$pair_chosen[j] = base::sample(1:set_selected$n_year_caps[j], 1)
      
      # pull out the info for the pairs for that fish
      pair_group <- pairs_df %>%
        filter(fish_id == set_selected$fish_id[j])
      
      # add the selected to the data frame
      pairs_2plus_set <- rbind(pairs_2plus_set, pair_group[set_selected$pair_chosen[j],])
      
    }
    
    # Remove that extra NA line at the top, join pair sets togethe
    all_pairs <- rbind(pairs_1_set, pairs_2plus_set[-1,])
    
    # Fit model, save output
    model_out <- lm(L2 ~ L1, data = all_pairs)
    
    out$intercept_est[i] = coef(summary(model_out))["(Intercept)", "Estimate"]
    out$slope_est[i] = coef(summary(model_out))["L1", "Estimate"]
    out$intercept_se[i] = coef(summary(model_out))["(Intercept)", "Std. Error"]
    out$slope_se[i] = coef(summary(model_out))["L1", "Std. Error"]
  }
  return(out)
}

#################### Running things: ####################

########## Finding recaptures and prepping the data ##########

##### Assign each individual fish marked and recaptured a fish_id (since one fish can have multiple tag_ids and a gen_id) - NOW JUST USING fish_indiv from Michelle code
# Pull out all fish that have been marked in some way, do a first-pass at assigning a fish_id (this and below is from the updated way of connecting fish, from the growth_newdata.R script in Growth repo)
# marked_fish <- allfish_caught %>%
#   dplyr::select(sample_id, color, size, sex, tag_id, gen_id, fish_table_id, anem_table_id, dive_table_id, date, site, year, fish_notes) %>%
#   filter(!is.na(gen_id) | (tag_id != 'NA' & !is.na(tag_id))) %>% # pull out fish "tagged" in any way, either PIT or via genetic sample
#   mutate(fish_id = case_when(!is.na(gen_id) ~ paste('gen', gen_id, sep=''),  # if fish has a gen_id, use that as the fish_id
#                              is.na(gen_id) ~ paste('tag', tag_id, sep='')))  # if it doesn't have a gen_id, use the tag_id
# 
# # Go through and match up fish ided in mulitple ways (now that gen_id is more stable than tag_id - because fish lose their tags - use that as the primary)
# # one tagged fish has two gen_ids - see GitHub issue, can directly assign id to that fish
# for(i in 1:length(marked_fish$fish_id)) {
#   if(substring(marked_fish$fish_id[i],1,3) == 'tag') {  # if it has a fish_id based on tag_id rather than gen_id...
#     tag_id_val <- marked_fish$tag_id[i]  # pull out the tag_id
#     
#     matches <- marked_fish %>%
#       filter(tag_id == tag_id_val)  # filter out any other cases that match that tag_id
#     
#     for(j in 1:length(matches$fish_id)) {
#       if(!is.na(matches$gen_id[j])) {  # if any of the captures with that tag_id have a gen_id
#         marked_fish$fish_id[i] <- matches$fish_id[j]  # update the fish_id to the fish_id based on the match with a gen_id - WHAT IF HAS MULTIPLE GEN_IDs?
#       }
#     }
#   }
# }

# Filter out just marked fish, using fish_indiv
marked_fish <- allfish_caught %>%  # when running on 
  filter(!is.na(fish_indiv))  # only use fish that have been "marked" in some way and assigned a fish_indiv id in fish_obs table generated by Michelle

# Make fish_indiv called fish_id so can use with rest of code (will change eventually)
marked_fish <- marked_fish %>%
  dplyr::rename(fish_id = fish_indiv)

# Remove any observations where size is not recorded - goes from XX (3321 distinct fish) to 4012 obs  (old info: goes from 4001 observation to 3994 observations (pre-2016-2018 genetic data numbers: goes from 3596 observations to 3592)
marked_fish <- marked_fish %>%  # 3053 individual fish on 10/15/19
  filter(size != "NA" & !is.na(size))

# Remove any observations where the date is not recorded (doesn't remove any)
marked_fish <- marked_fish %>% filter(!is.na(date))

##### Pull out recaptured fish, create recapture pairs
# Pull out fish captured more than once - 503 obs (old: 464 obs)
recap_fish <- marked_fish %>%
  group_by(fish_id) %>%
  summarize(ncaps = n(),
            year1 = min(year)) %>%
  filter(ncaps > 1)

# Find max number of recaptures for a fish
max_recaps <- max(recap_fish$ncaps)

# Pull out fish captured two through 8 (max_recaps) times
recaps_2X <- recap_fish %>% filter(ncaps == 2)
recaps_3X <- recap_fish %>% filter(ncaps == 3)
recaps_4X <- recap_fish %>% filter(ncaps == 4)
recaps_5X <- recap_fish %>% filter(ncaps == 5)
recaps_6X <- recap_fish %>% filter(ncaps == 6)
recaps_7X <- recap_fish %>% filter(ncaps == 7)
recaps_8X <- recap_fish %>% filter(ncaps == 8)

# Pull out sizes and time lapsed for the the various recapture pairs
recaps_2X_pair1 = findRecapPair(recaps_2X, marked_fish, 1, 2) #recap pair for fish caught 2X
recaps_3X_pair1 = findRecapPair(recaps_3X, marked_fish, 1, 2) #first recap pair for fish caught 3X (first and second encounters)
recaps_3X_pair2 = findRecapPair(recaps_3X, marked_fish, 2, 3) #second recap pair fish caught 3X (second and third encounters)
recaps_4X_pair1 = findRecapPair(recaps_4X, marked_fish, 1, 2) #first recap pair for fish caught 4X (first and second encounters)
recaps_4X_pair2 = findRecapPair(recaps_4X, marked_fish, 2, 3) #second recap pair for fish caught 4X (second and third encounters)
recaps_4X_pair3 = findRecapPair(recaps_4X, marked_fish, 3, 4) #third recap pair for fish caught 4X (third and fourth encounters)
recaps_5X_pair1 = findRecapPair(recaps_5X, marked_fish, 1, 2) #first and second encounters for fish caught 5X
recaps_5X_pair2 = findRecapPair(recaps_5X, marked_fish, 2, 3) #second and third encounters for fish caught 5X
recaps_5X_pair3 = findRecapPair(recaps_5X, marked_fish, 3, 4) #third and fourth encounters for fish caught 5X
recaps_5X_pair4 = findRecapPair(recaps_5X, marked_fish, 4, 5) #fourth and fifth encounters for fish caught 5X
recaps_6X_pair1 = findRecapPair(recaps_6X, marked_fish, 1, 2) #first and second encouters for fish caught 6X
recaps_6X_pair2 = findRecapPair(recaps_6X, marked_fish, 2, 3) #second and third encounters for fish caught 6X
recaps_6X_pair3 = findRecapPair(recaps_6X, marked_fish, 3, 4) #third and fourth encounters for fish caught 6X
recaps_6X_pair4 = findRecapPair(recaps_6X, marked_fish, 4, 5) #fourth and fifth encounters for fish caught 6X
recaps_6X_pair5 = findRecapPair(recaps_6X, marked_fish, 5, 6) #fifth and sixth encounters for fish caught 6X
recaps_7X_pair1 = findRecapPair(recaps_7X, marked_fish, 1, 2) #first and second encouters for fish caught 7X
recaps_7X_pair2 = findRecapPair(recaps_7X, marked_fish, 2, 3) #second and third encounters for fish caught 7X
recaps_7X_pair3 = findRecapPair(recaps_7X, marked_fish, 3, 4) #third and fourth encounters for fish caught 7X
recaps_7X_pair4 = findRecapPair(recaps_7X, marked_fish, 4, 5) #fourth and fifth encounters for fish caught 7X
recaps_7X_pair5 = findRecapPair(recaps_7X, marked_fish, 5, 6) #fifth and sixth encounters for fish caught 7X
recaps_7X_pair6 = findRecapPair(recaps_7X, marked_fish, 6, 7) #sixth and seventh encounters for fish caught 7X
recaps_8X_pair1 = findRecapPair(recaps_8X, marked_fish, 1, 2) #first and second encouters for fish caught 8X
recaps_8X_pair2 = findRecapPair(recaps_8X, marked_fish, 2, 3) #second and third encounters for fish caught 8X
recaps_8X_pair3 = findRecapPair(recaps_8X, marked_fish, 3, 4) #third and fourth encounters for fish caught 8X
recaps_8X_pair4 = findRecapPair(recaps_8X, marked_fish, 4, 5) #fourth and fifth encounters for fish caught 8X
recaps_8X_pair5 = findRecapPair(recaps_8X, marked_fish, 5, 6) #fifth and sixth encounters for fish caught 8X
recaps_8X_pair6 = findRecapPair(recaps_8X, marked_fish, 6, 7) #sixth and seventh encounters for fish caught 8X
recaps_8X_pair7 = findRecapPair(recaps_8X, marked_fish, 7, 8) #seventh and eighth encounters for fish caught 8X

# Join pairs together into one data frame - 697 pairs (earlier was 772 pairs? - maybe PIT typos were fixed?) (was 545 pairs before addition of 2016-2018 genetic data) 
recap_pairs <- rbind(recaps_2X_pair1, recaps_3X_pair1, recaps_3X_pair2,
                     recaps_4X_pair1, recaps_4X_pair2, recaps_4X_pair3,
                     recaps_5X_pair1, recaps_5X_pair2, recaps_5X_pair3, recaps_5X_pair4,
                     recaps_6X_pair1, recaps_6X_pair2, recaps_6X_pair3, recaps_6X_pair4, recaps_6X_pair5,
                     recaps_7X_pair1, recaps_7X_pair2, recaps_7X_pair3, 
                     recaps_7X_pair4, recaps_7X_pair5, recaps_7X_pair6,
                     recaps_8X_pair1, recaps_8X_pair2, recaps_8X_pair3,
                     recaps_8X_pair4, recaps_8X_pair5, recaps_8X_pair6, recaps_8X_pair7)

# Filter out those captured on the same day
recap_pairs <- recap_pairs %>%
  filter(tal != 0)

# Add in percentage-of-year column for time 
recap_pairs <- recap_pairs %>%
  dplyr::rename(tal_days = tal) %>%
  mutate(tal = as.numeric(tal_days/365))

# Add in column for growth, growth per year, make size columns numeric, add in projected growth 1 year later
recap_pairs <- recap_pairs %>%
  mutate(L1 = as.numeric(L1),
         L2 = as.numeric(L2),
         growth = L2-L1,
         growth_per_year = growth/tal,
         L2_one_year_later = round(L1 + growth_per_year, digits = 1))  # projected growth one year later, based on growth-per-year

# Filter out just those recaptured about a year later - 244
recap_pairs_year <- recap_pairs %>%
  filter(tal_days > year_lower_limit & tal_days <= year_upper_limit) 

# Find those fish that have multiple recaptures that are about a year apart
multiple_year_recap_fish <- recap_pairs_year %>%
  group_by(fish_id) %>%
  summarize(n_year_caps = n())

########## Fitting growth models ##########

# Set a size range for making predictions from models (to plot)
size_range_sequence <- data.frame(L1 = seq(1, max_size, length = max_size*10))

##### Hart et al. 2009 method - just fit a basic linear regression (not enough multi-year fish with about a year in between to do a random mixed effects model)
# If just do a linear regression between Lt and Lt+1, can calculate K and Linf from y-intercept (b) and slope (m) - but doesn't take into account random effect or individual effects

### Regular linear model regression, no random effects or anything, using all fish recaptures within the time frame (including the 23 fish that have 2-3 recaptures that fit that)
model_1 <- lm(L2 ~ L1, data =  (recap_pairs_year))

model_1_k <- as.numeric(-log(model_1$coefficients[2]))  # 0.919 if use all fish captured about a year later (even multiple pairs per fish) --  0.881 one time, 0.938 another, pretty similar to fishmethods lowest AIC estimate of 0.867
model_1_Linf <- as.numeric(model_1$coefficients[1]/(1 - model_1$coefficients[2]))  # 10.61 if use all fish captured about a year later (even multiple pairs per fish) --  10.69 one time, 10.58 another, pretty similar to fishmethods lowest AIC estimate of 10.64

model_1_predictions <- predict.lm(model_1, newdata = size_range_sequence, se.fit = TRUE)

#### Run a bunch, choosing one pair for each multiple-times-recaptured scripts each time 
model_1_runs <- chooseRecapFitModel_TimeLimitPreFiltered(n_runs, recap_pairs_year, multiple_year_recap_fish)

model_1_runs <- model_1_runs %>%
  mutate(k_est = -log(slope_est),
         Linf_est = intercept_est/(1 - slope_est))

### Also tried a gam, with random effects by fish_id but really don't have that many fish with multiple long recaptures, didn't fit well
# If come across a way to use varying time intervals, random effects by fish_id, measurement error at all measurement points, and uncertainty in Linf and k, would be great but maybe something to do down the road...

##### Comparing growth rates of fish alone on an anemone to those with larger fish

##### Bump up best estimate of K

#################### Plots: ####################

# Plot histogram of Linf and k estimates for n_runs with choosing different recap pairs
Linf_plot <- ggplot(data = model_1_runs, aes(x = Linf_est)) +
  geom_histogram(bins = 40, color = "gray", fill = "gray") +
  geom_vline(xintercept = mean(model_1_runs$Linf_est), color = "black") +
  xlab("Linf (cm) estimate") + ggtitle("Linf: choosing pairs, 1 year") +
  theme_bw()

k_plot <- ggplot(data = model_1_runs, aes(x = k_est)) +
  geom_histogram(bins = 40, color = "gray", fill = "gray") +
  geom_vline(xintercept = mean(model_1_runs$k_est)) +
  xlab("k estimate") + ggtitle("k: choosing pairs, 1 year") +
  theme_bw()

pdf(file = here::here("Plots/Growth", "Linf_and_k_histogram_runs_choosing_pairs.pdf"), width=6, height=3)
grid.arrange(Linf_plot, k_plot, nrow=1)
dev.off()

# Plot the data and linear model for fish caught about a year apart, models for only one recapture pair per fish, pairs selected randomly and models fit 1000x
pdf(file = here::here("Plots/Growth", "Data_L1_L2_with_fit_line_range_all_one_year_recap_fish.pdf"))
ggplot(data = recap_pairs_year, aes(x = L1, y = L2)) +
  geom_point(shape = 1) +
  geom_abline(aes(intercept = mean(model_1_runs$intercept_est), slope = mean(model_1_runs$slope_est)), color = "black") +
  #geom_abline(aes(intercept = min(model_1_runs$intercept_est), slope = min(model_1_runs$slope_est)), color = "gray") +
  #geom_abline(aes(intercept = max(model_1_runs$intercept_est), slope = max(model_1_runs$slope_est)), color = "gray") +
  geom_ribbon(aes(x=seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1)),
                  ymin = (min(model_1_runs$intercept_est) + min(model_1_runs$slope_est)*seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1))),
                  ymax = (max(model_1_runs$intercept_est) + max(model_1_runs$slope_est)*seq(from=2.5, to=12.5,length.out = length(recap_pairs_year$L1)))), fill = "light gray", alpha = 0.5) +
  xlab("Length (cm) at year t") + ylab("Length (cm) at year t+1") +
  ggtitle("Model fits across runs with data") +
  theme_bw()
dev.off()

#################### Saving output: ####################
growth_info_estimate = model_1_runs  # saving it with different name in case change input later...
save(recap_pairs_year, file = here::here("Data/Script_outputs", "recap_pairs_year.RData"))  # saving all recapture pairs caught about a year apart, for plotting purposes
save(growth_info_estimate, file = here::here("Data/Script_outputs", "growth_info_estimate.RData"))

#################### Old code: ####################
# 
# # Choose just one pair for fish caught multiple times
# pair_3X = chooseRecapPair_3X(recaps_3X_pair1, recaps_3X_pair2)
# pair_4X = chooseRecapPair_4X(recaps_4X_pair1, recaps_4X_pair2, recaps_4X_pair3)
# pair_5X = chooseRecapPair_5X(recaps_5X_pair1, recaps_5X_pair2, recaps_5X_pair3, recaps_5X_pair4)
# pair_6X = chooseRecapPair_6X(recaps_6X_pair1, recaps_6X_pair2, recaps_6X_pair3, recaps_6X_pair4, recaps_6X_pair5)
# pair_7X = chooseRecapPair_7X(recaps_7X_pair1, recaps_7X_pair2, recaps_7X_pair3, recaps_7X_pair4, recaps_7X_pair5, recaps_7X_pair6)
# pair_8X = chooseRecapPair_8X(recaps_8X_pair1, recaps_8X_pair2, recaps_8X_pair3, recaps_8X_pair4, recaps_8X_pair5, recaps_8X_pair6, recaps_8X_pair7)
# 
# # Join pairs together into one data frame where each fish has only one pair (even if caught multiple times) - 464 fish (was 507 fish earlier)
# recap_pairs_1perfish <- rbind(recaps_2X_pair1, pair_3X, pair_4X, pair_5X, pair_6X, pair_7X, pair_8X)
# 
# # Pull out any pairs where the fish was caught twice on the same day - cuts pairs down from 772 to 761 (pre-addition of 2016-2018 data, cut pairs down from 545 to 534)
# recap_pairs <- recap_pairs %>%
#   filter(tal != 0)
# 
# recap_pairs_1perfish <- recap_pairs_1perfish %>%
#   filter(tal != 0)
# 
# # Add in percentage-of-year column for time 
# recap_pairs <- recap_pairs %>%
#   dplyr::rename(tal_days = tal) %>%
#   mutate(tal = as.numeric(tal_days/365))
# 
# recap_pairs_1perfish <- recap_pairs_1perfish %>%
#   dplyr::rename(tal_days = tal) %>%
#   mutate(tal = as.numeric(tal_days/365))
# 
# # Add in column for growth, growth per year, make size columns numeric, add in projected growth 1 year later
# recap_pairs <- recap_pairs %>%
#   mutate(L1 = as.numeric(L1),
#          L2 = as.numeric(L2),
#          growth = L2-L1,
#          growth_per_year = growth/tal,
#          L2_one_year_later = round(L1 + growth_per_year, digits = 1))  # projected growth one year later, based on growth-per-year
# 
# recap_pairs_1perfish <- recap_pairs_1perfish %>%
#   mutate(L1 = as.numeric(L1),
#          L2 = as.numeric(L2),
#          growth = L2-L1,
#          growth_per_year = growth/tal,
#          L2_one_year_later = round(L1 + growth_per_year, digits = 1))
# 
# # ##### Pull out fish captured about a year later first, then choose among the pairs
# # recap_pairs_year <- recap_pairs %>%
# #   filter(tal_days > year_lower_limit & tal_days <= year_upper_limit) %>%
# #   group_by(fish_id) %>%
# #   summarize(ncaps_about1year = n(),
# #             year1_about1year = min(year)) %>%
# #   filter(ncaps_about1year > 1)



# # Just use recaptures more than one month apart
# recap_pairs_1monthplus <- recap_pairs %>%  # cuts pairs down from xx to 540 (old: 761 to 618)
#   filter(tal_days > 30)
# 
# recap_pairs_1perfish_1monthplus <- recap_pairs_1perfish %>% # cuts down to 374 (old: 412 pairs (from 502))
#   filter(tal_days > 30)
# 
# # Pull out a set where fish were captured about a year apart
# year_lower_limit = 345
# year_upper_limit = 385
# 
# recap_pairs_yearrecap <- recap_pairs %>% filter(tal_days > year_lower_limit & tal_days <= year_upper_limit)  # 244 pairs

# # Check out outlier growth rates... maybe think about some extreme outliers from the analysis and seeing if/how results change?
# recap_pairs %>% filter(growth_per_year < -10)
# recap_pairs %>% filter(growth_per_year > 10)
# recap_pairs_1perfish_1monthplus %>% filter(growth_per_year <= -10 | growth_per_year > 10)  # none 

# 
# Quite generally you want the vcov function which provides the complete parameter covariance matrix. To get the regular asymptotic standard errors reported by summary you can use
# 
# se <- sqrt(diag(vcov(model)))
# 
# # Testing standard error:
# k_test = length(model_1$coefficients) - 1
# SSE = sum(model_1$residuals**2)
# 
# # Try including random effects
# # filter out for at least 2 sets of recaps - only 50 pairs
# recaps_model_1a <- recap_pairs_yearrecap %>% 
#   group_by(fish_id) %>%
#   mutate(n_pairs = n()) %>%
#   filter(n_pairs >= 2) %>%
#   ungroup()
# 
# model_1a <- lmer(L2 ~ L1 + (L1 | fish_id), data = recap_pairs_yearrecap)
# ggplot(data = model_1_runs) +
#   geom_abline(aes(intercept = intercept_est[1], slope = slope_est[1])) 
# 
# 
# 
# geom_abline(aes(
#   intercept = b0, 
#   slope = b1,
#   linetype = if_else(metric == "pred", "", "dashed")),
# ) +
#   
#   # Plot the estimated gam for fish caught at least a month apart, only one recapture pair per fish
#   
#   
# # COULD DEFINITELY DO THIS IN A MORE STREAMLINED WAY!!
# # Function to choose one of the recap pairs for a fish caught 3 times 
# chooseRecapPair_3X <- function(pair1, pair2) {
#   # Set up output data frame
#   out <- data.frame(fish_id = pair1$fish_id)
#   out <- out %>%
#     mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
#   
#   for(i in 1:length(out$fish_id)) {
#     pair <- base::sample(1:2, 1)
#     
#     if(pair == 1) {
#       pair_choice <- pair1
#     } else if (pair == 2) {
#       pair_choice <- pair2
#     }
#     
#     out$L1[i] <- pair_choice$L1[i]
#     out$L2[i] <- pair_choice$L2[i]
#     out$tal[i] <- pair_choice$tal[i]
#   }
#   return(out)
# }
# 
# chooseRecapPair_4X <- function(pair1, pair2, pair3) {
#   # Set up output data frame
#   out <- data.frame(fish_id = pair1$fish_id)
#   out <- out %>%
#     mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
#   
#   for(i in 1:length(out$fish_id)) {
#     pair <- base::sample(1:3, 1)
#     
#     if(pair == 1) {
#       pair_choice <- pair1
#     } else if (pair == 2) {
#       pair_choice <- pair2
#     } else if (pair == 3) {
#       pair_choice <- pair3
#     }
#     
#     out$L1[i] <- pair_choice$L1[i]
#     out$L2[i] <- pair_choice$L2[i]
#     out$tal[i] <- pair_choice$tal[i]
#   }
#   return(out)
# }
# 
# chooseRecapPair_5X <- function(pair1, pair2, pair3, pair4) {
#   # Set up output data frame
#   out <- data.frame(fish_id = pair1$fish_id)
#   out <- out %>%
#     mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
#   
#   for(i in 1:length(out$fish_id)) {
#     pair <- base::sample(1:4, 1)
#     
#     if(pair == 1) {
#       pair_choice <- pair1
#     } else if (pair == 2) {
#       pair_choice <- pair2
#     } else if (pair == 3) {
#       pair_choice <- pair3
#     } else if (pair == 4) {
#       pair_choice <- pair4
#     }
#     
#     out$L1[i] <- pair_choice$L1[i]
#     out$L2[i] <- pair_choice$L2[i]
#     out$tal[i] <- pair_choice$tal[i]
#   }
#   return(out)
# }
# 
# chooseRecapPair_6X <- function(pair1, pair2, pair3, pair4, pair5) {
#   # Set up output data frame
#   out <- data.frame(fish_id = pair1$fish_id)
#   out <- out %>%
#     mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
#   
#   for(i in 1:length(out$fish_id)) {
#     pair <- base::sample(1:5, 1)
#     
#     if(pair == 1) {
#       pair_choice <- pair1
#     } else if (pair == 2) {
#       pair_choice <- pair2
#     } else if (pair == 3) {
#       pair_choice <- pair3
#     } else if (pair == 4) {
#       pair_choice <- pair4
#     } else if (pair == 5) {
#       pair_choice <- pair5
#     }
#     
#     out$L1[i] <- pair_choice$L1[i]
#     out$L2[i] <- pair_choice$L2[i]
#     out$tal[i] <- pair_choice$tal[i]
#   }
#   return(out)
# }
# 
# chooseRecapPair_7X <- function(pair1, pair2, pair3, pair4, pair5, pair6) {
#   # Set up output data frame
#   out <- data.frame(fish_id = pair1$fish_id)
#   out <- out %>%
#     mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
#   
#   for(i in 1:length(out$fish_id)) {
#     pair <- base::sample(1:6, 1)
#     
#     if(pair == 1) {
#       pair_choice <- pair1
#     } else if (pair == 2) {
#       pair_choice <- pair2
#     } else if (pair == 3) {
#       pair_choice <- pair3
#     } else if (pair == 4) {
#       pair_choice <- pair4
#     } else if (pair == 5) {
#       pair_choice <- pair5
#     } else if (pair == 6) {
#       pair_choice <- pair6
#     }
#     
#     out$L1[i] <- pair_choice$L1[i]
#     out$L2[i] <- pair_choice$L2[i]
#     out$tal[i] <- pair_choice$tal[i]
#   }
#   return(out)
# }
# 
# chooseRecapPair_8X <- function(pair1, pair2, pair3, pair4, pair5, pair6, pair7) {
#   # Set up output data frame
#   out <- data.frame(fish_id = pair1$fish_id)
#   out <- out %>%
#     mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
#   
#   for(i in 1:length(out$fish_id)) {
#     pair <- base::sample(1:7, 1)
#     
#     if(pair == 1) {
#       pair_choice <- pair1
#     } else if (pair == 2) {
#       pair_choice <- pair2
#     } else if (pair == 3) {
#       pair_choice <- pair3
#     } else if (pair == 4) {
#       pair_choice <- pair4
#     } else if (pair == 5) {
#       pair_choice <- pair5
#     } else if (pair == 6) {
#       pair_choice <- pair6
#     } else if (pair == 7) {
#       pair_choice <- pair7
#     }
#     
#     out$L1[i] <- pair_choice$L1[i]
#     out$L2[i] <- pair_choice$L2[i]
#     out$tal[i] <- pair_choice$tal[i]
#   }
#   return(out)
# }
# 
# # Choose one recap pair per fish, run lm, save coefficients
# chooseRecapFitModel <- function(n_runs, recap_pairs_list, time_lower_limit, time_upper_limit) {
#   
#   # Set up an output data frame
#   out = data.frame(run = seq(from=1, to=n_runs, by=1),
#                    intercept_est = NA,
#                    slope_est = NA,
#                    intercept_se = NA,
#                    slope_se = NA)
#   
#   for(i in 1:n_runs) {
#     # Choose a pair for fish caught mulitple times
#     pair_3X_func = chooseRecapPair_3X(recap_pairs_list$recaps_3X_pair1, recap_pairs_list$recaps_3X_pair2)
#     pair_4X_func = chooseRecapPair_4X(recap_pairs_list$recaps_4X_pair1, recap_pairs_list$recaps_4X_pair2, recap_pairs_list$recaps_4X_pair3)
#     pair_5X_func = chooseRecapPair_5X(recap_pairs_list$recaps_5X_pair1, recap_pairs_list$recaps_5X_pair2, recap_pairs_list$recaps_5X_pair3, recap_pairs_list$recaps_5X_pair4)
#     pair_6X_func = chooseRecapPair_6X(recap_pairs_list$recaps_6X_pair1, recap_pairs_list$recaps_6X_pair2, recap_pairs_list$recaps_6X_pair3, recap_pairs_list$recaps_6X_pair4, recap_pairs_list$recaps_6X_pair5)
#     pair_7X_func = chooseRecapPair_7X(recap_pairs_list$recaps_7X_pair1, recap_pairs_list$recaps_7X_pair2, recap_pairs_list$recaps_7X_pair3, recap_pairs_list$recaps_7X_pair4, recap_pairs_list$recaps_7X_pair5, recap_pairs_list$recaps_7X_pair6)
#     pair_8X_func = chooseRecapPair_8X(recap_pairs_list$recaps_8X_pair1, recap_pairs_list$recaps_8X_pair2, recap_pairs_list$recaps_8X_pair3, recap_pairs_list$recaps_8X_pair4, 
#                                       recap_pairs_list$recaps_8X_pair5, recap_pairs_list$recaps_8X_pair6, recap_pairs_list$recaps_8X_pair7)
#     
#     # Bind together
#     recap_pairs_func <- rbind(recap_pairs_list$recaps_2X_pair1, pair_3X_func, pair_4X_func, pair_5X_func, pair_6X_func, pair_7X_func, pair_8X_func)
#     
#     # Filter and process 
#     recap_pairs_func <- recap_pairs_func %>%
#       #filter(tal != 0) %>%  # filter out any recaps caught on the same day
#       filter(tal < time_lower_limit & tal >= time_upper_limit) %>%  # only use ones recaught within about a year
#       dplyr::rename(tal_days= tal) %>%
#       mutate(tal = as.numeric(tal_days/365)) %>%  # add percentage-of-year column
#       mutate(L1 = as.numeric(L1),
#              L2 = as.numeric(L2),
#              growth = L2 - L1,
#              growth_per_year = growth/tal)
#     
#     # Fit model 
#     model_out <- lm(L2 ~ L1, data = recap_pairs_func)
#     
#     out$intercept_est[i] = coef(summary(model_out))["(Intercept)", "Estimate"]
#     out$slope_est[i] = coef(summary(model_out))["L1", "Estimate"]
#     out$intercept_se[i] = coef(summary(model_out))["(Intercept)", "Std. Error"]
#     out$slope_se[i] = coef(summary(model_out))["L1", "Std. Error"]
#   }
#   
#   return(out)
# }
# # List of pairs, for feeding into function
# recap_pairs_list <- list(recaps_2X_pair1, recaps_3X_pair1 = recaps_3X_pair1, recaps_3X_pair2 = recaps_3X_pair2,
#                          recaps_4X_pair1 = recaps_4X_pair1, recaps_4X_pair2 = recaps_4X_pair2, recaps_4X_pair3 = recaps_4X_pair3,
#                          recaps_5X_pair1 = recaps_5X_pair1, recaps_5X_pair2 = recaps_5X_pair2, recaps_5X_pair3 = recaps_5X_pair3, recaps_5X_pair4 = recaps_5X_pair4,
#                          recaps_6X_pair1 = recaps_6X_pair1, recaps_6X_pair2 = recaps_6X_pair2, recaps_6X_pair3 = recaps_6X_pair3,
#                          recaps_6X_pair4 = recaps_6X_pair4, recaps_6X_pair5 = recaps_6X_pair5,
#                          recaps_7X_pair1 = recaps_7X_pair1, recaps_7X_pair2 = recaps_7X_pair2, recaps_7X_pair3 = recaps_7X_pair3,
#                          recaps_7X_pair4 = recaps_7X_pair4, recaps_7X_pair5 = recaps_7X_pair5, recaps_7X_pair6 = recaps_7X_pair6,
#                          recaps_8X_pair1 = recaps_8X_pair1, recaps_8X_pair2 = recaps_8X_pair2, recaps_8X_pair3 = recaps_8X_pair3,
#                          recaps_8X_pair4 = recaps_8X_pair4, recaps_8X_pair5 = recaps_8X_pair5, recaps_8X_pair6 = recaps_8X_pair6,
#                          recaps_8X_pair7 = recaps_8X_pair7)
# 
# 
# 
# plot(recap_pairs_1perfish_yearrecap$L1, recap_pairs_1perfish_yearrecap$L2, xlab = "length in year 1", ylab = "length in year 2")
# points(size_range_sequence$L1, model_1_predictions$fit, type = "l", col = "black", lty = 2, lwd = 3)
# #points(size_range_sequence$L1, model_1_predictions$fit + model_1_predictions$se.fit, type = "l", col = "blue", lty = 2, lwd = 3)
# #points(size_range_sequence$L1, model_1_predictions$fit - model_1_predictions$se.fit, type = "l", col = "blue", lty = 2, lwd = 3)
# 
# #points(size_range_sequence$L1, test_grow_predictions_2, type = "l", col = "red", lty = 2, lwd = 3)
# 
# ##### Run with all recaps (1 per fish), using growth per year to estimate length after a year
# model_2 <- lm(L2_one_year_later ~ L1, data = recap_pairs_1perfish)
# model_2_k <- as.numeric(-log(model_2$coefficients[2])) # (fishmethods lowest AIC estimate of 0.867)
# model_2_Linf <- as.numeric(model_2$coefficients[1]/(1 - model_2$coefficients[2]))  # (fishmethods lowest AIC estimate of 10.64)
# 

# ##### Rees et al 2014 gam method
# 
# # Fitting the data with a spline, from Rees et al. 2014 IPM paper, code modified from their supplementary script gamExample.R
# # Note from Rees script: Note: m=3 specifies 3rd derivative penalty so the fit can have
# # nonzero curvature at the endpoints of the data range.
# # With the default (m=2), curvature -> 0 at the endpoints.
# 
# # Just pull out fish recaptured within about a year 
# recap_pairs_1perfish_yearrecap <- recap_pairs_1perfish %>% filter(tal_days > 345 & tal_days <= 385)  # 180 fish
# test_grow <- gam(L2~s(L1, m=3), data = recap_pairs_1perfish_yearrecap)
# test_grow_2 <- gam(L2~s(L1, m=2), data = recap_pairs_1perfish_yearrecap)
# 
# # Predict new data
# size_range_sequence <- data.frame(L1 = seq(1, max_size, length = max_size*10))
# test_grow_predictions <- predict(test_grow, newdata = size_range_sequence, type = "response")
# test_grow_predictions_2 <- predict(test_grow_2, newdata = size_range_sequence, type = "response")
# 
# sse_test_grow_1 <- sum(test_grow$residuals^2)
# sd_predictions_1 <- sqrt(sse_test_grow_1/test_grow$df.residual)	
# sse_test_grow_2 <- sum(test_grow_2$residuals^2)
# sd_predictions_2 <- sqrt(sse_test_grow_2/test_grow_2$df.residual)	
# 
# #X <- data.frame(z=z,z1=z1)
# #gamGrow <- gam(z1~s(z,m=3),data=X);
# 
# plot(recap_pairs_1perfish_yearrecap$L1, recap_pairs_1perfish_yearrecap$L2, xlab = "length in year 1", ylab = "length in year 2")
# points(size_range_sequence$L1, test_grow_predictions, type = "l", col = "black", lty = 2, lwd = 3)
# points(size_range_sequence$L1, test_grow_predictions_2, type = "l", col = "red", lty = 2, lwd = 3)


