# Fitting a growth curve (pulling code from growth_newdata.R in Growth repository, figure out how to credit Michelle correctly)

#################### Set-up: ####################
source(here::here("Code", "Constants_database_common_functions.R"))

# These are already loaded in Constants_database_common_functions.R
# library(dplyr)
# library(tidyr)
# library(lubridate)
# library(stringr)

# library(fishmethods)
# library(readr)

library(ggplot2)
library(grid)
library(gridExtra)
library(mgcv)  # used by Rees et al. 2014
library(lme4)  # used by Hart et al. 2009

##### Load files from other scripts within this repository or source those scripts (below, commented out)

##### Set up parameters

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

# COULD DEFINITELY DO THIS IN A MORE STREAMLINED WAY!!
# Function to choose one of the recap pairs for a fish caught 3 times 
chooseRecapPair_3X <- function(pair1, pair2) {
  # Set up output data frame
  out <- data.frame(fish_id = pair1$fish_id)
  out <- out %>%
    mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
  
  for(i in 1:length(out$fish_id)) {
    pair <- base::sample(1:2, 1)
    
    if(pair == 1) {
      pair_choice <- pair1
    } else if (pair == 2) {
      pair_choice <- pair2
    }
    
    out$L1[i] <- pair_choice$L1[i]
    out$L2[i] <- pair_choice$L2[i]
    out$tal[i] <- pair_choice$tal[i]
  }
  return(out)
}

chooseRecapPair_4X <- function(pair1, pair2, pair3) {
  # Set up output data frame
  out <- data.frame(fish_id = pair1$fish_id)
  out <- out %>%
    mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
  
  for(i in 1:length(out$fish_id)) {
    pair <- base::sample(1:3, 1)
    
    if(pair == 1) {
      pair_choice <- pair1
    } else if (pair == 2) {
      pair_choice <- pair2
    } else if (pair == 3) {
      pair_choice <- pair3
    }
    
    out$L1[i] <- pair_choice$L1[i]
    out$L2[i] <- pair_choice$L2[i]
    out$tal[i] <- pair_choice$tal[i]
  }
  return(out)
}

chooseRecapPair_5X <- function(pair1, pair2, pair3, pair4) {
  # Set up output data frame
  out <- data.frame(fish_id = pair1$fish_id)
  out <- out %>%
    mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
  
  for(i in 1:length(out$fish_id)) {
    pair <- base::sample(1:4, 1)
    
    if(pair == 1) {
      pair_choice <- pair1
    } else if (pair == 2) {
      pair_choice <- pair2
    } else if (pair == 3) {
      pair_choice <- pair3
    } else if (pair == 4) {
      pair_choice <- pair4
    }
    
    out$L1[i] <- pair_choice$L1[i]
    out$L2[i] <- pair_choice$L2[i]
    out$tal[i] <- pair_choice$tal[i]
  }
  return(out)
}

chooseRecapPair_6X <- function(pair1, pair2, pair3, pair4, pair5) {
  # Set up output data frame
  out <- data.frame(fish_id = pair1$fish_id)
  out <- out %>%
    mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
  
  for(i in 1:length(out$fish_id)) {
    pair <- base::sample(1:5, 1)
    
    if(pair == 1) {
      pair_choice <- pair1
    } else if (pair == 2) {
      pair_choice <- pair2
    } else if (pair == 3) {
      pair_choice <- pair3
    } else if (pair == 4) {
      pair_choice <- pair4
    } else if (pair == 5) {
      pair_choice <- pair5
    }
    
    out$L1[i] <- pair_choice$L1[i]
    out$L2[i] <- pair_choice$L2[i]
    out$tal[i] <- pair_choice$tal[i]
  }
  return(out)
}

chooseRecapPair_7X <- function(pair1, pair2, pair3, pair4, pair5, pair6) {
  # Set up output data frame
  out <- data.frame(fish_id = pair1$fish_id)
  out <- out %>%
    mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
  
  for(i in 1:length(out$fish_id)) {
    pair <- base::sample(1:6, 1)
    
    if(pair == 1) {
      pair_choice <- pair1
    } else if (pair == 2) {
      pair_choice <- pair2
    } else if (pair == 3) {
      pair_choice <- pair3
    } else if (pair == 4) {
      pair_choice <- pair4
    } else if (pair == 5) {
      pair_choice <- pair5
    } else if (pair == 6) {
      pair_choice <- pair6
    }
    
    out$L1[i] <- pair_choice$L1[i]
    out$L2[i] <- pair_choice$L2[i]
    out$tal[i] <- pair_choice$tal[i]
  }
  return(out)
}

chooseRecapPair_8X <- function(pair1, pair2, pair3, pair4, pair5, pair6, pair7) {
  # Set up output data frame
  out <- data.frame(fish_id = pair1$fish_id)
  out <- out %>%
    mutate(L1 = rep(NA), L2 = rep(NA), tal = rep(NA))
  
  for(i in 1:length(out$fish_id)) {
    pair <- base::sample(1:7, 1)
    
    if(pair == 1) {
      pair_choice <- pair1
    } else if (pair == 2) {
      pair_choice <- pair2
    } else if (pair == 3) {
      pair_choice <- pair3
    } else if (pair == 4) {
      pair_choice <- pair4
    } else if (pair == 5) {
      pair_choice <- pair5
    } else if (pair == 6) {
      pair_choice <- pair6
    } else if (pair == 7) {
      pair_choice <- pair7
    }
    
    out$L1[i] <- pair_choice$L1[i]
    out$L2[i] <- pair_choice$L2[i]
    out$tal[i] <- pair_choice$tal[i]
  }
  return(out)
}

#################### Running things: ####################

########## Finding recaptures and prepping the data ##########

##### Assign each individual fish marked and recaptured a fish_id (since one fish can have multiple tag_ids and a gen_id)
# Pull out all fish that have been marked in some way, do a first-pass at assigning a fish_id (this and below is from the updated way of connecting fish, from the growth_newdata.R script in Growth repo)
marked_fish <- allfish_caught %>%
  dplyr::select(sample_id, color, size, sex, tag_id, gen_id, fish_table_id, anem_table_id, dive_table_id, date, site, year, fish_notes) %>%
  filter(!is.na(gen_id) | (tag_id != 'NA' & !is.na(tag_id))) %>% # pull out fish "tagged" in any way, either PIT or via genetic sample
  mutate(fish_id = case_when(!is.na(gen_id) ~ paste('gen', gen_id, sep=''),  # if fish has a gen_id, use that as the fish_id
                             is.na(gen_id) ~ paste('tag', tag_id, sep='')))  # if it doesn't have a gen_id, use the tag_id

# Go through and match up fish ided in mulitple ways (now that gen_id is more stable than tag_id - because fish lose their tags - use that as the primary)
# one tagged fish has two gen_ids - see GitHub issue, can directly assign id to that fish
for(i in 1:length(marked_fish$fish_id)) {
  if(substring(marked_fish$fish_id[i],1,3) == 'tag') {  # if it has a fish_id based on tag_id rather than gen_id...
    tag_id_val <- marked_fish$tag_id[i]  # pull out the tag_id
    
    matches <- marked_fish %>%
      filter(tag_id == tag_id_val)  # filter out any other cases that match that tag_id
    
    for(j in 1:length(matches$fish_id)) {
      if(!is.na(matches$gen_id[j])) {  # if any of the captures with that tag_id have a gen_id
        marked_fish$fish_id[i] <- matches$fish_id[j]  # update the fish_id to the fish_id based on the match with a gen_id - WHAT IF HAS MULTIPLE GEN_IDs?
      }
    }
  }
}

# Remove any observations where size is not recorded - goes from XX (3321 distinct fish) to 4012 obs  (old info: goes from 4001 observation to 3994 observations (pre-2016-2018 genetic data numbers: goes from 3596 observations to 3592)
marked_fish <- marked_fish %>%
  filter(size != "NA" & !is.na(size))

# Remove any observations where the date is not recorded (doesn't remove any)
marked_fish <- marked_fish %>% filter(!is.na(date))

##### Pull out recaptured fish, create recapture pairs
# Pull out fish captured more than once - 464 obs
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

# Choose just one pair for fish caught multiple times
pair_3X = chooseRecapPair_3X(recaps_3X_pair1, recaps_3X_pair2)
pair_4X = chooseRecapPair_4X(recaps_4X_pair1, recaps_4X_pair2, recaps_4X_pair3)
pair_5X = chooseRecapPair_5X(recaps_5X_pair1, recaps_5X_pair2, recaps_5X_pair3, recaps_5X_pair4)
pair_6X = chooseRecapPair_6X(recaps_6X_pair1, recaps_6X_pair2, recaps_6X_pair3, recaps_6X_pair4, recaps_6X_pair5)
pair_7X = chooseRecapPair_7X(recaps_7X_pair1, recaps_7X_pair2, recaps_7X_pair3, recaps_7X_pair4, recaps_7X_pair5, recaps_7X_pair6)
pair_8X = chooseRecapPair_8X(recaps_8X_pair1, recaps_8X_pair2, recaps_8X_pair3, recaps_8X_pair4, recaps_8X_pair5, recaps_8X_pair6, recaps_8X_pair7)

# Join pairs together into one data frame - 697 pairs (earlier was 772 pairs? - maybe PIT typos were fixed?) (was 545 pairs before addition of 2016-2018 genetic data) 
recap_pairs <- rbind(recaps_2X_pair1, recaps_3X_pair1, recaps_3X_pair2,
                     recaps_4X_pair1, recaps_4X_pair2, recaps_4X_pair3,
                     recaps_5X_pair1, recaps_5X_pair2, recaps_5X_pair3, recaps_5X_pair4,
                     recaps_6X_pair1, recaps_6X_pair2, recaps_6X_pair3, recaps_6X_pair4, recaps_6X_pair5,
                     recaps_7X_pair1, recaps_7X_pair2, recaps_7X_pair3, 
                     recaps_7X_pair4, recaps_7X_pair5, recaps_7X_pair6,
                     recaps_8X_pair1, recaps_8X_pair2, recaps_8X_pair3,
                     recaps_8X_pair4, recaps_8X_pair5, recaps_8X_pair6, recaps_8X_pair7)

# Join pairs together into one data frame where each fish has only one pair (even if caught multiple times) - 464 fish (was 507 fish earlier)
recap_pairs_1perfish <- rbind(recaps_2X_pair1, pair_3X, pair_4X, pair_5X, pair_6X, pair_7X, pair_8X)

# Pull out any pairs where the fish was caught twice on the same day - cuts pairs down from 772 to 761 (pre-addition of 2016-2018 data, cut pairs down from 545 to 534)
recap_pairs <- recap_pairs %>%
  filter(tal != 0)

recap_pairs_1perfish <- recap_pairs_1perfish %>%
  filter(tal != 0)

# Add in percentage-of-year column for time 
recap_pairs <- recap_pairs %>%
  dplyr::rename(tal_days = tal) %>%
  mutate(tal = as.numeric(tal_days/365))

recap_pairs_1perfish <- recap_pairs_1perfish %>%
  dplyr::rename(tal_days = tal) %>%
  mutate(tal = as.numeric(tal_days/365))

# Add in column for growth, growth per year, make size columns numeric
recap_pairs <- recap_pairs %>%
  mutate(L1 = as.numeric(L1),
         L2 = as.numeric(L2),
         growth = L2-L1,
         growth_per_year = growth/tal)

recap_pairs_1perfish <- recap_pairs_1perfish %>%
  mutate(L1 = as.numeric(L1),
         L2 = as.numeric(L2),
         growth = L2-L1,
         growth_per_year = growth/tal)

# Just use recaptures more than one month apart
recap_pairs_1monthplus <- recap_pairs %>%  # cuts pairs down from xx to 540 (old: 761 to 618)
  filter(tal_days > 30)

recap_pairs_1perfish_1monthplus <- recap_pairs_1perfish %>% # cuts down to 374 (old: 412 pairs (from 502))
  filter(tal_days > 30)

# Pull out a set where fish were captured about a year apart
year_lower_limit = 345
year_upper_limit = 385

recap_pairs_yearrecap <- recap_pairs %>% filter(tal_days > year_lower_limit & tal_days <= year_upper_limit)  # 244 pairs

# # Check out outlier growth rates... maybe think about some extreme outliers from the analysis and seeing if/how results change?
# recap_pairs %>% filter(growth_per_year < -10)
# recap_pairs %>% filter(growth_per_year > 10)
# recap_pairs_1perfish_1monthplus %>% filter(growth_per_year <= -10 | growth_per_year > 10)  # none 

########## Fitting growth models ##########

##### Hart et al. 2009 method
# If just do a linear regression between Lt and Lt+1, can calculate K and Linf from y-intercept (b) and slope (m) - but doesn't take into account random effect or individual effects
model_1 <- lm(L2 ~ L1, data =  recap_pairs_1perfish_yearrecap)

model_1_k <- as.numeric(-log(model_1$coefficients[2]))  # 0.881, pretty similar to fishmethods lowest AIC estimate of 0.867
model_1_Linf <- as.numeric(model_1$coefficients[1]/(1 - model_1$coefficients[2]))  # 10.69, pretty similar to fishmethods lowest AIC estimate of 10.64

size_range_sequence <- data.frame(L1 = seq(1, max_size, length = max_size*10))

model_1_predictions <- predict.lm(model_1, newdata = size_range_sequence, se.fit = TRUE)

# Try including random effects
# filter out for at least 2 sets of recaps - only 50 pairs
recaps_model_1a <- recap_pairs_yearrecap %>% 
  group_by(fish_id) %>%
  mutate(n_pairs = n()) %>%
  filter(n_pairs >= 2) %>%
  ungroup()

model_1a <- lmer(L2 ~ L1 + (L1 | fish_id), data = recap_pairs_yearrecap)

plot(recap_pairs_1perfish_yearrecap$L1, recap_pairs_1perfish_yearrecap$L2, xlab = "length in year 1", ylab = "length in year 2")
points(size_range_sequence$L1, model_1_predictions$fit, type = "l", col = "black", lty = 2, lwd = 3)
#points(size_range_sequence$L1, model_1_predictions$fit + model_1_predictions$se.fit, type = "l", col = "blue", lty = 2, lwd = 3)
#points(size_range_sequence$L1, model_1_predictions$fit - model_1_predictions$se.fit, type = "l", col = "blue", lty = 2, lwd = 3)

#points(size_range_sequence$L1, test_grow_predictions_2, type = "l", col = "red", lty = 2, lwd = 3)


##### Rees et al 2014 gam method

# Fitting the data with a spline, from Rees et al. 2014 IPM paper, code modified from their supplementary script gamExample.R
# Note from Rees script: Note: m=3 specifies 3rd derivative penalty so the fit can have
# nonzero curvature at the endpoints of the data range.
# With the default (m=2), curvature -> 0 at the endpoints.

# Just pull out fish recaptured within about a year 
recap_pairs_1perfish_yearrecap <- recap_pairs_1perfish %>% filter(tal_days > 345 & tal_days <= 385)  # 180 fish
test_grow <- gam(L2~s(L1, m=3), data = recap_pairs_1perfish_yearrecap)
test_grow_2 <- gam(L2~s(L1, m=2), data = recap_pairs_1perfish_yearrecap)

# Predict new data
size_range_sequence <- data.frame(L1 = seq(1, max_size, length = max_size*10))
test_grow_predictions <- predict(test_grow, newdata = size_range_sequence, type = "response")
test_grow_predictions_2 <- predict(test_grow_2, newdata = size_range_sequence, type = "response")

sse_test_grow_1 <- sum(test_grow$residuals^2)
sd_predictions_1 <- sqrt(sse_test_grow_1/test_grow$df.residual)	
sse_test_grow_2 <- sum(test_grow_2$residuals^2)
sd_predictions_2 <- sqrt(sse_test_grow_2/test_grow_2$df.residual)	

#X <- data.frame(z=z,z1=z1)
#gamGrow <- gam(z1~s(z,m=3),data=X);

plot(recap_pairs_1perfish_yearrecap$L1, recap_pairs_1perfish_yearrecap$L2, xlab = "length in year 1", ylab = "length in year 2")
points(size_range_sequence$L1, test_grow_predictions, type = "l", col = "black", lty = 2, lwd = 3)
points(size_range_sequence$L1, test_grow_predictions_2, type = "l", col = "red", lty = 2, lwd = 3)


#################### Plots: ####################

# Plot the estimated gam for fish caught at least a month apart, only one recapture pair per fish


#################### Saving output: ####################