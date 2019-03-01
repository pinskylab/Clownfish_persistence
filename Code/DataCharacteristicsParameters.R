# Find data characteristics

#################### Set-up: ####################

# Pull data, functions, constants 
source(here::here('Code', 'Constants_database_common_functions.R'))

#################### Running things: ####################
# Summarize size distributions by color
size_by_color_metrics <- allfish_caught %>%
  filter(color %in% c('YP','O','YR')) %>%  # just filter the main three for now (female, male, juvenile)
  filter(!is.na(size)) %>%  # take out any fish with NA size
  filter(dive_type == 'C') %>%  # only clownfish dives, to prevent recaptures and such
  group_by(color) %>%
  summarize(mean = mean(size),
            sd = sd(size),
            min = min(size),
            max = max(size))

# Pull out all sizes for breeding females
female_sizes <- allfish_caught %>%
  filter(color == 'YP' | (color == 'Y' & year == 2012 & size >= 6.0)) %>%
  filter(dive_type %in% c("0","C","D","E")) %>% # check that these are the dive types to use
  filter(!is.na(size))
  
#################### Plots: ####################
# Histograms of size by tail color - each color separate, only plotting YP/O/YR fish from C dives
pdf(file = here::here('Plots/DataCharacteristics', 'Tail_color_size_distributions.pdf'))
ggplot(data = (allfish_caught %>% filter(color %in% c('YP', 'YR', 'O'), !is.na(size), dive_type == 'C')), aes(x=size, fill=color)) +
  geom_histogram(binwidth = 0.5) +
  facet_wrap(~color) +
  theme_bw() +
  ggtitle('Size distributions by color')
dev.off()

# Histogram of size by tail color - all colors together, only plotting YP/O/YR fish from C dives
pdf(file = here::here('Plots/DataCharacteristics', 'Tail_color_size_distributions_all_together.pdf'))
ggplot(data = (allfish_caught %>% filter(color %in% c('YP', 'YR', 'O'), !is.na(size), dive_type == 'C')), aes(x=size, fill=color)) +
  geom_histogram(binwidth = 0.5) +
  theme_bw() +
  ggtitle('Size distributions by color')
dev.off()

# Histogram of just female sizes, dives from C, O, E, D, includes Y over 6cm in 2012
pdf(file = here::here('Plots/DataCharacteristics', 'Female_size_distributions.pdf'))
ggplot(data = female_sizes, aes(x=size)) +
  geom_histogram(binwidth = 0.5) +
  theme_bw() +
  ggtitle('Size distributions of females') 
dev.off()


#################### Saving output: ####################
save(size_by_color_metrics, file=here::here('Data','size_by_color_metrics.RData'))
save(female_sizes, file=here::here('Data','female_sizes.RData'))
