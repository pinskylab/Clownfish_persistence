# Find data characteristics

#################### Set-up: ####################

# Pull data, functions, constants 
source(here::here('Code', 'Constants_database_common_functions.R'))

#################### Running things: ####################
size_by_color_metrics <- allfish_caught %>%
  filter(color %in% c('YP','O','YR')) %>%  # just filter the main three for now (female, male, juvenile)
  filter(!is.na(size)) %>%  # take out any fish with NA size
  filter(dive_type == 'C') %>%  # only clownfish dives, to prevent recaptures and such
  group_by(color) %>%
  summarize(mean = mean(size),
            sd = sd(size),
            min = min(size),
            max = max(size))

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

#################### Saving output: ####################
save(size_by_color_metrics, file=here::here('Data','size_by_color_metrics.RData'))
