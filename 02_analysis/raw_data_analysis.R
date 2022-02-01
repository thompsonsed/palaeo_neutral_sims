## Analysis of the raw data for plotting some informative graphs

library(tidyverse)

tetrapod_df <- read.csv(file.path("input", "early_tetrapod_occurrences_10_04_19.csv"))
interval_data <- read.csv(file.path("input", "interval_data.csv")) %>% select(-X)
tetrapod_df %>% group_by(early_interval) %>%
	summarise(species_richness = sum(n_occs),
		  num_collections=length(unique(collection_no))) %>%
	rename(interval=early_interval) %>%
	inner_join(interval_data) %>%
	gather(key="metric", value="value", species_richnes, num_collections, )
	ggplot() +
	geom_line(aes(x=midpoint, y=species_richness))
