# Determine the best-fitting simulations

library(tidyverse)

results_dir_fragmented <- "results/Sim8"
parameters_df <- read.csv(file.path("input", "paleo_parameters.csv")) %>%
	mutate(interval = tools::toTitleCase(as.character(interval)),
	       tetrapod_group = tools::toTitleCase(as.character(tetrapod_group))) %>% select(-X) %>%
	mutate(deme = floor(deme))
best_sims <- read.csv(file.path(results_dir_fragmented, "single_param_best_5.csv")) %>%
	mutate(percent_cover = ifelse(scenario == "fragmented", percent_cover, 100))%>%
	inner_join(parameters_df) %>%
	select(speciation_rate, scenario, sigma, deme, percent_cover, interval, tetrapod_group,
	       seed, job_type)
write.csv(best_sims, file=file.path("input", "best_sim_parameters_seeds.csv"))

