## TODO improve or remove this
## Calculate the AIC values for the species richness of the best-fitting simulations for the fragmented and pristine scenarios
library(tidyverse)
library(broom)
library(viridis)
library(ggpubr)
library(ggthemr)

# library(ggthemr)
# Set the file paths for input and output
data_dir <- "input"
results_dir_pristine <- file.path("results", "Sim5")
results_dir_fragmented <- file.path("results", "Sim6")
fig_dir <- file.path("figures", "fragmented")
rdata_dir <- "rdata"
if(!dir.exists(fig_dir))
{
	dir.create(fig_dir)
}

load(file.path(rdata_dir, "fragmented_moving_results.RData"))
load(file.path(rdata_dir, "pristine_best.RData"))
fragmented_best <- best_5_results %>% mutate(free_vars=4) %>% ungroup() %>%
	select(actual_richness, deme, gof, interval, richness, sigma, simulated_individuals,
	       speciation_rate, tetrapod_group, midpoint, mean_gof) %>%
	rename(total_gof=mean_gof) %>% mutate(scenario="fragmented")
pristine_best <- single_param_result %>% mutate(free_vars=3) %>%
	select(actual_richness, deme, gof, interval, richness, sigma, simulated_individuals,
	       speciation_rate, tetrapod_group, midpoint, total_gof) %>% mutate(scenario="pristine")

best_df <- fragmented_best %>% full_join(pristine_best)

fit_lm <- best_df %>% group_by(tetrapod_group, scenario, sigma, deme, speciation_rate) %>%
	do(rich_model = lm(richness ~ actual_richness, data = .))
head(fit_lm)
summary_fit <- glance(fit_lm, rich_model) %>% select(tetrapod_group, sigma, deme, speciation_rate, scenario, r.squared)
