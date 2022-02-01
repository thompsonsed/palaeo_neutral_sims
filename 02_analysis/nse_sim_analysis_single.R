# Analyse the spatially implicit simulations assuming a single well-mixed community of individuals.

library(tidyverse)
library(viridis)

data_dir <- "input"
fig_dir <- file.path("figures", "nse")
results_dir <- file.path("results", "nse1_single")
rdata_dir <- "rdata"
if(!dir.exists(fig_dir))
{
	dir.create(fig_dir)
}

# Read in the real data
interval_data <- read.csv(file.path(data_dir, "interval_data.csv")) %>% select(-X) %>%
	mutate(interval=tolower(interval))

info_per_pcoord<- read.csv(file.path(data_dir, "info_per_pcoord_main.csv")) %>% select(-X) %>%
	rename(fragment=fragment_name, real_richness = individuals_total) %>%
	mutate(interval = tolower(interval), fragment=as.character(fragment)) %>%
	select(-species_total, -pcoords, -collections)

real_biodiversity<- read.csv(file.path(data_dir, "interval_biodiversity_metrics.csv")) %>%
	select(-X) %>% mutate(interval = tolower(interval)) %>%
	rename(alpha_real = alpha, beta_real = beta, richness_real = richness)

if(!file.exists(file.path(rdata_dir, "nse_single.RData")))
{
	nse_fragment_richness <- read.csv(file.path(results_dir, "nse_fragment_richness1.csv")) %>%
		select(-X) %>% mutate(fragment = as.character(fragment)) %>%
		rename(immigration_rate = speciation_rate,
		       speciation_rate = metacommunity_speciation_rate)
	nse_alpha_beta <- read.csv(file.path(results_dir, "nse_metrics1.csv")) %>% select(-X) %>%
		rename(immigration_rate = speciation_rate,
		       speciation_rate = metacommunity_speciation_rate) %>%
		filter(metacommunity_size > 0)


	# Fragment richness values first, then combined with the other metrics

	nse_combined <- nse_fragment_richness %>% full_join(info_per_pcoord) %>%
		mutate(richness_error = abs(richness - real_richness)/pmax(richness, real_richness)) %>%
		group_by(interval, tetrapod_group, seed,
			 speciation_rate, immigration_rate, metacommunity_size) %>%
		summarise(mean_error = mean(richness_error))


	nse_metrics <- nse_alpha_beta %>% full_join(real_biodiversity) %>%
		full_join(nse_combined) %>%
		rename(error_fragment_richness=mean_error) %>%
		mutate(error_richness =
		       	abs(species_richness-richness_real)/pmax(species_richness, richness_real),
		       error_alpha = abs(alpha-alpha_real)/pmax(alpha, alpha_real),
		       error_beta = abs(beta-beta_real)/pmax(beta, beta_real),
		       mean_error = (error_richness + error_alpha +
		       	      	error_beta+error_fragment_richness)/4) %>%
		group_by(interval, tetrapod_group, speciation_rate,
			 immigration_rate, metacommunity_size, alpha_real, beta_real) %>%
		# measure the total error across repeat sims
		summarise(mean_total_error = mean(mean_error),
			  mean_richness = mean(species_richness),
			  mean_alpha = mean(alpha),
			  mean_beta=mean(beta))


	nse_grouped_parameters <- nse_metrics %>%
		group_by(speciation_rate, immigration_rate, metacommunity_size) %>%
		summarise(mean_error = mean(mean_total_error))

	top_simulations <- nse_grouped_parameters %>% ungroup() %>%
		top_n(n=5, wt = 1-mean_error)

	# Plot diversity over time from the nse model
	best_simulations_richness <- top_simulations %>% inner_join(nse_metrics) %>%
		full_join(real_biodiversity) %>%
		full_join(interval_data)
	save(nse_metrics, nse_grouped_parameters, top_simulations, best_simulations_richness,
	     file=file.path(rdata_dir, "nse_single.RData"))

}else
{
	load(file.path(rdata_dir, "nse_single.RData"))
}


p_all <- best_simulations_richness %>% gather(key="metric", value="value",
					      mean_beta, mean_alpha, mean_richness) %>%
	mutate(value_real = ifelse(metric == "mean_alpha", alpha_real,
				   ifelse(metric == "mean_beta",
				          beta_real, richness_real))) %>%
	ggplot(aes(x=midpoint, y=value)) + theme_classic() +
	ylab("Biodiversity metric") + xlab("Ma") +
	stat_summary(fun.y = mean, geom = "line", aes(linetype="Simulated",
						      colour="Simulated")) +
	stat_summary(fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     colour=NA,
		     aes(fill="Simulated"), alpha=0.3)+
	geom_line(aes(x=midpoint, y=value_real, linetype="Fossil record",
		      colour="Fossil record", fill="Fossil record")) +
	scale_colour_viridis("Data source", discrete=TRUE, begin=0.2, end=0.65, option="plasma") +
	scale_fill_viridis("Data source", discrete=TRUE, begin=0.2, end=0.65, option="plasma") +
	scale_linetype_discrete("Data source") + scale_x_reverse() +
	facet_grid(metric~tetrapod_group, scales = "free_y",
		   labeller = labeller(metric = as_labeller(c("mean_alpha" = "Alpha diversity",
		   					   "mean_beta" = "Beta diversity",
		   					   "mean_richness" = "Species richness")),
		   		    tetrapod_group=as_labeller(c("amniote" = "Amniote",
		   		    			     "amphibian" = "Amphibian"))))
pdf(file.path(fig_dir, "nse_metrics_single_all.pdf"), 6, 5)
print(p_all)
dev.off()

# Generate some summary statistics

summary_df <- best_simulations_richness %>% ungroup() %>% group_by(tetrapod_group) %>%
	summarise(total_gof = 1-mean(mean_error))
