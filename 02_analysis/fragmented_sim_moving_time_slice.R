library(tidyverse)
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

# Import the coalescence simulation results


sim_results_pristine <- read.csv(file.path(results_dir_pristine, "results_.csv")) %>%
	mutate(percent_cover =100)
sim_results_fragmented <- read.csv(file.path(results_dir_fragmented, "results_fragmented.csv"))

# Import the real biodiversity data for each interval
real_biodiversity<- read.csv(file.path(data_dir, "interval_biodiversity_metrics.csv")) %>% select(-X)

# Read the interval data
interval_data <- read.csv(file.path(data_dir, "interval_data.csv")) %>% select(-X)

# Manipulate the data into the correct format
sim_results_pristine$interval <- tools::toTitleCase(as.character(sim_results_pristine$interval))
sim_results_fragmented$interval <- tools::toTitleCase(as.character(sim_results_fragmented$interval))
sim_results_fragmented$tetrapod_group <- tools::toTitleCase(
	as.character(sim_results_fragmented$tetrapod_group))
sim_results_pristine$tetrapod_group <- tools::toTitleCase(
	as.character(sim_results_pristine$tetrapod_group))

sim_results_pristine <- bind_rows(sim_results_pristine %>% mutate(percent_cover = 20)) %>%
	bind_rows(sim_results_pristine %>% mutate(percent_cover = 40)) %>%
	bind_rows(sim_results_pristine %>% mutate(percent_cover = 80)) %>%
	mutate(scenario = "pristine") %>% full_join(interval_data, by="interval") %>% na.omit()



sim_results_fragmented <- sim_results_fragmented %>% full_join(interval_data, by="interval") %>%
	na.omit() %>% mutate(scenario = "fragmented")

sim_results <- sim_results_pristine  %>% filter(midpoint >= 319) %>%
	bind_rows(sim_results_fragmented %>%
		  	filter(midpoint < 319)) %>%
	mutate(divider = 319)
for(i in seq(318, 273, -1))
{
	sim_results <- sim_results %>%
		bind_rows(sim_results_pristine  %>% filter(midpoint >= i) %>%
			  	bind_rows(sim_results_fragmented %>%
			  		  	filter(midpoint < i)) %>%
			  	mutate(divider = i))



}


# Modify the real biodiversity data
real_biodiversity <- real_biodiversity %>% full_join(interval_data, by="interval") %>%
	mutate(tetrapod_group = tools::toTitleCase(as.character(tetrapod_group)),
	       interval = tools::toTitleCase(as.character(interval)))

# Compute the mean GOF for each parameter combination, then find the top 5 parameter combinations
all_gofs <- sim_results %>% ungroup() %>%
	group_by(tetrapod_group, percent_cover, divider, speciation_rate, sigma) %>%
	summarise(mean_gof = mean(gof))

sim_results <- sim_results %>%
	group_by(tetrapod_group, percent_cover, divider, speciation_rate, sigma) %>%
	mutate(mean_gof = mean(gof))

top_gofs <- all_gofs %>% group_by(tetrapod_group, percent_cover, divider) %>%
	top_n(n = 5, wt = mean_gof) %>% rename(best_gof=mean_gof) %>%
	summarise(min_best_gof=min(best_gof))

sim_results <- sim_results %>% full_join(top_gofs)
best_5_results <- sim_results %>% na.omit(top_gofs) %>% filter(mean_gof >= min_best_gof)

save(best_5_results, file=file.path(rdata_dir, "fragmented_moving_results.RData"))

plot_top_gofs <- best_5_results %>% ungroup() %>%
	group_by(tetrapod_group, divider, percent_cover) %>%
	summarise(total_gof = mean(mean_gof),
		  min_gof = min(mean_gof),
		  max_gof = max(mean_gof))
p <- plot_top_gofs %>% ggplot() +
	geom_ribbon(aes(x=divider, ymin=min_gof, ymax=max_gof), alpha=0.5) +
	geom_line(aes(x=divider, y=total_gof))+
	facet_grid(percent_cover~tetrapod_group,
		   labeller = labeller(percent_cover=as_labeller(c("20"="80% habitat loss",
		   					      "40"="60% habitat loss",
		   					      "80" = "20% habitat loss")))) +
	theme_classic() +
	scale_y_continuous("Goodness-of-fit") +
	scale_x_reverse("Time of habitat loss (Ma)")
pdf(file.path(fig_dir, "gof_moving_slice.pdf"),  6.5, 4)
print(p)
dev.off()

# TODO finish this
best_5_results_307 <- sim_results %>% filter(divider==307) %>% na.omit() %>%
	ungroup() %>% filter(mean_gof >= min_best_gof) %>%
	gather(key="metric", value="simulation_value", richness, alpha_diversity, beta_diversity) %>%
	mutate(actual_value =
	       	ifelse(metric == "richness", actual_richness,
	       	       ifelse(metric == "alpha_diversity", actual_alpha, actual_beta))) %>%
	filter(percent_cover == 20)

ggthemr("light")
p <- best_5_results_307 %>% ggplot(aes(x=midpoint, y=simulation_value) )+
	stat_summary(fun.y = mean, geom = "line", aes(colour="Simulation\nwith CRC",
						      linetype="Simulation\nwith CRC",
						      fill="Simulation\nwith CRC")) +
	stat_summary(fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     aes(fill="Simulation\nwith CRC"),
		     colour=NA, alpha=0.3)+
	geom_line(aes(x=midpoint, y=actual_value,
		      linetype="Fossil record",
		      colour="Fossil record",
		      fill="Fossil record")) +
	theme_classic() +
	scale_linetype("")+
	scale_colour_discrete("")+
	scale_fill_discrete("")+
	# scale_colour_viridis("", option="plasma", discrete=TRUE, begin=0.2, end=0.6)+
	# scale_fill_viridis("",option="plasma", discrete=TRUE, begin=0.2, end=0.6)+
	scale_x_reverse("Time (Ma)") + scale_y_continuous("Biodiversity value") +
	geom_vline(aes(xintercept=307), linetype="dashed", colour="black")+
	facet_grid(metric~tetrapod_group, scales = "free_y",
		   labeller = labeller(metric = as_labeller(c("alpha_diversity" = "Alpha diversity",
		   					   "beta_diversity" = "Beta diversity",
		   					   "richness" = "Species richness"))))
pdf(file.path(fig_dir, "fragmented_best_all_metrics.pdf"),  6.5, 4)
print(p)
dev.off()

# Get the GOF for the single best simulation


summary_single_best <- best_5_results_307 %>% ungroup() %>%
	group_by(tetrapod_group, speciation_rate, percent_cover, deme, sigma) %>%
	summarise(best_gof = mean(mean_gof)) %>% ungroup() %>% group_by(tetrapod_group) %>%
	top_n(1, wt=best_gof)
summary_single_best

summary_single_best$best_gof
