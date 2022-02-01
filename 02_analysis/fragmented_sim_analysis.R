library(tidyverse)
library(viridis)
library(ggpubr)
library(ggthemr)

# library(ggthemr)
# Set the file paths for input and output
data_dir <- "input/"
results_dir_pristine <- file.path("results", "Sim7")
results_dir_fragmented <- file.path("results", "Sim8")
fig_dir <- file.path("figures", "fragmented")
fig_dir2 <- file.path("../../../Figures/Paleo")
rdata_dir <- "rdata"
for(folder in c(fig_dir, fig_dir2))
{
	if(!dir.exists(folder))
	{
		dir.create(folder)
	}
}
# Import the coalescence simulation results


sim_results_pristine <- read.csv(file.path(results_dir_pristine, "results_.csv")) %>%
	mutate(percent_cover =100)
sim_results_fragmented <- read.csv(file.path(results_dir_fragmented, "results_fragmented.csv")) %>%
	mutate(deme=floor(deme))
sim_results_pristine$interval <- tools::toTitleCase(as.character(sim_results_pristine$interval))
sim_results_fragmented$interval <- tools::toTitleCase(as.character(sim_results_fragmented$interval))
sim_results_fragmented$tetrapod_group <- tools::toTitleCase(as.character(sim_results_fragmented$tetrapod_group))
sim_results_pristine$tetrapod_group <- tools::toTitleCase(as.character(sim_results_pristine$tetrapod_group))

sim_results_pristine <- sim_results_pristine %>% mutate(scenario = "pristine")

# Import the real biodiversity data for each interval
real_biodiversity<- read.csv(file.path(data_dir, "interval_biodiversity_metrics.csv")) %>%
	select(-X) %>%
	mutate(tetrapod_group=tools::toTitleCase(as.character(tetrapod_group)))

# Read the interval data
# Read the interval data
interval_data <- read.csv(file.path(data_dir, "interval_data.csv")) %>% select(-X) %>%
	mutate(interval=factor(interval, levels=c("Bashkirian", "Moscovian", "Kasimovian",
						  "Gzhelian", "Asselian", "Sakmarian",
						  "Artinskian", "Kungurian")))
sim_results <- sim_results_fragmented %>% mutate(scenario="fragmented") %>%
	bind_rows(sim_results_pristine %>% mutate(percent_cover = 20)) %>%
	bind_rows(sim_results_pristine %>% mutate(percent_cover = 40)) %>%
	bind_rows(sim_results_pristine %>% mutate(percent_cover = 80)) %>%
	mutate(deme=round(deme)) %>%
	full_join(interval_data, by="interval") %>%
	full_join(real_biodiversity %>% select(interval, tetrapod_group, number_individuals)) %>%
	filter((midpoint >= 305 & scenario=="pristine") |
	       	(midpoint < 305 & scenario == "fragmented"), # this is 305 so that only times which are entirely post-CRC are included
	       # simulated_individuals < 5*number_individuals) %>%
	) %>%
	group_by(tetrapod_group, sigma, deme, speciation_rate, percent_cover) %>%
	mutate(number_sims = n()) %>% filter(number_sims >=7) %>%  ungroup()

write.csv(sim_results, file.path(results_dir_fragmented, "results_fragmented_finalised.csv"))

# Modify the real biodiversity data
real_biodiversity <- real_biodiversity %>% full_join(interval_data, by="interval") %>%
	mutate(tetrapod_group = tools::toTitleCase(as.character(tetrapod_group)),
	       interval = tools::toTitleCase(as.character(interval)))

## METHOD B
# Look at the top 5 fitting results for each tetrapod group
single_param_result_separate <- sim_results %>% ungroup() %>%
	group_by(sigma, speciation_rate, deme, tetrapod_group, percent_cover) %>%
	mutate(total_gof = mean(gof))

# Compute some stats about the parameters for the top 3 fitting simulations
summary_stats_separate <- single_param_result_separate  %>%
	select(sigma, speciation_rate, deme, total_gof, tetrapod_group, percent_cover) %>% distinct() %>%
	ungroup() %>% group_by(tetrapod_group) %>%
	top_n(total_gof, n=5) %>%
	gather(key="parameter", value="value", sigma, speciation_rate, deme, total_gof) %>%
	group_by(parameter, tetrapod_group) %>%
	summarise(mean_value = mean(value),
		  min_value = min(value),
		  max_value = max(value),
		  sd_value = sd(value),
		  total=n())
top_gofs <- single_param_result_separate %>% ungroup() %>% select(total_gof, tetrapod_group, percent_cover)%>%
	distinct() %>% group_by(tetrapod_group, percent_cover) %>% top_n(n=5, wt=total_gof) %>%
	summarise(min=min(total_gof))

# Get the top 5 best-fitting results
single_param_result_separate <- single_param_result_separate %>% ungroup() %>%
	group_by(tetrapod_group, percent_cover) %>% rowwise() %>%
	mutate(top_3 = total_gof >= top_gofs[top_gofs$tetrapod_group == tetrapod_group,]$min[1]) %>%
	filter(top_3)

# DECIDE WHICH METHOD
# Comment this out if you want to use the same parameters for amphibians and amniotes
# So with this line, use METHOD B for analysis
single_param_result <- single_param_result_separate


# Plot the landscape richness over time, compared to the real species richness over time
p <- ggplot(sim_results, aes(x=midpoint, y=richness, group=1))  +
	theme_classic() + ylab("Species richness") + xlab("Interval") +
	stat_summary(fun.y = mean, geom = "line", aes(linetype="Simulated",
						      colour=tetrapod_group,
						      group=tetrapod_group)) +
	stat_summary(fun.data = mean_se, geom = "errorbar")+geom_point()+
	geom_line(data=real_biodiversity, aes(x=midpoint, colour=tetrapod_group,
					      y=richness, group=tetrapod_group,linetype="Fossil record")) +
	scale_linetype_discrete("Data source") + scale_x_reverse()
# Combine all three into a single plot
# Define the plot labels
plot_labels <- data.frame(
	metric = factor("alpha_diversity", levels=c("richness",
						    "alpha_diversity",
						    "beta_diversity")),
	parameter="speciation_rate",
	interval=interval_data$interval,
	interval_abbr=fct_recode(interval_data$interval, "Ks" = "Kasimovian", "Gz" = "Gzhelian",
				 "As" = "Asselian", "Ba"="Bashkirian", "Mo" = "Moscovian", "Sa" = "Sakmarian",
				 "Ar" = "Artinskian", "Ku" = "Kungurian"),
	midpoint=interval_data$midpoint,
	interval_part_abbr=fct_recode(interval_data$interval, "Ks" = "Kasimovian", "Gz" = "Gzhelian","As" = "Asselian")
)
p_all <- single_param_result %>% gather(key="metric", value="value",
				       beta_diversity, alpha_diversity, richness) %>%
	mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
				   ifelse(metric == "beta_diversity",
				          actual_beta, actual_richness))) %>%
	ggplot() + theme_classic() +
	geom_rect(data=interval_data,
		  aes(x=NULL, y=NULL,
		      xmin=min_ma, xmax=max_ma, ymin=0, ymax=Inf,
		      group=interval), fill=rep(c("grey60", "grey80"), 24), colour=NA, alpha=0.5)+
	geom_text(data=plot_labels, aes(x=midpoint, y = 3.7, label=interval_abbr),
		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+

	ylab("Biodiversity value") + xlab("Time (Ma)") +
	stat_summary(fun.y = mean, geom = "line", aes(linetype="Simulated",
						      colour=as.factor(percent_cover),
						      group=as.factor(percent_cover),
						      x=midpoint, y=value)) +
	stat_summary(fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     colour=NA,
		     aes(fill=as.factor(percent_cover),
		         group=as.factor(percent_cover),
		         x=midpoint, y=value), alpha=0.3)+

	geom_line(aes(x=midpoint, y=value_real, linetype="Fossil record")) +
	scale_colour_viridis("Remaining\nhabitat %", discrete=TRUE, begin=0.2, end=0.65, option="plasma") +
	scale_fill_viridis("Remaining\nhabitat %", discrete=TRUE, begin=0.2, end=0.65, option="plasma") +
	scale_linetype_discrete("") + scale_x_reverse() +
	facet_grid(metric~tetrapod_group, scales = "free_y",
		   labeller = labeller(metric = as_labeller(c("alpha_diversity" = "Alpha diversity",
		   					   "beta_diversity" = "Beta diversity",
		   					   "richness" = "Species richness"))))+
	guides(linetype=guide_legend(override.aes = list(colour=c("#62bba5", "black"))))
print(p_all)
pdf(file.path(fig_dir, "fragmented_best_all_metrics_all_cover.pdf"), 6.5, 4)
print(p_all)
dev.off()
pdf(file.path(fig_dir2, "fragmented_best_all_metrics_all_cover.pdf"), 6.5, 4)
print(p_all)
dev.off()


# # Find the single one best parameter set across sims for output
write.csv(single_param_result, file.path(results_dir_fragmented, "single_param_best_5.csv"))
ggthemr_reset()
ggthemr("light")
pallete <- swatch()
pallete[3] <- pallete[9]
p_all <- single_param_result %>% filter(percent_cover==20) %>%
	gather(key="metric", value="value",
	       beta_diversity, alpha_diversity, richness) %>%
	mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
				   ifelse(metric == "beta_diversity",
				          actual_beta, actual_richness))) %>%
	ggplot() +
	theme_classic() +
	ylab("Biodiversity value") + xlab("Time (Ma)") +
	geom_rect(data=interval_data,
		  aes(x=NULL, y=NULL,
		      xmin=min_ma, xmax=max_ma, ymin=0, ymax=Inf,
		      group=interval), fill=rep(c("grey60", "grey80"), 24), colour=NA, alpha=0.5)+
	stat_summary(fun.y = mean, geom = "line",
		     aes(linetype="Simulations\n(with CRC)",
		         colour="Simulations\n(with CRC)",
		         fill="Simulations\n(with CRC)",
		         x=midpoint, y=value)) +
	stat_summary(fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     colour=NA,
		     aes(fill="Simulations\n(with CRC)",
		         linetype="Simulations\n(with CRC)",
		         colour="Simulations\n(with CRC)",
		         x=midpoint, y=value), alpha=0.3)+
	geom_line(aes(x=midpoint, y=value_real, colour="Fossil record",
		      linetype="Fossil record",
		      fill="Fossil record")) +
	scale_colour_manual("", values=pallete[2:3]) +
	scale_linetype("")+
	scale_fill_manual("", values=pallete[2:3]) +
	scale_x_reverse() +
	geom_text(data=plot_labels, aes(x=midpoint, y = 3.7, label=interval_abbr),
		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+
	facet_grid(metric~tetrapod_group, scales = "free_y",
		   labeller = labeller(metric = as_labeller(c("alpha_diversity" = "Alpha diversity",
		   					   "beta_diversity" = "Beta diversity",
		   					   "richness" = "Species richness")))) +
	geom_vline(linetype="dashed", colour="black", aes(xintercept=307))
pdf(file.path(fig_dir, "fragmented_best_all_metrics_20_cover.pdf"), 6.5, 4)
print(p_all)
dev.off()
pdf(file.path(fig_dir2, "fragmented_best_all_metrics_20_cover.pdf"), 6.5, 4)
print(p_all)
dev.off()

by_interval <- single_param_result %>% distinct(actual_alpha, actual_beta, actual_richness, simulated_individuals,
				 deme, number_individuals, interval, tetrapod_group, midpoint)


