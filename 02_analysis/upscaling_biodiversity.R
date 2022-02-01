# Up-scales biodiversity from the best-fitting simulation to see the detected level of change at within the
# sample.

library(tidyverse)
library(viridis)
library(ggthemr)

source("preston.R")

results_dir <- file.path("results")
data_dir <- "input"
fig_dir <- file.path("figures", "upscaling")
fig_dir2 <- file.path("../../../Figures/Paleo/")
for(folder in c(fig_dir, fig_dir2))
{
	if(!dir.exists(fig_dir))
	{
		dir.create(fig_dir)
	}
}
# Read the interval data
interval_data <- read.csv(file.path(data_dir, "interval_data.csv")) %>% select(-X) %>%
	mutate(interval=factor(interval, levels=c("Bashkirian", "Moscovian", "Kasimovian",
						  "Gzhelian", "Asselian", "Sakmarian",
						  "Artinskian", "Kungurian")))
params_best <- read.csv(file.path(results_dir, "Sim8", "single_param_best_5.csv")) %>%
	select(deme, sigma, tetrapod_group, interval, percent_cover)
all_sims <- read.csv(file.path(results_dir, "Sim8", "results_fragmented_finalised.csv"))

species_richness_df <- params_best %>% inner_join(all_sims)

p_richness <- species_richness_df %>%
	group_by(interval, tetrapod_group, speciation_rate, actual_richness, midpoint) %>%
	summarise(mean_richness = mean(richness),
		  sd_richness = sd(richness)) %>%
	ggplot(aes(x=midpoint, y=mean_richness)) +
	theme_classic() +
	geom_line(aes(colour=speciation_rate, group=speciation_rate, linetype="Simulated")) +
	geom_ribbon(aes(fill=speciation_rate, group=speciation_rate,
			ymin=mean_richness- sd_richness, ymax=mean_richness+sd_richness), alpha=0.3)+
	facet_grid(.~tetrapod_group) +
	scale_y_continuous("Species richness")+
	scale_x_reverse("Time (mya)") +
	scale_colour_viridis("Speciation rate", option="plasma", trans="log10", end=0.9) +
	scale_fill_viridis("Speciation rate", option="plasma", trans="log10", end=0.9)+
	geom_line(aes(x=midpoint, y=actual_richness, linetype="Fossil record")) +
	scale_linetype_discrete("")
pdf(file.path(fig_dir, "upscaling_speciation.pdf"),  5.5, 4)
print(p_richness)
dev.off()

base_local_richness <- species_richness_df %>% ungroup %>%
	group_by(interval, tetrapod_group, speciation_rate, actual_richness, midpoint) %>%
	summarise(mean_richness = mean(richness),
		  sd_richness = sd(richness)) %>%
	group_by(interval, tetrapod_group) %>%
	summarise(baseline_richness = min(mean_richness))

sampling_df <- species_richness_df %>% ungroup %>%
	group_by(interval, tetrapod_group, speciation_rate, actual_richness, midpoint) %>%
	summarise(mean_richness = mean(richness),
		  sd_richness = sd(richness)) %>%
	full_join(params_best) %>%
	full_join(base_local_richness) %>%
	mutate(a_max = 10000000,
	       area = ifelse(midpoint >= 307, a_max, a_max*(percent_cover/100)),
	       baseline_meta_richness = S_random_equilibrium(A_max = a_max,A =area,
	       					 nu = 1e-8, sigma_sq = (sigma^2)*deme),
	       meta_richness = S_random_equilibrium(A_max = a_max,A =area,
	       					 nu = speciation_rate, sigma_sq = (sigma^2)*deme),
	       meta_change = meta_richness/baseline_meta_richness,
	       local_change = mean_richness/baseline_richness,
	       detection_rate = local_change/meta_change) %>%
	rowwise() %>%
	mutate(missed_change_min = min(local_change, meta_change),
	       missed_change_max = meta_change,
	       detected_change_max = min(local_change, meta_change),
	       over_change_min = meta_change,
	       over_change_max = max(meta_change, local_change))# %>%
# Upscale with different metacommunity sizes

metacommunity_diversity <- expand.grid(tetrapod_group = c("Amniote", "Amphibian"),
				       interval = unique(species_richness_df$interval),
				       a_max = 10^seq(5, 9, 0.1))

sampling_df <- species_richness_df %>% ungroup %>%
	group_by(interval, tetrapod_group, speciation_rate, actual_richness, midpoint) %>%
	summarise(mean_richness = mean(richness),
		  sd_richness = sd(richness)) %>%
	full_join(params_best) %>%
	full_join(base_local_richness) %>%
	full_join(metacommunity_diversity) %>%
	mutate(area = ifelse(midpoint >= 307, a_max, a_max*(percent_cover/100)),
	       baseline_meta_richness = S_random_equilibrium(A_max = a_max,A =area,
	       					      nu = 1e-8, sigma_sq = (sigma^2)*deme),
	       meta_richness = S_random_equilibrium(A_max = a_max,A =area,
	       				     nu = speciation_rate, sigma_sq = (sigma^2)*deme),
	       meta_change = meta_richness/baseline_meta_richness,
	       local_change = mean_richness/baseline_richness,
	       detection_rate = local_change/meta_change) %>%
	rowwise() %>%
	mutate(missed_change_min = min(local_change, meta_change),
	       missed_change_max = meta_change,
	       detected_change_max = min(local_change, meta_change),
	       over_change_min = meta_change,
	       over_change_max = max(meta_change, local_change))

# Calculate the proportional increase in the metacommunity
metacommunity_stats <- sampling_df %>%
	select(a_max, interval, tetrapod_group, speciation_rate, meta_change, meta_richness, baseline_meta_richness) %>%
	distinct()
# Define the plot labels
plot_labels <- data.frame(
	metric = factor("alpha_diversity", levels=c("richness",
						    "alpha_diversity",
						    "beta_diversity")),
	parameter="speciation_rate",
	speciation_rate=1e-08,
	interval=interval_data$interval,
	interval_abbr=fct_recode(interval_data$interval, "Ks" = "Kasimovian", "Gz" = "Gzhelian",
				 "As" = "Asselian", "Ba"="Bashkirian", "Mo" = "Moscovian", "Sa" = "Sakmarian",
				 "Ar" = "Artinskian", "Ku" = "Kungurian"),
	midpoint=interval_data$midpoint,
	interval_part_abbr=fct_recode(interval_data$interval, "Ks" = "Kasimovian", "Gz" = "Gzhelian","As" = "Asselian")
)
ggthemr(palette="solarized")
p <- sampling_df %>% filter(a_max == 10000000,
			    speciation_rate %in% c(1e-08, 1e-06, 1e-04)) %>%
	group_by(speciation_rate, interval, tetrapod_group, midpoint) %>%
	summarise(mean_local_change = mean(local_change),
		  mean_detection_rate = mean(detection_rate),
		  sd_detection_rate = sd(detection_rate)) %>%
	ggplot()+
	geom_rect(data=interval_data,
		  aes(x=NULL, y=NULL,
		      xmin=min_ma, xmax=max_ma, ymin=0, ymax=Inf,
		      group=interval), fill=rep(c("grey60", "grey80"), 24), colour=NA, alpha=0.5)+
	geom_text(data=plot_labels, aes(x=midpoint, y = 140, label=interval_abbr),
		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+
	geom_line(aes(x=midpoint, y=mean_detection_rate*100, linetype="Fossil record\ndetection rate",
		      colour="Fossil record\ndetection rate"))+
	geom_ribbon(aes(x=midpoint, ymin=100*(mean_detection_rate - sd_detection_rate),
			ymax=100*(mean_detection_rate + sd_detection_rate)), alpha=0.2)+
	geom_hline(aes(yintercept=100.0,
		       linetype="Perfect detection\nof species richness",
		       colour="Perfect detection\nof species richness"))+
	scale_x_reverse("Time (Ma)") +
	scale_y_continuous("Sampling effectiveness of species richness (%)", breaks=c(0, 50, 100, 150)) +
	scale_linetype_discrete("")+
	scale_colour_discrete("")+
	facet_grid(speciation_rate~tetrapod_group,
		   labeller = labeller(speciation_rate = as_labeller(c("1e-08"=paste("baseline: nu==10^-8"),
		   						    "1e-06"=paste("nu==10^-6"),
		   						    "1e-04"=paste("nu==10^-4")), default=label_parsed))) +
	theme_classic()
# print(p)
pdf(file.path(fig_dir, "detection_rate.pdf"),  5.5, 4)
print(p)
dev.off()
pdf(file.path(fig_dir2, "paleo_detection_rate.pdf"),  6.5, 4)
print(p)
dev.off()
p <- sampling_df %>% filter(a_max %in% c(10^5, 10^6, 10^7, 10^8, 10^9),
			    speciation_rate %in% c(1e-08, 1e-06, 1e-04)) %>%
	group_by(a_max, speciation_rate, interval, tetrapod_group, midpoint) %>%
	summarise(mean_local_change = mean(local_change),
		  mean_detection_rate = mean(detection_rate),
		  sd_detection_rate = sd(detection_rate)) %>%
	ggplot()+
	geom_rect(data=interval_data,
		  aes(x=NULL, y=NULL,
		      xmin=min_ma, xmax=max_ma, ymin=0, ymax=Inf,
		      group=interval), fill=rep(c("grey60", "grey80"), 24), colour=NA, alpha=0.5)+
	geom_text(data=plot_labels, aes(x=midpoint, y = 175, label=interval_abbr),
		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+
	geom_line(aes(x=midpoint, y=mean_detection_rate*100,
		      colour=as.factor(a_max))) + theme_classic() +
	scale_color_viridis("Metacommunity size", discrete=TRUE,
			    option="plasma", end=0.8, labels=) +
	geom_hline(aes(yintercept=100.0,
		       linetype="Perfect detection\nof species richness"), colour="black")+
	facet_grid(speciation_rate~tetrapod_group,
		   labeller = labeller(speciation_rate = as_labeller(c("1e-08"=paste("baseline: nu==10^-8"),
		   						    "1e-06"=paste("nu==10^-6"),
		   						    "1e-04"=paste("nu==10^-4")), default=label_parsed))) +
	scale_linetype_manual("", values=c("dotted"))+
	theme_classic()+
	scale_x_reverse("Time (Ma)")+
	scale_y_continuous("Sampling effectiveness of species richness (%)", breaks=c(0, 50, 100, 150))
pdf(file.path(fig_dir, "detection_rate_all_sizes.pdf"),  5.5, 4)
print(p)
dev.off()
pdf(file.path(fig_dir2, "paleo_detection_rate_all_sizes.pdf"),  6.5, 4)
print(p)
dev.off()


