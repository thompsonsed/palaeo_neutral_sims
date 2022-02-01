library(tidyverse)
library(viridis)
library(akima)
library(interp)
library(ggpubr)
library(ggthemr)
# Import the coalescence simulation results
data_dir <- "input"
results_dir <- file.path("results", "Sim9")
rdata_dir <- "rdata"
fig_dir <- file.path("figures", "main")
fig_dir2 <- file.path("../../../Figures/Paleo/")
for(i in c(rdata_dir, fig_dir, fig_dir2))
{
	if(!dir.exists(i))
	{
		dir.create(i)
	}
}
sim_results <- read.csv(file.path(results_dir, "results__occ_sq.csv"))
sim_results$interval <- tools::toTitleCase(as.character(sim_results$interval))
sim_results$tetrapod_group <- tools::toTitleCase(as.character(sim_results$tetrapod_group))
# Import the real biodiversity data for each interval
real_biodiversity<- read.csv(file.path(data_dir, "interval_biodiversity_metrics_occ.csv")) %>% select(-X)

# Read the interval data
interval_data <- read.csv(file.path(data_dir, "interval_data.csv")) %>% select(-X) %>%
	mutate(interval=factor(interval, levels=c("Bashkirian", "Moscovian", "Kasimovian",
						  "Gzhelian", "Asselian", "Sakmarian",
						  "Artinskian", "Kungurian")))
sim_results <- sim_results %>% full_join(interval_data, by="interval")
# Compute the top 3 results for each interval and tetrapod group independently
result <- sim_results %>%
	group_by(interval, tetrapod_group) %>% top_n(n=3, wt=gof)


# METHOD A
# Find the best-fitting simulation with the same parameters across tetrapod groups and intervals
single_param_result <- sim_results %>% group_by(sigma, speciation_rate, deme) %>%
	mutate(total_gof = mean(gof))
# Compute some stats about the parameters for the top 10 fitting simulations
summary_stats <- single_param_result  %>%
	select(sigma, speciation_rate, deme, total_gof) %>% distinct() %>% ungroup() %>%
	top_n(total_gof, n=10) %>%
	gather(key="parameter", value="value", sigma, speciation_rate, deme) %>%
	group_by(parameter) %>%
	summarise(mean_value = mean(value),
		  min_value = min(value),
		  max_value = max(value),
		  sd_value = sd(value),
		  total=n())
write.csv(summary_stats, file=file.path(results_dir, "summary_stats_occ.csv"))
single_param_result <- single_param_result %>% ungroup() %>% top_n(n=1, wt=total_gof)


## METHOD B
# Look at the top 3 fitting results for each tetrapod group
single_param_result_separate <- sim_results %>%
	group_by(sigma, speciation_rate, deme, tetrapod_group) %>%
	mutate(total_gof = mean(gof))
# Compute some stats about the parameters for the top 3 fitting simulations
summary_stats_separate <- single_param_result_separate  %>%
	select(sigma, speciation_rate, deme, total_gof, tetrapod_group) %>% distinct() %>%
	ungroup() %>% group_by(tetrapod_group) %>%
	top_n(total_gof, n=5) %>%
	gather(key="parameter", value="value", sigma, speciation_rate, deme, total_gof) %>%
	group_by(parameter, tetrapod_group) %>%
	summarise(mean_value = mean(value),
		  min_value = min(value),
		  max_value = max(value),
		  sd_value = sd(value),
		  total=n())
write.csv(summary_stats_separate, file=file.path(results_dir, "summary_stats_separate_occ.csv"))
# Calculate the minimum goodness-of-fit in the top 5 fitting sims for each tetrapod group
top_gofs <- single_param_result_separate %>% ungroup() %>% select(total_gof, tetrapod_group)%>%
	distinct() %>% group_by(tetrapod_group) %>% top_n(n=5, wt=total_gof) %>%
	summarise(min=min(total_gof))

# Get the top 5 best-fitting results
single_param_result_separate <- single_param_result_separate %>% ungroup() %>%
	group_by(tetrapod_group) %>% rowwise() %>%
	mutate(top_3 = total_gof >= top_gofs[top_gofs$tetrapod_group == tetrapod_group,]$min[1]) %>%
	filter(top_3)

# DECIDE WHICH METHOD
# Comment this out if you want to use the same parameters for amphibians and amniotes
# So with this line, use METHOD B for analysis
single_param_result <- single_param_result_separate


save(single_param_result, file=file.path(rdata_dir, "pristine_best_occ.RData"))

# Modify the real biodiversity data
real_biodiversity <- real_biodiversity %>% full_join(interval_data, by="interval") %>%
	mutate(tetrapod_group = tools::toTitleCase(as.character(tetrapod_group)),
	       interval = tools::toTitleCase(as.character(interval)))

# Read in the fragment abundances
fragment_richness_rdata_path <- file.path(rdata_dir, "fragment_richness_occ.RData")
if(file.exists(fragment_richness_rdata_path))
{
	load(fragment_richness_rdata_path)
}else{
	fragment_richness <- read_csv(file.path(results_dir, "results_fragment_abundances__occ_sq.csv"),
				      col_types = cols(fragment = col_character())) %>%
		select(-X1) %>% mutate(interval=tools::toTitleCase(interval),
				       tetrapod_group=tools::toTitleCase(tetrapod_group))
	save(fragment_richness, file=fragment_richness_rdata_path)
}

real_fragment_richness <- read_csv(file.path(data_dir, "info_per_pcoord_main.csv"),
				   col_types = cols(fragment_name = col_character())) %>%
	select(-X1, -pcoords, -collections) %>%
	mutate(tetrapod_group = tools::toTitleCase(tetrapod_group)) %>%
	rename(real_richness = species_total, fragment = fragment_name)
fragment_richness <- fragment_richness %>% full_join(real_fragment_richness)

best_sim_fragment_richness <-fragment_richness %>%
	full_join(single_param_result %>% select(interval, deme, sigma, speciation_rate,
						 tetrapod_group, total_gof)) %>% na.omit()

# Species-area relationship
p4 <- best_sim_fragment_richness %>% ggplot() +
	xlab("No. individuals") + ylab("Species richness") + theme_classic()+
	facet_grid(.~tetrapod_group) +
	geom_point(alpha=0.2, aes(x=individuals_total, y=richness, colour="Simulation"))+
	geom_point(alpha=0.2, aes(x=individuals_total, y=real_richness, colour="Fossil record"))+
	geom_smooth(fill=NA, aes(x=individuals_total, y=richness, colour="Simulation"), method="lm")+
	geom_smooth(fill=NA, aes(x=individuals_total, y=real_richness, colour="Fossil record"), method="lm")+
	scale_colour_viridis("", discrete=TRUE,
			     begin=0.1, end=0.7, option="plasma") +
	# scale_linetype_discrete("Data source",
				# guide = guide_legend(override.aes = list(color = "black"))) +
	theme(aspect.ratio = 1)

pdf(file.path(fig_dir, "species_area_best_occ.pdf"), 6, 3)
print(p4)
dev.off()

# Create a data frame containing information about the two tetrapod groups combined

all_tet_groups <- single_param_result%>%
	gather(key="metric", value="value",
	       beta_diversity, alpha_diversity, richness) %>%
	mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
				   ifelse(metric == "beta_diversity",
				          actual_beta, actual_richness))) %>%
	group_by(interval, tetrapod_group, metric, value_real) %>%
	summarise(mean_val=mean(value),
		  min_val = min(value),
		  max_val = max(value))


# Plot all the metrics
ggthemr_reset()
ggthemr("light")
pallete <- swatch()
pallete[3] <- pallete[7]
#
# p_all <- single_param_result %>% gather(key="metric", value="value",
# 					beta_diversity, alpha_diversity, richness) %>%
# 	mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
# 				   ifelse(metric == "beta_diversity",
# 				          actual_beta, actual_richness))) %>%
# 	ggplot(aes(x=midpoint, y=value)) + theme_classic() +
# 	ylab("Biodiversity value") + xlab("Time (Ma)") +
# 	stat_summary(fun.y = mean, geom = "line",
# 		     aes(linetype="Simulation\n(without CRC)", colour="Simulation\n(without CRC)",
# 		         fill="Simulation\n(without CRC)")) +
# 	stat_summary(fun.ymin = min,
# 		     fun.ymax = max,
# 		     geom = "ribbon",
# 		     colour=NA,
# 		     aes(fill="Simulation\n(without CRC)", linetype="Simulation\n(without CRC)",
# 		         colour="Simulation\nwithout CRC"), alpha=0.3)+
# 	geom_line(aes(x=midpoint, y=value_real, colour="Fossil record",
# 		      linetype="Fossil record",
# 		      fill="Fossil record")) +
# 	scale_colour_manual("", values = pallete[2:3]) +
# 	scale_linetype("")+
# 	scale_fill_manual("", values = pallete[2:3])+
# 	scale_x_reverse() +
# 	facet_grid(metric~tetrapod_group, scales = "free_y",
# 		   labeller = labeller(metric = as_labeller(c("alpha_diversity" = "Alpha diversity",
# 		   					   "beta_diversity" = "Beta diversity",
# 		   					   "richness" = "Species richness"))))
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
	ggplot() +
	geom_rect(data=interval_data,
		  aes(x=NULL, y=NULL,
		      xmin=min_ma, xmax=max_ma, ymin=0, ymax=Inf,
		      group=interval), fill=rep(c("grey60", "grey80"), 24), colour=NA, alpha=0.5)+
	theme_classic() +
	ylab("Biodiversity value") + xlab("Time (Ma)") +
	stat_summary(fun.y = mean, geom = "line",
		     aes(linetype="Simulation\n(without CRC)", colour="Simulation\n(without CRC)",
		         fill="Simulation\n(without CRC)",x=midpoint, y=value)) +
	stat_summary(fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     colour=NA,
		     aes(x=midpoint, y=value,
		         fill="Simulation\n(without CRC)", linetype="Simulation\n(without CRC)",
		         colour="Simulation\nwithout CRC"), alpha=0.3)+
	geom_text(data=plot_labels, aes(x=midpoint, y = 3.7, label=interval_abbr),
		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+
	geom_line(aes(x=midpoint, y=value_real, colour="Fossil record",
		      linetype="Fossil record",
		      fill="Fossil record")) +
	scale_colour_manual("", values = pallete[2:3]) +
	scale_linetype("")+
	scale_fill_manual("", values = pallete[2:3])+
	scale_x_reverse() +
	facet_grid(metric~tetrapod_group, scales = "free_y",
		   labeller = labeller(metric = as_labeller(c("alpha_diversity" = "Alpha diversity",
		   					   "beta_diversity" = "Beta diversity",
		   					   "richness" = "Species richness"))))
print(p_all)
pdf(file.path(fig_dir, "main_metrics_all_occ.pdf"),  5.5, 4)
print(p_all)
dev.off()
pdf(file.path(fig_dir2, "main_metrics_all_occ.pdf"),  6.5, 4)
print(p_all)
dev.off()


# Analyse just the late-stage parameters for best-fitting simulations

late_results <- sim_results %>% mutate(post_307 = midpoint < 307) %>%
	group_by(sigma, speciation_rate, deme, tetrapod_group, post_307) %>%
	mutate(total_gof = mean(gof))

# Compute some stats about the parameters for the top 3 fitting simulations
summary_stats_late <- late_results  %>%
	select(sigma, speciation_rate, deme, total_gof, tetrapod_group, post_307) %>% distinct() %>%
	ungroup() %>% group_by(tetrapod_group, post_307) %>%
	top_n(total_gof, n=5) %>%
	gather(key="parameter", value="value", sigma, speciation_rate, deme, total_gof) %>%
	group_by(parameter, tetrapod_group, post_307) %>%
	summarise(mean_value = mean(value),
		  min_value = min(value),
		  max_value = max(value),
		  sd_value = sd(value),
		  total=n())
write.csv(summary_stats_late, file=file.path(results_dir, "summary_stats_late_occ.csv"))
# Calculate the minimum goodness-of-fit in the top 5 fitting sims for each tetrapod group
top_gofs <- late_results %>% ungroup() %>% select(total_gof, tetrapod_group, post_307)%>%
	distinct() %>% group_by(tetrapod_group, post_307) %>% top_n(n=3, wt=total_gof) %>%
	summarise(min=min(total_gof))
late_results <- late_results %>% ungroup() %>%
	group_by(tetrapod_group) %>% rowwise() %>%
	mutate(top_5 = total_gof >= top_gofs[top_gofs$tetrapod_group == tetrapod_group &
					     	top_gofs$post_307 == post_307,]$min[1]) %>%
	filter(top_5) %>%
	filter((post_307 & midpoint < 307) | (!post_307 & midpoint >= 307))

p <- late_results %>% gather(key="metric", value="value",
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
		     aes(linetype="Simulation predictions\nfrom Carboniferous",
		         colour="Simulation predictions\nfrom Carboniferous",
		         fill="Simulation predictions\nfrom Carboniferous",
		         x=midpoint, y=value)) +
	stat_summary(fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     colour=NA,
		     aes(fill="Simulation predictions\nfrom Carboniferous",
		         linetype="Simulation predictions\nfrom Carboniferous",
		         colour="Simulation predictions\nfrom Carboniferous",
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

# p
# Plot all the metrics


# p_all <- late_results %>% gather(key="metric", value="value",
# 					beta_diversity, alpha_diversity, richness) %>%
# 	mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
# 				   ifelse(metric == "beta_diversity",
# 				          actual_beta, actual_richness))) %>%
# 	ggplot(aes()) + theme_classic() +
# 	ylab("Biodiversity value") + xlab("Time (Ma)") +
# 	geom_rect(data=interval_data,
# 		  aes(x=NULL, y=NULL,
# 		      xmin=min_ma, xmax=max_ma, ymin=0, ymax=Inf,
# 		      group=interval), fill=rep(c("grey60", "grey80"), 24), colour=NA, alpha=0.5)+
# 	stat_summary(fun.y = mean, geom = "line",
# 		     aes(linetype="Simulation\n(split parameters)", colour="Simulation\n(split parameters)",
# 		         fill="Simulation\n(split parameters)",x=midpoint, y=value)) +
# 	stat_summary(fun.ymin = min,
# 		     fun.ymax = max,
# 		     geom = "ribbon",
# 		     colour=NA,
# 		     aes(fill="Simulation\n(split parameters)", linetype="Simulation\n(split parameters)",
# 		         colour="Simulation\n(split parameters)", x=midpoint, y=value), alpha=0.3)+
# 	geom_line(aes(x=midpoint, y=value_real, colour="Fossil record",
# 		      linetype="Fossil record",
# 		      fill="Fossil record")) +
# 	geom_text(data=plot_labels, aes(x=midpoint, y = 3.7, label=interval_abbr),
# 		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+
# 	scale_colour_manual("", values=pallete[2:3]) +
# 	scale_linetype("")+
# 	scale_fill_manual("", values=pallete[2:3])+
# 	scale_x_reverse() +
# 	facet_grid(metric~tetrapod_group, scales = "free_y",
# 		   labeller = labeller(metric = as_labeller(c("alpha_diversity" = "Alpha diversity",
# 		   					   "beta_diversity" = "Beta diversity",
# 		   					   "richness" = "Species richness")))) +
# 	geom_vline(linetype="dashed", colour="black", aes(xintercept=307))
#
# # print(p_all)
# pdf(file.path(fig_dir, "main_metrics_all_split.pdf"),  5.5, 4)
# print(p_all)
# dev.off()


## Parameterise our model based on early tetrapod diversity, and use this to predict late tetrapod diversity

early_params <- sim_results %>% group_by(sigma, speciation_rate, deme, tetrapod_group) %>%
	filter(midpoint > 307) %>%
	mutate(total_gof = mean(gof))%>%
	ungroup() %>% select(deme, sigma, speciation_rate, total_gof, tetrapod_group) %>%
	distinct()

# Calculate the minimum goodness-of-fit in the top 5 fitting sims for each tetrapod group
top_gofs_early <- early_params %>%
	select(total_gof, tetrapod_group)%>%
	distinct() %>% group_by(tetrapod_group) %>% top_n(n=5, wt=total_gof) %>%
	summarise(min=min(total_gof))

early_results <- early_params %>% rowwise() %>%
	mutate(top_3 = total_gof >= top_gofs_early[top_gofs_early$tetrapod_group == tetrapod_group,]$min[1]) %>%
	filter(top_3)  %>%
	inner_join(sim_results)

actual_gofs <- early_results %>% group_by(sigma, speciation_rate, deme, tetrapod_group) %>%
	summarise(total_gof = mean(gof))
# Get new colours for the plot
ggthemr_reset()
ggthemr("light")
pallete <- swatch()
pallete[3] <- pallete[5]

p_all <- early_results %>% gather(key="metric", value="value",
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
		     aes(linetype="Simulation predictions\nfrom Carboniferous",
		         colour="Simulation predictions\nfrom Carboniferous",
		         fill="Simulation predictions\nfrom Carboniferous",
		         x=midpoint, y=value)) +
	stat_summary(fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     colour=NA,
		     aes(fill="Simulation predictions\nfrom Carboniferous",
		         linetype="Simulation predictions\nfrom Carboniferous",
		         colour="Simulation predictions\nfrom Carboniferous",
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
pdf(file.path(fig_dir, "main_metrics_early_occ.pdf"),  5.5, 4)
print(p_all)
dev.off()

pdf(file.path(fig_dir2, "main_metrics_early_occ.pdf"),  6.5, 4)
print(p_all)
dev.off()


## Parameterise our model based on late tetrapod diversity, and use this to predict late tetrapod diversity

late_params <- sim_results %>% group_by(sigma, speciation_rate, deme, tetrapod_group) %>%
	filter(midpoint < 307) %>%
	mutate(total_gof = mean(gof))%>%
	ungroup() %>% select(deme, sigma, speciation_rate, total_gof, tetrapod_group) %>%
	distinct()

# Calculate the minimum goodness-of-fit in the top 5 fitting sims for each tetrapod group

variable_results <- early_results %>% filter(midpoint > 307) %>%
	bind_rows(late_results %>% filter(midpoint < 307))

# Get new colours for the plot
ggthemr_reset()
ggthemr("light")
pallete <- swatch()
pallete[3] <- pallete[8]
p_all <- variable_results %>%
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
		     aes(linetype="Simulations using\nvariable parameters",
		         colour="Simulations using\nvariable parameters",
		         fill="Simulations using\nvariable parameters",
		         x=midpoint, y=value)) +
	stat_summary(fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     colour=NA,
		     aes(fill="Simulations using\nvariable parameters",
		         linetype="Simulations using\nvariable parameters",
		         colour="Simulations using\nvariable parameters",
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
pdf(file.path(fig_dir, "main_metrics_variable_occ.pdf"),  5.5, 4)
print(p_all)
dev.off()

pdf(file.path(fig_dir2, "main_metrics_variable_occ.pdf"),  6.5, 4)
print(p_all)
dev.off()
# p_all
#
#
# ## Calculate the best-fitting parameter sets over time
#
# sim_results_variable <- sim_results %>% group_by(sigma, speciation_rate, deme, tetrapod_group, interval) %>%
# 	mutate(total_gof= mean(gof))
#
# top_gofs_variable <- sim_results_variable %>% ungroup() %>%
# 	group_by(interval, tetrapod_group, sigma, speciation_rate, deme) %>%
# 	summarise(total_gof=mean(gof)) %>% group_by(tetrapod_group, interval) %>%
# 	top_n(n=5, wt=total_gof)
#
# sim_results_variable <- sim_results %>% inner_join(top_gofs_variable)
#
# p_all <- sim_results_variable %>% gather(key="metric", value="value",
# 				 beta_diversity, alpha_diversity, richness) %>%
# 	mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
# 				   ifelse(metric == "beta_diversity",
# 				          actual_beta, actual_richness))) %>%
# 	ggplot() +
# 	theme_classic() +
# 	ylab("Biodiversity value") + xlab("Time (Ma)") +
# 	geom_rect(data=interval_data,
# 		  aes(x=NULL, y=NULL,
# 		      xmin=min_ma, xmax=max_ma, ymin=0, ymax=Inf,
# 		      group=interval), fill=rep(c("grey60", "grey80"), 24), colour=NA, alpha=0.5)+
# 	stat_summary(fun.y = mean, geom = "line",
# 		     aes(linetype="Simulation predictions\nfrom Carboniferous",
# 		         colour="Simulation predictions\nfrom Carboniferous",
# 		         fill="Simulation predictions\nfrom Carboniferous",
# 		         x=midpoint, y=value)) +
# 	stat_summary(fun.ymin = min,
# 		     fun.ymax = max,
# 		     geom = "ribbon",
# 		     colour=NA,
# 		     aes(fill="Simulation predictions\nfrom Carboniferous",
# 		         linetype="Simulation predictions\nfrom Carboniferous",
# 		         colour="Simulation predictions\nfrom Carboniferous",
# 		         x=midpoint, y=value), alpha=0.3)+
# 	geom_line(aes(x=midpoint, y=value_real, colour="Fossil record",
# 		      linetype="Fossil record",
# 		      fill="Fossil record")) +
# 	scale_colour_manual("", values=pallete[2:3]) +
# 	scale_linetype("")+
# 	scale_fill_manual("", values=pallete[2:3]) +
# 	scale_x_reverse() +
# 	geom_text(data=plot_labels, aes(x=midpoint, y = 3.7, label=interval),
# 		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+
# 	facet_grid(metric~tetrapod_group, scales = "free_y",
# 		   labeller = labeller(metric = as_labeller(c("alpha_diversity" = "Alpha diversity",
# 		   					   "beta_diversity" = "Beta diversity",
# 		   					   "richness" = "Species richness")))) +
# 	geom_vline(linetype="dashed", colour="black", aes(xintercept=307))
# p_all

# Plot the best-fitting parameters over time
#
# p <- top_gofs_variable %>%
# 	gather(key="parameter", value="value", sigma, speciation_rate, deme) %>%
# 	mutate(parameter=factor(parameter, levels=c("speciation_rate", "deme", "sigma"))) %>%
# 	full_join(interval_data) %>%
# 	ggplot()+
# 	geom_rect(data=interval_data,
# 		  aes(x=NULL, y=NULL,
# 		      xmin=min_ma, xmax=max_ma, ymin=0, ymax=Inf,
# 		      group=interval), fill=rep(c("grey60", "grey80"), 12), colour=NA, alpha=0.5)+
# 	geom_text(data=plot_labels, aes(x=midpoint, y = 1.2e-4, label=interval_part_abbr),
# 		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+
# 	geom_point(aes(x=midpoint, y=value), alpha=0.6)+
# 	geom_smooth(aes(x=midpoint, y=value), method="lm") +
# 	facet_grid(parameter~., scales="free_y", labeller = labeller(parameter=as_labeller(c("deme"="Density",
# 											 "sigma"="Dispersal", "speciation_rate"="Speciation rate"))))+
# 	theme_classic()+
# 	scale_y_log10()+
# 	scale_y_continuous("Parameter value")+
# 	scale_x_reverse("Time (Ma)")
# p

all_gofs_variable <-sim_results %>%
	mutate(scaled_sigma = sigma*sqrt(deme)) %>%
	group_by(scaled_sigma, speciation_rate, tetrapod_group, interval) %>%
	summarise(total_gof= mean(gof)) %>% ungroup() %>%
	mutate(speciation_rate = log10(speciation_rate),
	       scaled_sigma=log10(scaled_sigma))
interp_gofs <- data.frame()
for(interval_sel in unique(all_gofs_variable$interval)){
	for(tetrapod_group_sel in unique(all_gofs_variable$tetrapod_group))
	{
	all_gofs_variable_subset <- all_gofs_variable %>% filter(interval == interval_sel, tetrapod_group==tetrapod_group_sel)
	subset <-  interp(x=all_gofs_variable_subset$scaled_sigma, y=all_gofs_variable_subset$speciation_rate,
			  z=all_gofs_variable_subset$total_gof)
	df <- as.data.frame(subset) %>% select(-x, -y)
	maxnum <- ncol(df)
	names(df) <- 1:maxnum

	df <- df %>% mutate(scaled_sigma=seq(min(all_gofs_variable$scaled_sigma), max(all_gofs_variable$scaled_sigma), length.out = 40)) %>%
		gather(key="speciation_rate", value="gof", -scaled_sigma) %>%
		mutate(speciation_rate=(as.numeric(speciation_rate)*(max(all_gofs_variable$speciation_rate) -
								     	min(all_gofs_variable$speciation_rate))/40) +
		       	min(all_gofs_variable$speciation_rate),
		       interval=interval_sel, tetrapod_group=tetrapod_group_sel)
	# browser()
	interp_gofs <- interp_gofs %>% bind_rows(df)
	}
}

interp_gofs <- interp_gofs %>% na.omit() %>%rowwise() %>%
	mutate(gof=min(gof, 1),
	       interval=factor(interval, levels=levels(interval_data$interval))) %>% na.omit()%>%
	mutate(scaled_sigma=10^scaled_sigma, speciation_rate=10^speciation_rate)
p <- interp_gofs %>%
	ggplot()+
	geom_tile(aes(x=speciation_rate, y=scaled_sigma, fill=gof))+
	scale_x_continuous("Speciation rate", trans="log10",
			   labels=scales::trans_format("log10", scales::math_format()),
			   breaks=c(10^-8, 10^-6, 10^-4))+
	scale_y_log10("Rescaled dispersal")+
	facet_grid(interval~tetrapod_group)+
	theme_classic()+
	scale_fill_viridis("Goodness\nof fit", option="inferno", end=0.95)+
	theme(aspect.ratio=1.0,
	      panel.spacing.x = unit(0.5, "cm"))
pdf(file.path(fig_dir2, "paleo_parameter_space_occ.pdf"), 4.5, 8.5)
print(p)
dev.off()

