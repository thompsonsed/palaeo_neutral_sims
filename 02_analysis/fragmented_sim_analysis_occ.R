library(tidyverse)
library(viridis)
library(ggpubr)
library(ggthemr)

# library(ggthemr)
# Set the file paths for input and output
data_dir <- "input/"
results_dir <- file.path("results", "PaleoMainOcc2")
fig_dir <- file.path("figures", "fragmented")
rdata_dir <- "rdata"
for (folder in c(fig_dir))
{
  if (!dir.exists(folder)) {
    dir.create(folder)
  }
}


# Import the real biodiversity data for each interval
real_biodiversity <- read.csv(file.path(data_dir, "interval_biodiversity_metrics_occ.csv")) %>%
  select(-X) %>%
  mutate(tetrapod_group = tools::toTitleCase(as.character(tetrapod_group)))

# Read the interval data
interval_data <- read.csv(file.path(data_dir, "interval_data.csv")) %>%
	select(-X) %>%
	mutate(interval = factor(interval, levels = c(
		"Bashkirian", "Moscovian", "Kasimovian",
		"Gzhelian", "Asselian", "Sakmarian",
		"Artinskian", "Kungurian"
	)))

# Import the coalescence simulation results


sim_results_all <- read.csv(file.path(results_dir, "results_all_occ_sq.csv")) %>%
  mutate(deme = floor(deme)) %>%
  mutate(
    interval = tools::toTitleCase(as.character(interval)),
    tetrapod_group = tools::toTitleCase(as.character(tetrapod_group))
  ) %>%
  full_join(interval_data, by = "interval") %>%
  full_join(real_biodiversity %>% select(interval, tetrapod_group, number_individuals)) %>%
	filter(!(midpoint >= 305 & scenario == "fragmented"))


sim_results_pristine <- sim_results_all %>% filter(scenario == "pristine")
sim_results_fragmented <- sim_results_all %>%
  filter(scenario == "fragmented") %>%
  bind_rows(sim_results_pristine %>% mutate(percent_cover = 0.2)) %>%
  bind_rows(sim_results_pristine %>% mutate(percent_cover = 0.4)) %>%
  bind_rows(sim_results_pristine %>% mutate(percent_cover = 0.8)) %>%
  filter(
    (midpoint >= 305 & scenario == "pristine") |
      (midpoint < 305 & scenario == "fragmented"), # this is 305 so that only times which are entirely post-CRC are included
    # simulated_individuals < 5*number_individuals) %>%
  ) %>%
  mutate(scenario = "fragmented")
sim_results_clustered <- sim_results_all %>%
  filter(scenario == "clustered") %>%
  bind_rows(sim_results_pristine) %>%
  filter(
    (midpoint >= 305 & scenario == "pristine") |
      (midpoint < 305 & scenario == "clustered"), # this is 305 so that only times which are entirely post-CRC are included
    # simulated_individuals < 5*number_individuals) %>%
  ) %>%
  mutate(scenario = "clustered")


counts_all <- sim_results_all %>%
  group_by(interval, tetrapod_group, scenario) %>%
  summarise(total = n()) %>%
  left_join(interval_data)


sim_results <- sim_results_pristine %>%
  bind_rows(sim_results_clustered) %>%
  bind_rows(sim_results_fragmented) %>%
  mutate(deme = round(deme)) %>%
  group_by(scenario, tetrapod_group, sigma, deme, speciation_rate, percent_cover) %>%
  mutate(number_sims = n()) %>%
  filter(number_sims >= 7) %>%
  ungroup()


# Find the number of complete sims for each parameter set
complete_sims <- sim_results %>%
  ungroup() %>%
  select(scenario, tetrapod_group, sigma, deme, percent_cover, number_sims) %>%
  distinct()

# Modify the real biodiversity data
real_biodiversity <- real_biodiversity %>%
  full_join(interval_data, by = "interval") %>%
  mutate(
    tetrapod_group = tools::toTitleCase(as.character(tetrapod_group)),
    interval = tools::toTitleCase(as.character(interval))
  )



# Add the bashkirian (which only has a single individual in the fossil record)
deme_and_sigma <- sim_results %>% select(sigma, deme) %>% distinct() %>% mutate(t=row_number())
bashkirian_result <- expand.grid(
	actual_alpha = 1, actual_beta = 1, actual_richness = 1,
	alpha_diversity = 1, beta_diversity = 1, t=seq(length(deme_and_sigma$sigma)),
	gof = 1.0, interval = "Bashkirian", percent_cover = c(0.2, 0.4, 0.8),
	richness = 1, scenario = c("fragmented"), simulated_individuals = 1,
	speciation_rate = unique(sim_results$speciation_rate), tetrapod_group = "Amniote",
	min_ma = 315.2, max_ma = 323.2, midpoint = 319.20, map_file = "320m", number_individuals = 1, number_sims = 240, total_gof = 1,
	top_3 = TRUE
) %>% full_join(deme_and_sigma) %>% select(-t)

sim_results <- sim_results %>% bind_rows(bashkirian_result)

write.csv(sim_results, file.path(results_dir, "results_fragmented_finalised_occ.csv"))

# Count results by parameters

parameter_sets <- sim_results %>% ungroup() %>% group_by(tetrapod_group, sigma, deme, scenario) %>%
	summarise(total=n())

## METHOD B
# Look at the top 5 fitting results for each tetrapod group
single_param_result_separate <- sim_results %>%
  ungroup() %>%
  group_by(scenario, sigma, speciation_rate, deme, tetrapod_group, percent_cover) %>%
  mutate(total_gof = mean(gof))

# Compute some stats about the parameters for the top 3 fitting simulations
summary_stats_separate <- single_param_result_separate %>%
  select(scenario, sigma, speciation_rate, deme, total_gof, tetrapod_group, percent_cover) %>%
  distinct() %>%
  ungroup() %>%
  group_by(tetrapod_group) %>%
  top_n(total_gof, n = 15) %>%
  gather(key = "parameter", value = "value", sigma, speciation_rate, deme, total_gof) %>%
  group_by(scenario, parameter, tetrapod_group) %>%
  summarise(
    mean_value = mean(value),
    min_value = min(value),
    max_value = max(value),
    sd_value = sd(value),
    total = n()
  )
top_gofs <- single_param_result_separate %>%
  ungroup() %>%
  select(scenario, total_gof, tetrapod_group, percent_cover) %>%
  distinct() %>%
  group_by(scenario, tetrapod_group, percent_cover) %>%
  top_n(n = 3, wt = total_gof) %>%
  summarise(
    min = min(total_gof),
    max = max(total_gof)
  )

# Get the top 5 best-fitting results
single_param_result_separate <- single_param_result_separate %>%
  ungroup() %>%
  group_by(scenario, tetrapod_group, percent_cover) %>%
  rowwise() %>%
  mutate(top_3 = total_gof >= top_gofs[(top_gofs$tetrapod_group == tetrapod_group & top_gofs$scenario == scenario), ]$min[1]) %>%
  filter(top_3)

# DECIDE WHICH METHOD
# Comment this out if you want to use the same parameters for amphibians and amniotes
# So with this line, use METHOD B for analysis
single_param_result <- single_param_result_separate



single_param_result <- single_param_result %>% bind_rows(bashkirian_result)
# Plot the landscape richness over time, compared to the real species richness over time
# p <- ggplot(sim_results, aes(x=midpoint, y=richness, group=1))  +
# 	theme_classic() + ylab("Species richness") + xlab("Interval") +
# 	stat_summary(fun.y = mean, geom = "line", aes(linetype="Simulated",
# 						      colour=tetrapod_group,
# 						      group=tetrapod_group)) +
# 	stat_summary(fun.data = mean_se, geom = "errorbar")+geom_point()+
# 	geom_line(data=real_biodiversity, aes(x=midpoint, colour=tetrapod_group,
# 					      y=richness, group=tetrapod_group,linetype="Fossil record")) +
# 	scale_linetype_discrete("Data source") + scale_x_reverse()
# p
# Combine all three into a single plot
# Define the plot labels
plot_labels <- data.frame(
  metric = factor("alpha_diversity", levels = c(
    "richness",
    "alpha_diversity",
    "beta_diversity"
  )),
  parameter = "speciation_rate",
  interval = interval_data$interval,
  interval_abbr = fct_recode(interval_data$interval,
    "Ks" = "Kasimovian", "Gz" = "Gzhelian",
    "As" = "Asselian", "Ba" = "Bashkirian", "Mo" = "Moscovian", "Sa" = "Sakmarian",
    "Ar" = "Artinskian", "Ku" = "Kungurian"
  ),
  midpoint = interval_data$midpoint,
  interval_part_abbr = fct_recode(interval_data$interval, "Ks" = "Kasimovian", "Gz" = "Gzhelian", "As" = "Asselian")
)
p_all <- single_param_result %>%
  filter(scenario == "fragmented") %>%
  gather(
    key = "metric", value = "value",
    beta_diversity, alpha_diversity, richness
  ) %>%
  mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
    ifelse(metric == "beta_diversity",
      actual_beta, actual_richness
    )
  )) %>%
  ggplot() + theme_classic() +
  geom_rect(
    data = interval_data,
    aes(
      x = NULL, y = NULL,
      xmin = min_ma, xmax = max_ma, ymin = 0, ymax = Inf,
      group = interval
    ), fill = rep(c("grey60", "grey80"), 24), colour = NA, alpha = 0.5
  ) + geom_text(
    data = plot_labels, aes(x = midpoint, y = 4.2, label = interval_abbr),
    colour = "black", size = 2, angle = 0, vjust = 1, hjust = 0.5
  ) +

  ylab("Biodiversity value") + xlab("Time (Ma)") +
  stat_summary(fun.y = mean, geom = "line", aes(
    linetype = "Simulated",
    colour = as.factor(percent_cover),
    group = as.factor(percent_cover),
    x = midpoint, y = value
  )) +
  stat_summary(
    fun.ymin = min,
    fun.ymax = max,
    geom = "ribbon",
    colour = NA,
    aes(
      fill = as.factor(percent_cover),
      group = as.factor(percent_cover),
      x = midpoint, y = value
    ), alpha = 0.3
  ) +

  geom_line(aes(x = midpoint, y = value_real, linetype = "Fossil record")) +
  scale_colour_viridis("Remaining\nhabitat %", discrete = TRUE, begin = 0.2, end = 0.65, option = "plasma",
  		     breaks=c(1.0, 0.8, 0.4, 0.2),
  		     labels=c("100","80","40","20")) +
  scale_fill_viridis("Remaining\nhabitat %", discrete = TRUE, begin = 0.2, end = 0.65, option = "plasma",
  		   breaks=c(1.0, 0.8, 0.4, 0.2),
  		   labels=c("100","80","40","20")) +
  scale_linetype_discrete("") + scale_x_reverse() +
  facet_grid(metric ~ tetrapod_group,
    scales = "free_y",
    labeller = labeller(metric = as_labeller(c(
      "alpha_diversity" = "Alpha diversity",
      "beta_diversity" = "Beta diversity",
      "richness" = "Species richness"
    )))
  ) +
	geom_vline(linetype = "dashed", colour = "black", aes(xintercept = 307))+

  guides(linetype = guide_legend(override.aes = list(colour = c("#62bba5", "black"))))
print(p_all)
pdf(file.path(fig_dir, "fragmented_best_all_metrics_all_cover_occ.pdf"), 6.5, 4)
print(p_all)
dev.off()


# Plot the number of individuals
p_number_individuals <- single_param_result %>%
  filter(scenario == "fragmented") %>%
  gather(
    key = "metric", value = "value",
    simulated_individuals, number_individuals
  ) %>%
  # mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
  # 			   ifelse(metric == "beta_diversity",
  # 			          actual_beta, actual_richness))) %>%
  ggplot() + theme_classic() +
  geom_rect(
    data = interval_data,
    aes(
      x = NULL, y = NULL,
      xmin = min_ma, xmax = max_ma, ymin = 0, ymax = Inf,
      group = interval
    ), fill = rep(c("grey60", "grey80"), 4), colour = NA, alpha = 0.5
  ) +
  geom_text(
    data = plot_labels, aes(x = midpoint, y = 3.7, label = interval_abbr),
    colour = "black", size = 2, angle = 0, vjust = 1, hjust = 0.5
  ) +
  geom_line(aes(x = midpoint, y = value, linetype = metric, colour = tetrapod_group))

# p_number_individuals
# # Find the single one best parameter set across sims for output
# write.csv(single_param_result, file.path(results_dir_fragmented, "single_param_best_5_occ.csv"))
ggthemr_reset()
ggthemr("light")
pallete <- swatch()
pallete[3] <- pallete[9]
p_all <- single_param_result %>%
  filter(percent_cover == 0.2) %>%
  gather(
    key = "metric", value = "value",
    beta_diversity, alpha_diversity, richness
  ) %>%
  mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
    ifelse(metric == "beta_diversity",
      actual_beta, actual_richness
    )
  )) %>%
  ggplot() +
  theme_classic() +
  ylab("Biodiversity value") + xlab("Time (Ma)") +
  geom_rect(
    data = interval_data,
    aes(
      x = NULL, y = NULL,
      xmin = min_ma, xmax = max_ma, ymin = 0, ymax = Inf,
      group = interval
    ), fill = rep(c("grey60", "grey80"), 24), colour = NA, alpha = 0.5
  ) +
  stat_summary(
    fun.y = mean, geom = "line",
    aes(
      linetype = "Simulations\n(with CRC)",
      colour = "Simulations\n(with CRC)",
      fill = "Simulations\n(with CRC)",
      x = midpoint, y = value
    )
  ) +
  stat_summary(
    fun.ymin = min,
    fun.ymax = max,
    geom = "ribbon",
    colour = NA,
    aes(
      fill = "Simulations\n(with CRC)",
      linetype = "Simulations\n(with CRC)",
      colour = "Simulations\n(with CRC)",
      x = midpoint, y = value
    ), alpha = 0.3
  ) +
  geom_line(aes(
    x = midpoint, y = value_real, colour = "Fossil record",
    linetype = "Fossil record",
    fill = "Fossil record"
  )) +
  scale_colour_manual("", values = pallete[2:3]) +
  scale_linetype("") +
  scale_fill_manual("", values = pallete[2:3]) +
  scale_x_reverse() +
  geom_text(
    data = plot_labels, aes(x = midpoint, y = 4.2, label = interval_abbr),
    colour = "black", size = 2, angle = 0, vjust = 1, hjust = 0.5
  ) +
  facet_grid(metric ~ tetrapod_group,
    scales = "free_y",
    labeller = labeller(metric = as_labeller(c(
      "alpha_diversity" = "Alpha diversity",
      "beta_diversity" = "Beta diversity",
      "richness" = "Species richness"
    )))
  ) +
  geom_vline(linetype = "dashed", colour = "black", aes(xintercept = 307))
print(p_all)
pdf(file.path(fig_dir, "fragmented_best_all_metrics_20_cover_occ.pdf"), 6.5, 4)
print(p_all)
dev.off()

p_all_pristine <- single_param_result %>%
  filter(scenario == "pristine") %>%
  gather(
    key = "metric", value = "value",
    beta_diversity, alpha_diversity, richness
  ) %>%
  mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
    ifelse(metric == "beta_diversity",
      actual_beta, actual_richness
    )
  )) %>%
  ggplot() +
  theme_classic() +
  ylab("Biodiversity value") + xlab("Time (Ma)") +
  geom_rect(
    data = interval_data,
    aes(
      x = NULL, y = NULL,
      xmin = min_ma, xmax = max_ma, ymin = 0, ymax = Inf,
      group = interval
    ), fill = rep(c("grey60", "grey80"), 24), colour = NA, alpha = 0.5
  ) +
  stat_summary(
    fun.y = mean, geom = "line",
    aes(
      linetype = "Simulations\n(without CRC)",
      colour = "Simulations\n(without CRC)",
      fill = "Simulations\n(without CRC)",
      x = midpoint, y = value
    )
  ) +
  stat_summary(
    fun.ymin = min,
    fun.ymax = max,
    geom = "ribbon",
    colour = NA,
    aes(
      fill = "Simulations\n(without CRC)",
      linetype = "Simulations\n(without CRC)",
      colour = "Simulations\n(without CRC)",
      x = midpoint, y = value
    ), alpha = 0.3
  ) +
  geom_line(aes(
    x = midpoint, y = value_real, colour = "Fossil record",
    linetype = "Fossil record",
    fill = "Fossil record"
  )) +
  scale_colour_manual("", values = pallete[2:3]) +
  scale_linetype("") +
  scale_fill_manual("", values = pallete[2:3]) +
  scale_x_reverse() +
  geom_text(
    data = plot_labels, aes(x = midpoint, y = 4.2, label = interval_abbr),
    colour = "black", size = 2, angle = 0, vjust = 1, hjust = 0.5
  ) +
  facet_grid(metric ~ tetrapod_group,
    scales = "free_y",
    labeller = labeller(metric = as_labeller(c(
      "alpha_diversity" = "Alpha diversity",
      "beta_diversity" = "Beta diversity",
      "richness" = "Species richness"
    )))
  ) +
  geom_vline(linetype = "dashed", colour = "black", aes(xintercept = 307))
print(p_all_pristine)
pdf(file.path(fig_dir, "pristine_best_all_metrics_occ.pdf"), 6.5, 4)
print(p_all_pristine)
dev.off()

ggthemr_reset()
ggthemr("light")
pallete <- swatch()
pallete[3] <- pallete[7]

p_all_clustered <- single_param_result %>%
	filter(scenario == "clustered") %>%
	gather(
		key = "metric", value = "value",
		beta_diversity, alpha_diversity, richness
	) %>%
	mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
				   ifelse(metric == "beta_diversity",
				          actual_beta, actual_richness
				   )
	)) %>%
	ggplot() +
	theme_classic() +
	ylab("Biodiversity value") + xlab("Time (Ma)") +
	geom_rect(
		data = interval_data,
		aes(
			x = NULL, y = NULL,
			xmin = min_ma, xmax = max_ma, ymin = 0, ymax = Inf,
			group = interval
		), fill = rep(c("grey60", "grey80"), 24), colour = NA, alpha = 0.5
	) +
	stat_summary(
		fun.y = mean, geom = "line",
		aes(
			linetype = "Simulations\n(with CRC - clusters)",
			colour = "Simulations\n(with CRC - clusters)",
			fill = "Simulations\n(with CRC - clusters)",
			x = midpoint, y = value
		)
	) +
	stat_summary(
		fun.ymin = min,
		fun.ymax = max,
		geom = "ribbon",
		colour = NA,
		aes(
			fill = "Simulations\n(with CRC - clusters)",
			linetype = "Simulations\n(with CRC - clusters)",
			colour = "Simulations\n(with CRC - clusters)",
			x = midpoint, y = value
		), alpha = 0.3
	) +
	geom_line(aes(
		x = midpoint, y = value_real, colour = "Fossil record",
		linetype = "Fossil record",
		fill = "Fossil record"
	)) +
	scale_colour_manual("", values = pallete[2:3]) +
	scale_linetype("") +
	scale_fill_manual("", values = pallete[2:3]) +
	scale_x_reverse() +
	geom_text(
		data = plot_labels, aes(x = midpoint, y = 4.2, label = interval_abbr),
		colour = "black", size = 2, angle = 0, vjust = 1, hjust = 0.5
	) +
	facet_grid(metric ~ tetrapod_group,
		   scales = "free_y",
		   labeller = labeller(metric = as_labeller(c(
		   	"alpha_diversity" = "Alpha diversity",
		   	"beta_diversity" = "Beta diversity",
		   	"richness" = "Species richness"
		   )))
	) +
	geom_vline(linetype = "dashed", colour = "black", aes(xintercept = 307))
pdf(file.path(fig_dir, "clustered_best_all_metrics_occ.pdf"), 6.5, 4)
print(p_all_clustered)
dev.off()


# Project early results later into time


# Analyse just the late-stage parameters for best-fitting simulations

time_split_results <- sim_results %>%
	filter(scenario == "pristine") %>%
	mutate(post_307 = midpoint < 307) %>%
	group_by(sigma, speciation_rate, deme, tetrapod_group, post_307) %>%
	mutate(total_gof = mean(gof))

# Compute some stats about the parameters for the top 3 fitting simulations
summary_stats_split <- time_split_results  %>%
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

# Calculate the minimum goodness-of-fit in the top 5 fitting sims for each tetrapod group
top_gofs_split <- time_split_results %>% ungroup() %>% select(total_gof, tetrapod_group, post_307)%>%
	distinct() %>% group_by(tetrapod_group, post_307) %>% top_n(n=3, wt=total_gof) %>%
	summarise(min=min(total_gof))
early_min_results <- (top_gofs_split %>% filter(!post_307))
late_min_results <- (top_gofs_split %>% filter(post_307))
early_results <- time_split_results %>% ungroup() %>%
	mutate(total_gof = ifelse(post_307, 0.0, total_gof)) %>%
	group_by(sigma, speciation_rate, deme, tetrapod_group) %>%
	mutate(total_gof = max(total_gof)) %>% ungroup() %>%
	rowwise() %>%
	mutate(top_5 = total_gof >= early_min_results[early_min_results$tetrapod_group == tetrapod_group,]$min[1]) %>%
	filter(top_5)
late_results <- time_split_results %>%
	ungroup() %>%
	mutate(total_gof = ifelse(!post_307, 0.0, total_gof)) %>%
	group_by(sigma, speciation_rate, deme, tetrapod_group) %>%
	mutate(total_gof = max(total_gof)) %>% ungroup() %>%
	rowwise() %>%
	mutate(top_5 = total_gof >= late_min_results[late_min_results$tetrapod_group == tetrapod_group,]$min[1]) %>%
	filter(top_5)
split_results <- time_split_results %>% ungroup() %>%
	rowwise() %>%
	mutate(top_5 = total_gof >= top_gofs_split[top_gofs_split$tetrapod_group == tetrapod_group &
						   	top_gofs_split$post_307 == post_307,]$min[1]) %>%
	filter(top_5, (post_307 & midpoint < 307) | (!post_307 & midpoint > 307))

p <- split_results  %>%
	gather(key="metric", value="value",beta_diversity, alpha_diversity, richness) %>%
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
print(p)
pdf(file.path(fig_dir, "pristine_split_best_occ.pdf"), 6.5, 4)
print(p)
dev.off()

# Plot the early results
ggthemr_reset()
ggthemr("light")
pallete <- swatch()
pallete[3] <- pallete[5]
p <- early_results %>%
	gather(key="metric", value="value",beta_diversity, alpha_diversity, richness) %>%
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
	geom_text(data=plot_labels, aes(x=midpoint, y = 4.6, label=interval_abbr),
		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+
	facet_grid(metric~tetrapod_group, scales = "free_y",
		   labeller = labeller(metric = as_labeller(c("alpha_diversity" = "Alpha diversity",
		   					   "beta_diversity" = "Beta diversity",
		   					   "richness" = "Species richness")))) +
	geom_vline(linetype="dashed", colour="black", aes(xintercept=307))
pdf(file.path(fig_dir, "pristine_early_best_occ.pdf"), 6.5, 4)
print(p)
dev.off()


# Plot the late results

p <- late_results %>%
	gather(key="metric", value="value",beta_diversity, alpha_diversity, richness) %>%
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
pdf(file.path(fig_dir, "pristine_late_best_occ.pdf"), 6.5, 4)
print(p)
dev.off()


# Plot the impact of fragmentation itself, rather than finding the best-fitting model

pristine_fragmented_results <- sim_results %>%
	group_by(sigma, speciation_rate, deme, tetrapod_group, scenario, percent_cover) %>%
	mutate(total_gof = mean(gof)) %>%
	ungroup() %>%
	mutate(total_gof = ifelse(scenario == "pristine", total_gof, 0.0)) %>%
	group_by(sigma, speciation_rate, deme, tetrapod_group) %>%
	mutate(total_gof = max(total_gof))

# Compute some stats about the parameters for the top 3 fitting simulations
summary_stats_pf <- pristine_fragmented_results  %>%
	select(sigma, speciation_rate, deme, total_gof, tetrapod_group, scenario) %>% distinct() %>%
	ungroup() %>% group_by(tetrapod_group) %>%
	top_n(total_gof, n=5) %>%
	gather(key="parameter", value="value", sigma, speciation_rate, deme, total_gof) %>%
	group_by(parameter, tetrapod_group, scenario) %>%
	summarise(mean_value = mean(value),
		  min_value = min(value),
		  max_value = max(value),
		  sd_value = sd(value),
		  total=n())

# Calculate the minimum goodness-of-fit in the top 5 fitting sims for each tetrapod group
top_gofs_pf <- pristine_fragmented_results %>% ungroup() %>%
	select(total_gof, tetrapod_group, scenario, percent_cover) %>%
	distinct() %>%
	# mutate(c = n()) %>%
	# filter(c == 3) %>%
	group_by(tetrapod_group, percent_cover) %>%
	# mutate(number_scenarios = length(unique(scenario))) %>%
	# filter(number_scenarios > 2) %>%
	top_n(n=10, wt=total_gof) %>%
	summarise(min=min(total_gof))


pristine_fragmented_results_plot <- pristine_fragmented_results %>%
	rowwise() %>%
	mutate(top_5 = total_gof >= top_gofs_pf[top_gofs_pf$tetrapod_group == tetrapod_group,]$min[1],
	       use_this = ifelse(percent_cover != 1.0, ifelse(midpoint > 307, FALSE, TRUE), TRUE)) %>%
	filter(top_5, use_this)
pristine_fragmented_results <- pristine_fragmented_results_plot %>% mutate(post307 = midpoint < 307) %>%
	filter(scenario == "fragmented") %>%
	bind_rows(single_param_result %>% filter(scenario == "pristine"))
pristine_fragmented_results <- pristine_fragmented_results %>%
	bind_rows(pristine_fragmented_results %>%
		  	filter(percent_cover == 1.0, midpoint >= 307) %>%
		  	mutate(percent_cover = 0.2)) %>%
	bind_rows(pristine_fragmented_results %>%
		  	filter(percent_cover == 1.0, midpoint >= 307) %>%
		  	mutate(percent_cover = 0.4)) %>%
	bind_rows(pristine_fragmented_results %>%
		  	filter(percent_cover == 1.0, midpoint >= 307) %>%
		  	mutate(percent_cover = 0.8))
ggthemr_reset()
ggthemr("light")
pallete <- swatch()
pallete[3] <- pallete[5]
p_all <- pristine_fragmented_results %>%
	gather(
		key = "metric", value = "value",
		beta_diversity, alpha_diversity, richness
	) %>%
	mutate(value_real = ifelse(metric == "alpha_diversity", actual_alpha,
				   ifelse(metric == "beta_diversity",
				          actual_beta, actual_richness
				   )),
	   percent_cover = as.character(percent_cover)) %>%
	ggplot() + theme_classic() +
	geom_rect(
		data = interval_data,
		aes(
			x = NULL, y = NULL,
			xmin = min_ma, xmax = max_ma, ymin = 0, ymax = Inf,
			group = interval
		), fill = rep(c("grey60", "grey80"), 24), colour = NA, alpha = 0.5
	) + geom_text(
		data = plot_labels, aes(x = midpoint, y = 4.2, label = interval_abbr),
		colour = "black", size = 2, angle = 0, vjust = 1, hjust = 0.5
	) +

	ylab("Biodiversity value") + xlab("Time (Ma)") +
	stat_summary(fun.y = mean, geom = "line", aes(
		linetype = "Simulated",
		colour = as.factor(percent_cover),
		group = as.factor(percent_cover),
		x = midpoint, y = value
	)) +
	stat_summary(
		fun.ymin = min,
		fun.ymax = max,
		geom = "ribbon",
		colour = NA,
		aes(
			fill = as.factor(percent_cover),
			group = as.factor(percent_cover),
			x = midpoint, y = value
		), alpha = 0.3
	) +

	geom_line(aes(x = midpoint, y = value_real, linetype = "Fossil record")) +
	scale_colour_viridis("Remaining\nhabitat %", discrete = TRUE, begin = 0.2, end = 0.65, option = "plasma",
			     breaks=c(1.0, 0.8, 0.4, 0.2),
			     labels=c("100","80","40","20")) +
	scale_fill_viridis("Remaining\nhabitat %", discrete = TRUE, begin = 0.2, end = 0.65, option = "plasma",
			   breaks=c(1.0, 0.8, 0.4, 0.2),
			   labels=c("100","80","40","20")) +
	scale_linetype_discrete("") + scale_x_reverse() +
	facet_grid(metric ~ tetrapod_group,
		   scales = "free_y",
		   labeller = labeller(metric = as_labeller(c(
		   	"alpha_diversity" = "Alpha diversity",
		   	"beta_diversity" = "Beta diversity",
		   	"richness" = "Species richness"
		   )))
	) +
	geom_vline(linetype = "dashed", colour = "black", aes(xintercept = 307))+
	guides(linetype = guide_legend(override.aes = list(colour = c("#62bba5", "black"))))
print(p_all)
pdf(file.path(fig_dir, "fragmented_predictions_best_all_metrics_all_cover_occ.pdf"), 6.5, 4)
print(p_all)
dev.off()



# Check the change in parameters over time

parameter_results <- sim_results %>% filter(gof> 0.8) %>%
	group_by(scenario, interval, tetrapod_group) %>%
	pivot_longer(cols = c(deme, sigma), names_to = "parameter")
parameter_results %>% filter(percent_cover == 1.0) %>%
	ggplot(aes(x=midpoint, y=value))+ theme_classic()+
	stat_summary(fun.y = mean, geom = "line", colour="grey50", aes(

	# colour = as.factor(percent_cover),
	# group = as.factor(percent_cover),
	x = midpoint, y = value
)) +
	geom_point(aes(colour=gof))+
	stat_summary(
		fun.ymin = min,
		fun.ymax = max,
		geom = "ribbon",
		colour=NA,
		fill="grey50",
		aes(
			# fill = as.factor(percent_cover),
			# group = as.factor(percent_cover),
			x = midpoint, y = value
		), alpha = 0.3
	)+
	facet_grid(parameter~., scales = "free_y")+
	geom_vline(linetype = "dashed", colour = "black", aes(xintercept = 307))+
	geom_rect(
		data = interval_data ,
		aes(
			x = NULL, y = NULL,
			xmin = min_ma, xmax = max_ma, ymin = 0, ymax = Inf,
			group = interval
		), fill = rep(c("grey60", "grey80"), 8), colour = NA, alpha = 0.5
	) + geom_text(
		data = plot_labels %>% mutate(parameter = "deme"), aes(x = midpoint, y = 1000, label = interval_abbr),
		colour = "black", size = 2, angle = 0, vjust = 1, hjust = 0.5
	) + scale_x_reverse()

