library(tidyverse)
library(ggthemr)
fig_dir <- "figures/upscaled"
fig_dir2 <- "../../../Figures/Paleo"
for(folder in c(fig_dir, fig_dir2))
{
	if(!dir.exists(folder))
	{
		dir.create(folder)
	}
}
result_dir <- "results"
data_dir <- "input/"
results_dir_pristine <- file.path("results", "Sim7")
results_dir_fragmented <- file.path("results", "Sim8")
best_sims <- read.csv(file.path(results_dir_fragmented, "single_param_best_5.csv"))
upscaled_diversity <- read.csv(file.path(result_dir, "upscaled_best_fitting.csv")) %>%
	filter(percent_cover %in% c(20, 100))
equal_interval_diversity <- read.csv(file.path(result_dir, "upscaled_sampled_best_fitting.csv"))
equal_interval_diversity <- equal_interval_diversity %>%
	bind_rows(equal_interval_diversity %>%
		  	filter(tetrapod_group == "Amniote", interval == "Moscovian") %>%
		  	mutate(interval = "Bashkirian"))%>%
	filter(percent_cover %in% c(20, 100))
interval_data <- read.csv(file.path(data_dir, "interval_data.csv")) %>% select(-X) %>%
	mutate(interval=factor(interval, levels=c("Bashkirian", "Moscovian", "Kasimovian",
						  "Gzhelian", "Asselian", "Sakmarian",
						  "Artinskian", "Kungurian")))
real_biodiversity<- read.csv(file.path(data_dir, "interval_biodiversity_metrics.csv")) %>%
	select(-X) %>%
	mutate(tetrapod_group=tools::toTitleCase(as.character(tetrapod_group)))
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
ggthemr("dust")
pallete <- swatch()
pallete[3] <- pallete[6]
pallete <- c(0, "#62bba5", "#2053A2", "#373634",
	     "#EFA86E","#9A8A76", "#F3C57B","#7A6752","#2A91A2","#87F28A","#6EDCEF")
best_sims_reform <- best_sims %>% full_join(interval_data) %>%
	gather(key="metric", value="value", alpha_diversity, beta_diversity, richness) %>%
	na.omit(tetrapod_group)
p <- upscaled_diversity %>% full_join(interval_data) %>%
	gather(key="metric", value="value", alpha_diversity, beta_diversity, richness) %>%
	na.omit(tetrapod_group) %>%
	ggplot(aes(x=midpoint, y=value))+
	geom_rect(data=interval_data,
		  aes(x=NULL, y=NULL,
		      xmin=min_ma, xmax=max_ma, ymin=0, ymax=Inf,
		      group=interval), fill=rep(c("grey60", "grey80"), 24), colour=NA, alpha=0.5)+
	stat_summary(fun.y = mean, geom = "line",
		     aes(linetype="Simulations\n(10x upscaled\n sampling)",
		         colour="Simulations\n(10x upscaled\n sampling)",
		         fill="Simulations\n(10x upscaled\n sampling)",
		         x=midpoint, y=value)) +
	stat_summary(fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     colour=NA,
		     aes(fill="Simulations\n(10x upscaled\n sampling)",
		         linetype="Simulations\n(10x upscaled\n sampling)",
		         colour="Simulations\n(10x upscaled\n sampling)",
		         x=midpoint, y=value), alpha=0.3)+
	stat_summary(data=best_sims_reform,
		     fun.y = mean, geom = "line",
		     aes(linetype="Simulations\n(equal sampling\n to fossil record)",
		         colour="Simulations\n(equal sampling\n to fossil record)",
		         fill="Simulations\n(equal sampling\n to fossil record)",
		         x=midpoint, y=value)) +
	stat_summary(data=best_sims_reform,
		     fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     colour=NA,
		     aes(fill="Simulations\n(equal sampling\n to fossil record)",
		         linetype="Simulations\n(equal sampling\n to fossil record)",
		         colour="Simulations\n(equal sampling\n to fossil record)",
		         x=midpoint, y=value), alpha=0.3)+
	geom_line(data=real_biodiversity %>%
		  	rename(alpha_diversity=alpha,beta_diversity=beta) %>%
		  	inner_join(interval_data) %>%
		  	gather(key="metric", value="value",
		  	       alpha_diversity, beta_diversity, richness),
		  aes(x=midpoint, y=value, colour="Fossil record",
		      linetype="Fossil record",
		      fill="Fossil record"))+
	theme_classic()+
	ylab("Biodiversity value") + xlab("Time (Ma)") +
	facet_grid(metric~tetrapod_group, scales = "free_y",
		   labeller = labeller(metric = as_labeller(c("alpha_diversity" = "Alpha diversity",
		   					   "beta_diversity" = "Beta diversity",
		   					   "richness" = "Species richness"))))+
	scale_colour_manual("", values=pallete[2:4],
			    breaks=c("Fossil record",
			    	 "Simulations\n(equal sampling\n to fossil record)",
			    	 "Simulations\n(10x upscaled\n sampling)")) +
	scale_linetype_manual("",
		       breaks=c("Fossil record",
				   "Simulations\n(equal sampling\n to fossil record)",
				   "Simulations\n(10x upscaled\n sampling)"),
				   values=c("solid", "dotted", "dashed"))+
	scale_fill_manual("", values=pallete[2:4],
			  breaks=c("Fossil record",
			  	 "Simulations\n(equal sampling\n to fossil record)",
			  	 "Simulations\n(10x upscaled\n sampling)")) +
	scale_x_reverse() +
	geom_text(data=plot_labels, aes(x=midpoint, y = 35, label=interval_abbr),
		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+
	geom_vline(linetype="dashed", colour="black", aes(xintercept=307))+
	theme(legend.key.size=unit(2, unit="lines"))
# print(p)
pdf(file.path(fig_dir, "upscaled_biodiversity.pdf"), 6.5, 4)
print(p)
dev.off()

pdf(file.path(fig_dir2, "paleo_upscaled_biodiversity.pdf"), 6.5, 4)
print(p)
dev.off()



p <- equal_interval_diversity %>% full_join(interval_data) %>%
	gather(key="metric", value="value", alpha_diversity,
	       beta_diversity, richness) %>%
	na.omit(tetrapod_group) %>%
	ggplot(aes(x=midpoint, y=value))+
	geom_rect(data=interval_data,
		  aes(x=NULL, y=NULL,
		      xmin=min_ma, xmax=max_ma, ymin=0, ymax=Inf,
		      group=interval), fill=rep(c("grey60", "grey80"), 24), colour=NA, alpha=0.5)+
	stat_summary(fun.y = mean, geom = "line",
		     aes(linetype="Simulations\n(equal sampling\n over time)",
		         colour="Simulations\n(equal sampling\n over time)",
		         fill="Simulations\n(equal sampling\n over time)",
		         x=midpoint, y=value)) +
	stat_summary(fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     colour=NA,
		     aes(fill="Simulations\n(equal sampling\n over time)",
		         linetype="Simulations\n(equal sampling\n over time)",
		         colour="Simulations\n(equal sampling\n over time)",
		         x=midpoint, y=value), alpha=0.3)+
	stat_summary(data=best_sims %>% gather(key="metric", value="value", alpha_diversity,
					       beta_diversity, richness),
		     fun.y = mean, geom = "line",
		     aes(linetype="Simulations\n(equal sampling\n to fossil record)",
		         colour="Simulations\n(equal sampling\n to fossil record)",
		         fill="Simulations\n(equal sampling\n to fossil record)",
		         x=midpoint, y=value)) +
	stat_summary(data=best_sims %>%
		     	gather(key="metric", value="value", alpha_diversity,
		     	       beta_diversity, richness),
		     fun.ymin = min,
		     fun.ymax = max,
		     geom = "ribbon",
		     colour=NA,
		     aes(fill="Simulations\n(equal sampling\n to fossil record)",
		         linetype="Simulations\n(equal sampling\n to fossil record)",
		         colour="Simulations\n(equal sampling\n to fossil record)",
		         x=midpoint, y=value), alpha=0.3)+
	geom_line(data=real_biodiversity %>%
		  	inner_join(interval_data) %>%
		  	rename(alpha_diversity = "alpha",
		  	       beta_diversity = "beta") %>%
		  	gather(key="metric", value="value", alpha_diversity,
		  	       beta_diversity, richness),
		  aes(x=midpoint, y=value, colour="Fossil record",
		      linetype="Fossil record",
		      fill="Fossil record"))+
	theme_classic()+
	ylab("Biodiversity value") + xlab("Time (Ma)") +
	facet_grid(metric~tetrapod_group, scales = "free_y",
		   labeller=labeller(metric=as_labeller(c("richness" = "Species richness",
		   				     "alpha_diversity" = "Alpha diversity",
		   				     "beta_diversity" = "Beta diversity"))))+
	scale_colour_manual("", values=pallete[2:4],
			    breaks=c("Fossil record",
			    	 "Simulations\n(equal sampling\n to fossil record)",
			    	 "Simulations\n(equal sampling\n over time)")) +
	scale_linetype_manual("", breaks=c("Fossil record",
				    "Simulations\n(equal sampling\n to fossil record)",
				    "Simulations\n(equal sampling\n over time)"),
				    values=c("solid", "dotted", "dashed"))+
	scale_fill_manual("", values=pallete[2:4],
			  breaks=c("Fossil record",
			  	 "Simulations\n(equal sampling\n to fossil record)",
			  	 "Simulations\n(equal sampling\n over time)")) +
	scale_x_reverse() +
	geom_text(data=plot_labels, aes(x=midpoint, y = 20, label=interval_abbr),
		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+
	geom_vline(linetype="dashed", colour="black", aes(xintercept=307))+
	theme(legend.key.size=unit(2, unit="lines"))
pdf(file.path(fig_dir, "upscaled_equal_interval_biodiversity.pdf"), 6.5, 4)
print(p)
dev.off()

pdf(file.path(fig_dir2, "paleo_upscaled_biodiversity_equal_interval.pdf"), 6.5, 4)
print(p)
dev.off()

