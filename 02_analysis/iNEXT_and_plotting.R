# ______________________________________________________
#
#                           iNEXT
#
#                Hsieh, Ma, & Chao (2016) MEE
#
# ______________________________________________________


# # install iNEXT package from CRAN
# install.packages("iNEXT")

# # install the latest version from github
# install.packages('devtools')
# library(devtools)
# install_github('JohnsonHsieh/iNEXT')

library(iNEXT)
library(tidyverse)
library(pbmcapply)
library(viridis)

# The directory containing the fragmented best-fitting results
results_dir <- "results"
data_dir <- "input"
rdata_dir <- "rdata"
figure_dir <- file.path("figures", "inext")
for(i in c(rdata_dir, fig_dir))
{
	if(!dir.exists(i))
	{
		dir.create(i)
	}
}
# Read the interval data
intervals_df <- read.csv(file.path(data_dir, "interval_data.csv")) %>% select(-X)

# Read the best-fitting simulation result
best_df_separate <- read.csv(file.path(results_dir, "single_param_best_sad.csv")) %>%
	select(-X, -sim_path) %>%
	filter(abundance>0, species_id != 0) %>%
	rename(tetrapod_group = tet_group) %>%
	mutate(interval=tools::toTitleCase(as.character(interval)),
	       tetrapod_group = tools::toTitleCase(as.character(tetrapod_group))) %>%
	select(tetrapod_group, abundance, interval) %>%
	group_by(tetrapod_group, interval) %>%
	summarise(total= n(),
		  abundances=list(abundance)) %>%
	rowwise() %>%
	mutate(abundances = list(prepend(abundances, total))) %>%
	select(-total) %>%
	spread(tetrapod_group, abundances)
amniote_abundances <- best_df_separate$Amniote
amphibian_abundances <- best_df_separate$Amphibian
names(amphibian_abundances) <- best_df_separate$interval
names(amniote_abundances) <- best_df_separate$interval

# Generate a single abundances vector combining the amniotes and amphibians
best_df_combined <- read.csv(file.path(results_dir, "single_param_best_sad.csv")) %>%
	select(-X, -sim_path) %>%
	filter(abundance>0, species_id != 0) %>%
	rename(tetrapod_group = tet_group) %>%
	mutate(interval=tools::toTitleCase(as.character(interval)),
	       tetrapod_group = tools::toTitleCase(as.character(tetrapod_group))) %>%
	select(abundance, interval) %>%
	group_by(interval) %>%
	summarise(total= n(),
		  abundances=list(abundance)) %>%
	rowwise() %>%
	mutate(abundances = list(prepend(abundances, total))) %>%
	select(-total)
all_abundances <- best_df_combined$abundances
names(all_abundances) <- best_df_combined$interval

#=====< estimateD >======================================================================================================================================================================================================
quorum_levels <- round(seq(from = 0.1, to = 0.9, by = 0.2), 1) #quorum levels


if(file.exists(file.path(rdata_dir, "estD.RData")))
{
	load(file.path(rdata_dir, "estD.RData"))
}else
{
	# First separate into the two tetrapod groups
	estD_amniote <- pbmclapply(1:length(quorum_levels), function(i) {

	  estD <- estimateD(amniote_abundances, datatype="incidence_freq",
	  		  base="coverage", level=quorum_levels[i])

	  # estD <- estD[estD$order == 0, ] #filter to just species richness

	  estD$reference_t <- sapply(amniote_abundances, sum) #tally total occurrences in each bin

	  # estD[which(estD$t >= 2 * estD$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) #no more than twice ref sample size

	  estD$quorum_level <- quorum_levels[i] #new column

	  estD #returns objects (for apply functions)

	}, mc.cores = detectCores(), mc.preschedule = T, mc.cleanup = T) #for parallel processing

	estD_amphibian <- pbmclapply(1:length(quorum_levels), function(i) {

		estD <- estimateD(amphibian_abundances, datatype="incidence_freq",
				  base="coverage", level=quorum_levels[i])

		# estD <- estD[estD$order == 0,figu ] #filter to just species richness

		estD$reference_t <- sapply(amphibian_abundances, sum) #tally total occurrences in each bin

		# estD[which(estD$t >= 2 * estD$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) #no more than twice ref sample size

		estD$quorum_level <- quorum_levels[i] #new column

		estD #returns objects (for apply functions)

	}, mc.cores = detectCores(), mc.preschedule = T, mc.cleanup = T) #for parallel processing
	estD_combined <- pbmclapply(1:length(quorum_levels), function(i) {

		estD <- estimateD(all_abundances, datatype="incidence_freq",
				  base="coverage", level=quorum_levels[i])

		# estD <- estD[estD$order == 0, ] #filter to just species richness

		estD$reference_t <- sapply(all_abundances, sum) #tally total occurrences in each bin

		# estD[which(estD$t >= 2 * estD$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) #no more than twice ref sample size

		estD$quorum_level <- quorum_levels[i] #new column

		estD #returns objects (for apply functions)

	}, mc.cores = detectCores(), mc.preschedule = T, mc.cleanup = T) #for parallel processing

	estD_amniote_df <- bind_rows(estD_amniote) %>% mutate(tetrapod_group="Amniote")
	estD_amphibian_df <- bind_rows(estD_amphibian) %>% mutate(tetrapod_group="Amphibian")
	estD_combined_df <- bind_rows(estD_combined) %>% mutate(tetrapod_group="Combined")

	estD_df <- estD_amniote_df %>% bind_rows(estD_amphibian_df) %>% bind_rows(estD_combined_df) %>%
		rename(interval=site) %>%
		full_join(intervals_df)

	save(estD_df, estD_amniote_df, estD_amphibian_df, estD_combined_df, file=file.path(rdata_dir, "estD.RData"))
}

p1 <- estD_df %>% filter(order==0) %>%
	ggplot(aes(x=midpoint, y = qD, ymin = qD.LCL, ymax = qD.UCL,
		   colour = factor(quorum_level, levels=c(0.9, 0.7, 0.5, 0.3, 0.1)),
		   fill=factor(quorum_level, levels=c(0.9, 0.7, 0.5, 0.3, 0.1))))+
	geom_line() + theme_classic() + facet_grid(.~tetrapod_group) +
	geom_ribbon(alpha=0.5)+
	scale_fill_viridis("Quorum level", option="plasma", discrete=TRUE, end=0.9) +
	scale_colour_viridis("Quorum level", option="plasma", discrete=TRUE, end=0.9)+
	scale_x_reverse("Time (Ma)") +
	scale_y_log10("Coverage rarified richness (log scale)")
pdf(file.path(figure_dir, "crr_intervals.pdf"), 5.5, 3)
print(p1)
dev.off()

#
#
# # ###### TODO remove/finish
# p <- estD_df %>% filter(order==0) %>%
# 	ggplot(aes(x=quorum_level, y=qD, colour=interval, linetype=method)) +
# 	geom_line() + geom_point()+
# 	facet_grid(. ~ tetrapod_group) + theme_classic()
#
# p
# estD0 <- bind_rows(estD) #binds lists of dataframes
#
# estD0 <- estD0 %>% rename(bin = site) %>%
# 	full_join(., intervals, by = "bin") %>% rename(stage_int = bin) #NB: warning = ok!
#
# estD0[which(estD0$stage_int %in% intervals$bin[1:4]), "period_int"] <- "Permian"
# estD0[which(estD0$stage_int %in% intervals$bin[5:11]), "period_int"] <- "Carboniferous"
#
# estD0$stages <- as.factor(estD0$stages)
# estD0$quorum_level <- as.factor(estD0$quorum_level) #to avoid mystifying errors
#
#
# estD0_0.4 <- subset(estD0, quorum_level == 0.4)
# estD0_0.5 <- subset(estD0, quorum_level == 0.5)
# estD0_0.6 <- subset(estD0, quorum_level == 0.6)
# estD0_0.7 <- subset(estD0, quorum_level == 0.7)
#
#
# ## plot currently only shows Bashkirian-Kungurian estimates
#
# est_plot <- ggplot(filter(estD0, quorum_level %in%
# 			  	quorum_levels[4:7]), aes(x = midpoint, y = qD, ymin = qD.LCL, ymax = qD.UCL, colour = quorum_level)) +
#   coord_fixed(ratio = 1/3) +
#   geom_segment(aes(x = 298.9, xend = 298.9, y = 0, yend = Inf), linetype = "longdash", colour = "grey80", size = 0.7) +
#   #geom_errorbar(width = 0.7) +
#   geom_ribbon(data=estD0_0.4, aes(x = midpoint, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "darkorange1", alpha = 0.3) +
#   geom_ribbon(data=estD0_0.5, aes(x = midpoint, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "palegreen2", alpha = 0.5) +
#   geom_ribbon(data=estD0_0.6, aes(x = midpoint, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "dodgerblue2", alpha = 0.3) +
#   geom_ribbon(data=estD0_0.7, aes(x = midpoint, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "purple1", alpha = 0.3) +
#   geom_line(size = 1.1) +
#   geom_point(aes(pch = method), size = 4.5) +
#   scale_shape_manual(values=c(15, 16, 17)) +
#   scale_colour_manual(values = rev(c("purple3", "royalblue3", "chartreuse4", "darkorange3"))) +
#   scale_x_reverse(expand=c(0,0), limits = c(320.2, 275), breaks = c(315.2, 307, 303.7, 298.9, 295.5, 290.1, 279.2)) +
#   scale_y_continuous(trans = "log10", limits = c(30, 300), breaks = c(10, 50, 100, 150, 200), expand=c(0,0)) +
#   labs(x = "", y = "Coverage rarified richness (log scale)") +
#   theme(panel.background = element_blank(),
#         #legend.position="none",
#         #plot.margin = margin(2, 2, 2, 2, "cm"),
#         panel.grid.minor.y = element_line(colour = "grey90"),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_line(colour = "grey90"),
#         panel.grid.major.x = element_line(colour = "grey90"),
#         panel.border = element_rect(colour = "black", fill = NA),
#         axis.text.x = element_text(size=14, angle=0, hjust=0.5),
#         axis.text.y = element_text(size=14),
#         axis.title = element_text(size=12),
#         aspect.ratio=1)
# est_plot
#
# # y axis = coverage rarified richness
#
#
#
#
#
# #=====< plot with interval bands >======================================================================================================================================================================================================
#
#
# ## (Figure 2 in paper)
#
# plotSTT <- ggplot() +
#   geom_rect(aes(xmin = 356, xmax = 346.7, ymin=0, ymax=Inf), alpha = 0.2, fill = "grey80") +
#   geom_rect(aes(xmin = 330.9, xmax = 323.2, ymin=0, ymax=Inf), alpha = 0.2, fill = "grey80") +
#   geom_rect(aes(xmin = 315.2, xmax = 307.0, ymin=0, ymax=Inf), alpha = 0.2, fill = "grey80") +
#   geom_rect(aes(xmin = 303.7, xmax = 298.9, ymin=0, ymax=Inf), alpha = 0.2, fill = "grey80") +
#   geom_rect(aes(xmin = 295.5, xmax = 290.1, ymin=0, ymax=Inf), alpha = 0.2, fill = "grey80") +
#   geom_rect(aes(xmin = 279.3, xmax = 272, ymin=0, ymax=Inf), alpha = 0.2, fill = "grey80") +
#   geom_segment(aes(x = 298.9, xend = 298.9, y = 0, yend = Inf), linetype = "longdash",colour = "grey80", size = 0.8) +
#   geom_line(data = STT, aes(ages, CeP_collections.all), colour = 'darkcyan', size = 1.2, linetype = "dashed")  +
#   geom_point(data = STT, aes(ages, CeP_collections.all), colour = "darkcyan", size = 4, shape = 17) +
#   geom_line(data = STT, aes(ages, CeP_formations.all), colour = 'peru', size = 1.2, linetype = "twodash")  +
#   geom_point(data = STT, aes(ages, CeP_formations.all), colour = "peru", size = 4, shape = 15) +
#   geom_line(data = STT, aes(ages, CeP_grid_cells50.all), colour = 'chartreuse4', size = 1.2, linetype = "dotted")  +
#   geom_point(data = STT, aes(ages, CeP_grid_cells50.all), colour = "chartreuse4", size = 4, shape = 18) +
#   geom_line(data = STT, aes(ages, CeP_species), colour = 'grey40', size = 1.2)  +
#   geom_point(data = STT, aes(ages, CeP_species), colour = "grey40", size = 3) +
#   theme(panel.background = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA),
#         axis.text.x = element_text(size=14, angle=0, hjust=0.5),
#         axis.text.y = element_text(size=14),
#         axis.title = element_text(size=12, face="bold")) +
#   labs(x = "", y = "") +
#   scale_x_reverse(limits = c(356, 272), breaks=c(350, 330, 310, 290, 270), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, 20), expand = c(0, 0))
# plotSTT




