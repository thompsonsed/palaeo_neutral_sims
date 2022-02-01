# Data analysis to generate the new info_per_pcoord_occ file for generating simulation data.
# This file uses the occurrences-squared metric for determining the number of individuals that might be present in a simulation.

library(tidyverse)
library(ggthemr)
library(RSQLite)

# Csvs are written to ./input
# Databases are written to data_dir, which is set for the whole file
# data_dir <- "/Volumes/Seagate 3TB/Paleo/Data"
data_dir <-
	file.path("/run", "media", "sam", "Media", "Paleo", "Data")
figure_dir <- file.path("figures")
figure_dir2 <- file.path("../../../Figures/Paleo")
for (folder in c(figure_dir, figure_dir2))
{
	if (!dir.exists(folder)) {
		dir.create(folder)
	}
}
# Read the midpoints for the interval data
interval_data <-
	read.csv(file = file.path("input", "interval_data.csv")) %>% select(-X) %>%
	mutate(interval = factor(
		interval,
		levels = c(
			"Bashkirian",
			"Moscovian",
			"Kasimovian",
			"Gzhelian",
			"Asselian",
			"Sakmarian",
			"Artinskian",
			"Kungurian"
		)
	))

# full_amphibians <- read.csv(file=file.path("input", "early_amphibian_occs_cleaned.csv")) %>%
# 	mutate(tetrapod_group = "amphibian")
# full_amniotes <- read.csv(file=file.path("input", "early_amniote_occs_cleaned.csv")) %>%
# 	mutate(tetrapod_group = "amniote")

# tet_dat <- full_amphibians %>% bind_rows(full_amniotes)

tet_dat <-
	read.csv(file = file.path("input", "tetrapod_data_cleaned.csv")) %>%
	full_join(read.csv(file = file.path("input", "species_per_locality.csv")),
		  by = c("collection_no",
		         "tetrapod_group")) %>%
	full_join(interval_data) %>%
	mutate(mid_ma = (max_ma - min_ma) / 2 + min_ma) %>%
	select(
		accepted_name,
		collection_no,
		paleolat_rounded,
		paleolng_rounded,
		mid_ma,
		tetrapod_group
	) %>%
	filter(mid_ma <= max(interval_data$midpoint),
	       mid_ma >= min(interval_data$midpoint)) %>% na.omit(paleolat) %>%
	rowwise() %>%
	mutate(interval = interval_data[interval_data$max_ma > mid_ma &
						interval_data$min_ma < mid_ma, ]$interval) %>%
	mutate(
		coord_lat = round(paleolat_rounded, digits = 2),
		coord_long = round(paleolng_rounded, digits = 2)
	) %>%
	group_by(coord_lat, coord_long) %>%
	mutate(fragment = paste(sort(unique(collection_no)), sep = "-", collapse =
					"")) %>%
	filter(collection_no != 181085)


info_per_collection <- tet_dat %>% ungroup() %>%
	group_by(fragment,
		 coord_lat,
		 coord_long,
		 tetrapod_group,
		 interval,
		 accepted_name) %>%
	summarise(occurrences = n()) %>% ungroup() %>%
	group_by(fragment, coord_lat, coord_long, tetrapod_group, interval) %>%
	summarise(
		individuals_total = sum(occurrences),
		individuals_total_sq = min
		(individuals_total ^ 2, individuals_total * 5)
	)
# Should be 887 individuals in total
#
info_per_pcoord <-
	tet_dat %>% ungroup() %>% group_by(coord_lat, coord_long, tetrapod_group, interval) %>%
	summarise(species_richness = length(unique(accepted_name))) %>%
	full_join(info_per_collection) %>%
	ungroup() %>% group_by(interval, tetrapod_group) %>%
	mutate(max_ind = max(individuals_total_sq),
	       prop_ind = individuals_total_sq / max_ind)
write.csv(info_per_pcoord, file = file.path("input", "info_per_pcoord_occ.csv"))
info_per_interval <-
	tet_dat %>% group_by(interval, tetrapod_group) %>%
	mutate(pcoord = paste(coord_lat, coord_long)) %>%
	summarise(
		species_richness = length(unique(accepted_name)),
		occurrences = length(unique(accepted_name)),
		collections_no = length(unique(collection_no)),
		localities_no = length(unique(pcoord))
	) %>%
	full_join(
		info_per_collection %>%
			group_by(tetrapod_group, interval) %>%
			summarise(total_individuals_sq = sum(individuals_total_sq))
	)
# write.csv(info_per_pcoord, file=file.path("input", "info_per_interval.csv"))

# Create the csv of global biodiversity metrics
# Find the maximum number of individuals for any one site, per interval
max_number_individuals <- info_per_pcoord %>%
	group_by(interval, tetrapod_group) %>%
	summarise(max_no = max(individuals_total))
# Find the landscape richness per interval
landscape_richness <- tet_dat %>%
	group_by(interval, tetrapod_group) %>%
	summarise(richness = length(unique(accepted_name))) %>%
	full_join(
		info_per_collection %>%
			group_by(tetrapod_group, interval) %>%
			summarise(number_individuals = sum(individuals_total_sq))
	)
# Find the alpha diversity per interval (mean of each site's richness)
alpha <-
	tet_dat %>% group_by(coord_lat, coord_long, interval, tetrapod_group) %>%
	summarise(alpha_ind = length(unique(accepted_name))) %>%
	group_by(interval, tetrapod_group) %>%
	summarise(alpha_diversity = mean(alpha_ind))
# Create a csv of these biodiversity metrics per interval
landscape_richness <- landscape_richness %>% full_join(alpha)
landscape_richness$beta <-
	landscape_richness$richness / landscape_richness$alpha_diversity
write.csv(landscape_richness,
	  file = file.path("input", "interval_biodiversity_metrics_occ.csv"))
# Define the colours for the plot
bw_colours <- rep(c("grey60", "grey80"), 5)
ggthemr("flat dark")
colours <- ggthemr::swatch()
colours[3] <- colours[5]
# Define the plot labels
plot_labels <- data.frame(
	metric = factor(
		"species_richness",
		levels = c(
			"species_richness",
			"alpha_diversity",
			"occurrences",
			"collections_no"
		)
	),
	interval = fct_recode(
		interval_data$interval,
		"Ks" = "Kasimovian",
		"Gz" = "Gzhelian",
		"As" = "Asselian"
	),
	midpoint = interval_data$midpoint
)
# Plot face-value data
p <- info_per_interval %>%
	full_join(interval_data) %>% full_join(alpha) %>%
	gather(
		key = "metric",
		value = "value",
		species_richness,
		occurrences,
		collections_no,
		alpha_diversity
	) %>%
	mutate(
		tetrapod_group = tools::toTitleCase(as.character(tetrapod_group)),
		metric = factor(
			metric,
			levels = c(
				"species_richness",
				"alpha_diversity",
				"occurrences",
				"collections_no"
			)
		)
	) %>%
	ggplot() +
	theme_classic() +
	scale_x_reverse("Time (Ma)") +
	scale_y_continuous("Biodiversity value") +
	facet_grid(metric ~ .,
		   labeller = labeller(metric = as_labeller(
		   	c(
		   		"species_richness" = "Species richness",
		   		"alpha_diversity" = "Alpha diversity",
		   		"occurrences" = "No. individuals",
		   		"collections_no" = "No. collections"
		   	)
		   )),
		   scales = "free_y") +
	scale_linetype("Tetrapod group") +
	scale_colour_manual("Tetrapod group", values = colours[2:6]) +

	geom_rect(
		data = interval_data,
		aes(
			x = NULL,
			y = NULL,
			xmin = min_ma,
			xmax = max_ma,
			ymin = 0,
			ymax = Inf,
			group = interval,
			fill = interval
		),
		colour = NA,
		alpha = 0.5
	) +
	geom_text(
		data = plot_labels,
		aes(x = midpoint, y = 120, label = interval),
		colour = "black",
		size = 2,
		angle = 0,
		vjust = 1,
		hjust = 0.5
	) +
	geom_line(aes(
		x = midpoint,
		y = value,
		linetype = tetrapod_group,
		colour = tetrapod_group
	)) +
	geom_vline(aes(xintercept = 307), colour = "grey10", linetype = "dotted") +
	scale_fill_manual("", values = bw_colours[2:9], guide = FALSE)
pdf(file.path(figure_dir, "raw_diversity_2.pdf"), 6.5, 5)
print(p)
dev.off()
# Save a copy to thesis dir as well
pdf(file.path(figure_dir2, "paleo_raw_diversity.pdf"), 6.5, 5)
print(p)
dev.off()
# Generate the SQLite databases for comparing simulation outputs to
intervals <- unique(landscape_richness$interval)
tetrapod_groups <- unique(landscape_richness$tetrapod_group)
if (!dir.exists(file.path(data_dir, "databases")))
{
	dir.create(file.path(data_dir, "databases"))
}
for (each in intervals)
{
	for (tet_group in tetrapod_groups) {
		interval_df <-
			info_per_pcoord %>% filter(interval == each, tetrapod_group == tet_group)
		## TODO maybe remove this
		fragment_abundances <- tet_dat %>%
			filter(interval == each,
			       tetrapod_group == tet_group) %>%
			group_by(fragment, accepted_name) %>%
			summarise(no_individuals = n() ^ 2) %>%
			mutate(community_reference = 0) %>%
			full_join(
				interval_df %>% select(coord_lat, coord_long, fragment) %>% distinct()
			) %>% ungroup() %>%
			select(accepted_name,
			       no_individuals,
			       community_reference,
			       fragment) %>%
			rename(species_id = accepted_name)
		fragment_richness <- interval_df %>% ungroup() %>%
			select(fragment,
			       species_richness,
			       individuals_total_sq) %>%
			rename(richness = species_richness,
			       no_individuals = individuals_total_sq) %>%
			mutate(community_reference = 0)

		fragment_richness_out <-
			fragment_richness %>% select(fragment, richness) %>%
			mutate(community_reference = 0)
		species_abundances <-
			fragment_abundances %>% group_by(species_id) %>%
			summarise(no_individuals = sum(no_individuals)) %>%
			mutate(community_reference = 0)

		plot_data <- info_per_pcoord %>% filter(interval == interval,
							tetrapod_group == tet_group) %>% ungroup() %>%
			select(fragment, individuals_total) %>%
			rename(no_individuals = individuals_total)
		db_name = paste(
			file.path(data_dir, "databases/"),
			tolower(gsub("\ ", "_", each)),
			"_",
			tet_group,
			"_occ_sq.db",
			sep = ""
		)
		if (file.exists(db_name)) {
			file.remove(db_name)
		}
		db = dbConnect(SQLite(), dbname = db_name)

		dbWriteTable(db,
			     "FRAGMENT_RICHNESS",
			     fragment_richness_out)
		dbWriteTable(db,
			     "FRAGMENT_ABUNDANCES",
			     fragment_abundances)
		dbWriteTable(db, "SPECIES_ABUNDANCES", species_abundances)
		dbWriteTable(db, "PLOT_DATA", plot_data)
		fragment_richness <- fragment_richness %>%
			mutate(metric = "fragment_richness") %>%
			rename(value = richness)
		species_richness <-
			length(unique((
				tet_dat %>% filter(interval == each, tetrapod_group == tet_group)
			)$accepted_name))
		total_individuals <- sum(interval_df$individuals_total_sq)
		biodiversity_metrics <-
			fragment_richness %>% select(-community_reference)
		biodiversity_metrics$metric <- "fragment_richness"
		spec_richness <- data.frame(species_richness,
					    total_individuals)
		spec_richness$fragment <- "whole"
		spec_richness$metric <- "fragment_richness"
		spec_richness <-
			spec_richness %>% rename(value = species_richness, no_individuals = total_individuals)
		names(spec_richness) <-
			c("value",
			  "no_individuals",
			  "fragment",
			  "metric")
		biodiversity_metrics <- rbind(biodiversity_metrics, spec_richness)
		this_landscape <- landscape_richness %>% filter(interval == each,
								tetrapod_group == tet_group)
		alpha <-
			data.frame(
				"alpha_diversity",
				spec_richness$no_individuals[1],
				this_landscape$alpha_diversity,
				"whole"
			)
		beta <-
			data.frame(
				"beta_diversity",
				spec_richness$no_individuals[1],
				this_landscape$beta,
				"whole"
			)
		names(alpha) <- c("metric",
				  "no_individuals",
				  "value",
				  "fragment")
		names(beta) <- names(alpha)
		biodiversity_metrics <- rbind(biodiversity_metrics, alpha, beta)
		dbWriteTable(db,
			     "BIODIVERSITY_METRICS",
			     biodiversity_metrics)
		dbDisconnect(db)
	}
}

# Now calculate the max density per coordinate

max_density_per_pcoord <- info_per_pcoord %>% ungroup() %>%
	group_by(tetrapod_group, interval) %>%
	summarise(max_number_individuals = max(individuals_total_sq)) %>% ungroup() %>%
	mutate_if(is.factor, as.character)
