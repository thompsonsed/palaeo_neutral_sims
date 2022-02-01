library(tidyverse)
library(ggthemr)
library(RSQLite)

# Csvs are written to ./input
# Databases are written to data_dir, which is set for the whole file
# data_dir <- "/Volumes/Seagate 3TB/Paleo/Data"
data_dir <- "../../Data/"
# data_dir <- file.path("/run", "media", "sam", "Media", "Paleo", "Data")
figure_dir <- file.path("figures")
figure_dir2 <- file.path("../../../Figures/Paleo")
for(folder in c(figure_dir, figure_dir2))
{
	if(!dir.exists(folder)){
		dir.create(folder)
	}
}

# Tetrapod data import and cleaning
tetdat <-  read.csv(file=file.path("input", "tetrapod_data_cleaned.csv"))
tetdat <- unite(tetdat, pcoords, paleolat_rounded, paleolng_rounded, sep = ",")
site_biodiversity <- read.csv(file=file.path("input", "species_per_locality.csv"))
interval_data <- read.csv(file=file.path("input", "interval_data.csv"))%>%
	mutate(interval=factor(interval, levels=c("Bashkirian", "Moscovian", "Kasimovian",
						  "Gzhelian", "Asselian", "Sakmarian",
						  "Artinskian", "Kungurian")))

# Combine the data - this is used for generating the comparison databases
site_biodiversity <- site_biodiversity %>% full_join(tetdat, by=c("collection_no",
								  "tetrapod_group"))

# These are the ones with NAs in them
res_nas <- site_biodiversity %>% filter(is.na(species_total)) %>%
	dplyr::select(collection_no, accepted_name, tetrapod_group)
write.csv(res_nas, file.path("input", "missing_data_tetrapods.csv"))
# other_nas <- site_biodiversity %>% filter(is.na(accepted_name))
# write.csv(other_nas, "./input/missing_data_accepted_name.csv")
site_biodiversity <- site_biodiversity %>% na.omit() %>% select(-X) %>%
	group_by(pcoords, interval, tetrapod_group) %>%
	mutate(fragment_name = paste(sort(unique(collection_no)), sep="-", collapse=""))



# Generate the aggregated biodiversity metrics per pcoord and interval - this is used for generation
# of the map files. This tallies counts of individuals and species per site and interval
info_per_pcoord <- site_biodiversity %>% group_by(fragment_name, pcoords, interval, tetrapod_group) %>%
	summarise(collections = length(unique(collection_no)),
		  species_total = length(unique(accepted_name)),
		  individuals_total = n())
info_per_pcoord <- info_per_pcoord %>% ungroup() %>% group_by(interval, tetrapod_group) %>%
	mutate(max_ind = max(individuals_total),
	       prop_ind = individuals_total/max_ind)
write.csv(info_per_pcoord, file=file.path("input", "info_per_pcoord_main.csv"))
info_per_interval <- site_biodiversity %>% group_by(interval, tetrapod_group) %>%
	summarise(species_richness = length(unique(accepted_name)),
		  occurrences = length(unique(accepted_name)),
		  collections_no=length(unique(collection_no)),
		  localities_no=length(unique(fragment_name)))
# Find the abundance of each species in each interval, site and tetrapod group
site_biodiversity <- site_biodiversity %>% left_join(info_per_pcoord) %>%
	group_by(fragment_name, interval, accepted_name, tetrapod_group) %>%
	summarise(species_abundance = n())
# Create the csv of global biodiversity metrics
# Find the maximum number of individuals for any one site, per interval
max_number_individuals <- info_per_pcoord %>% group_by(interval, tetrapod_group) %>%
	summarise(max_no = max(individuals_total))
# Find the landscape richness per interval
landscape_richness <- site_biodiversity %>% group_by(interval, tetrapod_group) %>%
	summarise(richness = length(unique(accepted_name)),
		  number_individuals = n()) %>%
	full_join(info_per_pcoord %>% group_by(interval, tetrapod_group) %>%
		  	summarise(no_collections=sum(collections)))

# Find the alpha diversity per interval (mean of each site's richness)
alpha <- site_biodiversity %>% group_by(fragment_name, interval, tetrapod_group) %>%
	summarise(alpha_ind=length(unique(accepted_name))) %>%
	group_by(interval, tetrapod_group) %>% summarise(alpha = mean(alpha_ind))

# Create a csv of these biodiversity metrics per interval
landscape_richness <- landscape_richness %>% full_join(alpha)
landscape_richness$beta <- landscape_richness$richness/landscape_richness$alpha
write.csv(landscape_richness, file = file.path("input", "interval_biodiversity_metrics.csv"))
# Define the colours for the plot
bw_colours <- rep(c("grey60", "grey80"), 5)
ggthemr("flat dark")
colours <- ggthemr::swatch()
colours[3] <- colours[5]
# Define the plot labels
plot_labels <- data.frame(
	metric = factor("richness", levels=c("richness", "alpha",
						     "number_individuals", "no_collections")),
	interval=fct_recode(interval_data$interval, "Ks" = "Kasimovian", "Gz" = "Gzhelian",
			    "As" = "Asselian"),
	midpoint=interval_data$midpoint
)
p <- landscape_richness %>%
	full_join(interval_data) %>%
	gather(key="metric", value="value", richness, number_individuals, no_collections, alpha,
	       beta) %>%
	mutate(tetrapod_group = tools::toTitleCase(as.character(tetrapod_group)),
	       metric=factor(metric, levels=c("alpha", "beta","richness",
	       			       "number_individuals", "no_collections"))) %>%
	ggplot()+
	theme_classic() +
	scale_x_reverse("Time (Ma)")+
	scale_y_continuous("Value")+
	facet_grid(metric~.,
		   labeller = labeller(metric=as_labeller(c("richness" = "Species richness",
		   					 "alpha" = "Alpha diversity",
		   					 "beta" = "Beta diversity",
		   					 "number_individuals" = "No. individuals",
		   					 "no_collections" = "No. collections"))),
		   scales = "free_y") +
	scale_linetype("Tetrapod group")+
	scale_colour_manual("Tetrapod group", values=colours[2:6])+

	geom_rect(data=interval_data,
		  aes(x=NULL, y=NULL,
		      xmin=min_ma, xmax=max_ma, ymin=0, ymax=Inf,
		      group=interval, fill=interval), colour=NA, alpha=0.5)+
	geom_text(data=plot_labels, aes(x=midpoint, y = 120, label=interval),
		  colour="black", size=2,  angle=0, vjust=1, hjust=0.5)+
	geom_line(aes(x=midpoint, y=value, linetype=tetrapod_group, colour=tetrapod_group))+
	geom_vline(aes(xintercept=307), colour="grey10", linetype="dotted")+
	scale_fill_manual("", values=bw_colours[2:9], guide=FALSE)
pdf(file.path(figure_dir, "raw_diversity.pdf"), 6.5, 6)
print(p)
dev.off()
# Save a copy to thesis dir as well
pdf(file.path(figure_dir2, "paleo_raw_diversity.pdf"), 6.5, 6)
print(p)
dev.off()
# Generate the SQLite databases for comparing simulation outputs to
intervals <- unique(landscape_richness$interval)
tetrapod_groups <- unique(landscape_richness$tetrapod_group)
if(!dir.exists(file.path(data_dir, "databases")))
{
	dir.create(file.path(data_dir, "databases"))
}
for(each in intervals)
{
	for(tet_group in tetrapod_groups){
		interval_df <- site_biodiversity %>% filter(interval == each,
							    tetrapod_group == tet_group)
		fragment_abundances <- subset(interval_df, select=c(fragment_name,
								    accepted_name,
								    species_abundance))
		names(fragment_abundances) <- c("fragment", "species_id", "no_individuals")
		fragment_abundances$community_reference <- 0
		fragment_richness <- interval_df %>% group_by(fragment_name) %>%
			summarise(richness=n(),
				  no_individuals=sum(species_abundance))
		names(fragment_richness) <- c("fragment", "richness", "no_individuals")
		fragment_richness_out <- subset(fragment_richness, select=c("fragment", "richness"))
		fragment_richness$community_reference <- 0
		species_abundances <- interval_df %>% group_by(accepted_name) %>%
			summarise(no_individuals = sum(species_abundance))
		names(species_abundances) <- c("species_id", "no_individuals")
		species_abundances$community_reference <- 0
		plot_data <- interval_df %>% ungroup() %>%
			dplyr::select(fragment_name, species_abundance) %>%
			rename(no_individuals = species_abundance)
		names(plot_data) <- c("fragment", "no_individuals")
		db_name <- paste(file.path(data_dir, "databases/"),
				 tolower(gsub("\ ", "_", each)), "_",
				 tet_group,".db", sep="")
		if(file.exists(db_name))
		{
			file.remove(db_name)
		}
		db = dbConnect(SQLite(), dbname=db_name)
		dbWriteTable(db, "FRAGMENT_RICHNESS", fragment_richness_out)
		dbWriteTable(db, "FRAGMENT_ABUNDANCES", fragment_abundances)
		dbWriteTable(db, "SPECIES_ABUNDANCES", species_abundances)
		dbWriteTable(db, "PLOT_DATA", plot_data)
		fragment_richness$metric <- "fragment_richness"
		colnames(fragment_richness)[colnames(fragment_richness) == 'richness'] <- 'value'
		biodiversity_metrics <-  dplyr::select(fragment_richness, -community_reference)
		# browser()
		spec_richness <- data.frame(length(unique(species_abundances$species_id)),
					    sum(species_abundances$no_individuals))
		spec_richness$fragment <- "whole"
		spec_richness$metric <- "fragment_richness"
		names(spec_richness) <- c("value", "no_individuals", "fragment", "metric")
		biodiversity_metrics <- rbind(biodiversity_metrics, spec_richness)
		this_landscape <- landscape_richness %>% filter(interval == each,
								tetrapod_group == tet_group)
		# browser()
		alpha <- data.frame("alpha_diversity", spec_richness$no_individuals[1],
				    this_landscape$alpha, "whole")
		beta <- data.frame("beta_diversity", spec_richness$no_individuals[1],
				   this_landscape$beta, "whole")
		names(alpha) <- c("metric", "no_individuals", "value", "fragment")
		names(beta) <- names(alpha)
		biodiversity_metrics <- rbind(biodiversity_metrics, alpha, beta)
		dbWriteTable(db, "BIODIVERSITY_METRICS", biodiversity_metrics)
		dbDisconnect(db)
	}
}
