
##### PLOT PALEOZOIC TETRAPOD DATA ONTO SCOTESE PALEOMAPS FOR SENM PROJECT WITH EMMA DUNNE, SAM THOMPSON AND JAMES ROSINDELL ####

# SCRIPT WRITTEN BY ROGER CLOSE ON 22 SEPTEMBER 2017
# Last modified: April 10th 2019 by Sam Thompson


# install.packages(c("tidyverse","rgdal","paleobioDB","mapproj)) #install packages if necessary
library(tidyverse)
library(mapproj)
library(rgdal)
library(paleobioDB)


reload.scotese = T #reload Scotese spatial polygons from .shp files? (I left this code in so you can see how to load spatial polygon objects from shapefiles)
if (reload.scotese == F) {
	dirs <- list.dirs("input/scotese", recursive = FALSE, full.names = TRUE)
	scotese.ages <- list.dirs("input/scotese", recursive = FALSE, full.names = FALSE)

	scotese.shapes <- list()

	for (i in 1:length(dirs)) {

		scotese.shapes[[i]] <- readOGR(dsn = dirs[i],layer="wcnt")

	}; names(scotese.shapes) <- paste("scotese_", scotese.ages, sep = "")

	save(scotese.shapes, scotese.ages, file = "./input/scotese.shapes.RData")

} else {load("./input/scotese.shapes.RData")}

#make a named vector of numeric ages for Scotese maps
scotese.ages.numeric <- gsub(pattern = "m", replacement = "", x = scotese.ages) %>% as.integer
names(scotese.ages.numeric) <- names(scotese.shapes)


#read Emma's data file containing a list of collections per interval for each tetrapod group (amphibians and amniotes)
tetdat <- read_csv("./input/tetrapod_localities_Apr30.csv")


#read pbdb download that has Scotese paleocoordinates, then merge them into tetdat
scotese_pcoords <- read_csv("./input/scotese_emma_pbdb_data.csv", skip = 21)
scotese_pcoords <- distinct(scotese_pcoords, collection_no, paleolng, paleolat)
scotese_pcoords <- filter(scotese_pcoords, collection_no %in% tetdat$collection_no)
scotese_pcoords <- dplyr::rename(scotese_pcoords, paleolng_scotese = paleolng, paleolat_scotese = paleolat)
tetdat <- filter(tetdat, collection_no %in% scotese_pcoords$collection_no)
tetdat <- left_join(tetdat, scotese_pcoords, by = "collection_no")

# Also need the raw scotese data to know the species IDS at each time point
scotese_raw <- read_csv("./input/scotese_emma_pbdb_data.csv", skip = 21) %>%
	dplyr::select(collection_no, identified_no, accepted_no, paleolng, paleolat)

#make rounded palaeocoordinate variables
tetdat$paleolat_rounded <- plyr::round_any(tetdat$paleolat_scotese, 0.01) #number is nearest n degrees to round to
tetdat$paleolng_rounded <- plyr::round_any(tetdat$paleolng_scotese, 0.01)

#define intervals for plotting data onto maps - 8 intervals in total (C/P boundary is at 298.9)
ints <- tibble(interval = c("Bashkirian","Moscovian","Kasimovian", "Gzhelian", "Asselian", "Sakmarian", "Artinskian", "Kungurian"),
	       min_ma = c(315.2, 307.0, 303.4, 298.9, 295.5, 290.1, 279.3, 272.3),
	       max_ma = c(323.2, 315.2, 307.0, 303.4, 298.9, 295.5, 290.1, 279.3))
ints$midpoint <- (ints$max_ma + ints$min_ma) / 2


#subset tetdat into amniotes and ampbibians
amphibian_data <- filter(tetdat, tetrapod_group == "amphibian")
amniote_data <- filter(tetdat, tetrapod_group == "amniote")

write.csv(tetdat, file="./input/tetrapod_data_cleaned.csv")

# Generate lookup table of maps to use for each interval by finding the scotese map closest in age
# to time interval midpoint
intervals_data <- ints %>% rowwise() %>%
	mutate(map_file=scotese.ages[which(abs(scotese.ages.numeric - midpoint) ==
			   	min(abs(scotese.ages.numeric - midpoint)))])
write.csv(intervals_data, "./input/interval_data.csv")
#plot data onto maps - amphibians
if(FALSE){
	for (i in 1:nrow(ints)) {
		#find the scotese map closest in age to time interval midpoint
		use.me <- names(which(abs(scotese.ages.numeric - ints$midpoint[i]) == min(abs(scotese.ages.numeric - ints$midpoint[i]))))
		#filter amphibian_data to contain collections from this time interval only
		tmp <- filter(amphibian_data, ma_min >= ints$min_ma[i] & ma_max <= ints$max_ma[i])
		#get median lat and long for plotting on sphere
		medianlat <- median(tmp$paleolat_scotese)
		medianlng <- median(tmp$paleolng_scotese)

		ggplot() +
			geom_polygon(data = scotese.shapes[[use.me]], aes(x = long, y = lat, group = group), fill = "black", color = "black") +
			# coord_map("mercator") + #comment out all coord_map calls if you want to just plot maps using cartesian coordinates
			# coord_map("orthographic", orientation = c(medianlat, medianlng, 0)) + #spherical projection
			coord_map("mollweide") + #oooh pretty
			# coord_map("azequalarea", orientation = c(medianlat, medianlng, 0)) + #ooooh fancy but pointless
			geom_point(data = tmp, aes(x = paleolng_rounded, y = paleolat_rounded), col = "red", size = 0.75) +
			scale_y_continuous(breaks = seq(from = -90, to = 90, by = 10), limits = c(-90,90)) + #uncomment if lines are undesirable
			scale_x_continuous(breaks = seq(from = -180, to = 180, by = 10), limits = c(-180,180)) + #uncomment if lines are undesirable
			xlab('') + ylab('') +
			theme_minimal() +
			ggtitle(paste("Paleomap for ", ints$interval[i], sep = ""))
		file_name <- paste("./Paleomap for ", ints$interval[i], " amphibians.pdf", sep = "")
		ggsave(file_name, width = 10, height = 8)
	}


	#plot data onto maps - amniotes
	for (i in 1:nrow(ints)) {
		#find the scotese map closest in age to time interval midpoint
		use.me <- names(which(abs(scotese.ages.numeric - ints$midpoint[i]) == min(abs(scotese.ages.numeric - ints$midpoint[i]))))
		#filter amniote_data to contain collections from this time interval only
		tmp <- filter(amniote_data, ma_min >= ints$min_ma[i] & ma_max <= ints$max_ma[i])
		#get median lat and long for plotting on sphere
		medianlat <- median(tmp$paleolat_scotese)
		medianlng <- median(tmp$paleolng_scotese)

		ggplot() +
			geom_polygon(data = scotese.shapes[[use.me]], aes(x = long, y = lat, group = group), fill = "black", color = "black") +
			coord_map("mollweide") + #oooh pretty
			geom_point(data = tmp, aes(x = paleolng_rounded, y = paleolat_rounded), col = "red", size = 0.75) +
			scale_y_continuous(breaks = seq(from = -90, to = 90, by = 10), limits = c(-90,90)) + #uncomment if lines are undesirable
			scale_x_continuous(breaks = seq(from = -180, to = 180, by = 10), limits = c(-180,180)) + #uncomment if lines are undesirable
			xlab('') + ylab('') +
			theme_minimal() +
			ggtitle(paste("Paleomap for ", ints$interval[i], sep = ""))
		file_name <- paste("./Paleomap for ", ints$interval[i], " amniotes.pdf", sep = "")
		ggsave(file_name, width = 10, height = 8)
	}



	#plot data onto maps - all tetrapods (amphibians and amniotes together)
	for (i in 1:nrow(ints)) {
		#find the scotese map closest in age to time interval midpoint
		use.me <- names(which(abs(scotese.ages.numeric - ints$midpoint[i]) == min(abs(scotese.ages.numeric - ints$midpoint[i]))))
		#filter tetdat to contain collections from this time interval only
		tmp <- filter(tetdat, ma_min >= ints$min_ma[i] & ma_max <= ints$max_ma[i])
		#get median lat and long for plotting on sphere
		medianlat <- median(tmp$paleolat_scotese)
		medianlng <- median(tmp$paleolng_scotese)

		ggplot() +
			geom_polygon(data = scotese.shapes[[use.me]], aes(x = long, y = lat, group = group), fill = "black", color = "black") +
			coord_map("mollweide") + #oooh pretty
			geom_point(data = tmp, aes(x = paleolng_rounded, y = paleolat_rounded), col = "red", size = 0.75) +
			scale_y_continuous(breaks = seq(from = -90, to = 90, by = 10), limits = c(-90,90)) + #uncomment if lines are undesirable
			scale_x_continuous(breaks = seq(from = -180, to = 180, by = 10), limits = c(-180,180)) + #uncomment if lines are undesirable
			xlab('') + ylab('') +
			theme_minimal() +
			ggtitle(paste("Paleomap for ", ints$interval[i], sep = ""))
		file_name <- paste(".palaeomaps/Paleomap for ", ints$interval[i], " tetrapods.pdf", sep = "")
		ggsave(file_name, width = 10, height = 8)
	}
}



