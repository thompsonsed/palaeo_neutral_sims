# Maps of tetrapod samples from the Carboniferous and Permian.

library(tidyverse)
library(rgdal)
library(sp)
library(rgeos)
library(maptools) # Requires rgeos
library(viridis)
library(ggpubr)
map_input_dir <- "input/scotese/"
fig_dir <- "../../../Figures/Paleo/"
interval_data <- read.csv(file.path("input", "interval_data.csv"))
info_per_pcoord <- read.csv(file.path("input", "info_per_pcoord_main.csv")) %>% select(-X) %>%
	separate(pcoords, c("lat", "long"), sep = ",") %>%
	full_join(interval_data) %>%
	mutate(lat=as.numeric(lat),
	       long=as.numeric(long))

ps <- vector("list", ncol(interval_data))
for(i in 1:nrow(interval_data))
{
	cur_interval <- interval_data[i,]$interval


	tmp <- readOGR(dsn=file.path(map_input_dir, interval_data[i,]$map_file, "wcnt.shp"))
	tmp <- gBuffer(tmp, byid=TRUE, width=0)

	tmp@data$ID = rownames(tmp@data)
	tmp.points <- fortify(tmp, region="ID") %>% rename(ID=id)
	tmp.df <- tmp@data %>% inner_join(tmp.points, by="ID") %>%
		mutate(colour_code = as.numeric(PLATE_CODE)%%2)
	info_interval <- info_per_pcoord %>% filter(interval==cur_interval)
	p <- tmp.df %>% ggplot() +
		geom_polygon(aes(x=long,y=lat,group=group),
			     fill="grey70", colour="grey50",
			     show.legend = FALSE) +
		geom_point(data=info_interval, aes(x=long, y=lat,
						   colour=individuals_total,
						   size=individuals_total), alpha=0.6)+
		scale_colour_viridis("Number of fossils",  limits=c(0, 20),
				     breaks=seq(0, 20, 5), begin = 0.2, end=0.9)+
		scale_size_continuous("Number of fossils", limits=c(0, 20),
				      breaks=seq(0, 20, 5))+
		scale_x_continuous(breaks = seq(-180, 180, 45))+
				   # expand = expand_scale(mult = c(5, 5)))+
				   # seq(-180, 180, 45))+

		guides(colour = guide_legend(), size = guide_legend()) +
		# scale_fill_grey(start=0.4, end=0.7)+
		theme_bw()+
		xlab(NULL) + ylab(NULL)+
		coord_map("mollweide", xlim = c(-180, 180), ylim=c(-90, 90), clip = "off")+
		ggtitle(paste(cur_interval, ": ", interval_data[i,]$max_ma, "-",
			      interval_data[i,]$min_ma, " Ma", sep=""))+
		theme(
			# panel.background = element_rect(fill = "#BFD5E3",
			# 			      size = 2, linetype = "solid"),
		      plot.title = element_text(hjust = 0),
		      panel.grid=element_line(),
		      plot.margin = margin(unit(2, "lines"),
		      		     unit(2, "lines"),
		      		     unit(2, "lines"),
		      		     unit(2, "lines")),
		      panel.border=element_blank(),
		      axis.ticks=element_blank(),
		      axis.text=element_blank())
	# p
	ps[[i]] <- p
}


gga1 <- ggarrange(plotlist = ps, ncol = 2, common.legend = TRUE, legend = "bottom",
		  nrow=4)
# gga1
pdf(file.path(fig_dir, "paleo_tetrapod_maps.pdf"), 6.5, 9)
print(gga1)
dev.off()
