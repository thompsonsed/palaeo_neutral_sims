# Plotting of the conceptual figure

library(tidyverse)
library(viridis)
library(ggpubr)


input_folder <- "input"

# interval_biodiversity_metrics <- read.csv(file.path(input_folder, "interval_biodiversity_metrics.csv"))
interval_data <- read.csv(file.path(input_folder, "interval_data.csv")) %>%
	select(-X) %>%
	arrange(midpoint)
colour_1 <- "mediumseagreen"
colour_2 <- "sienna3"

p1 <- interval_data %>% ggplot() +
	geom_rect(aes(ymin=0, ymax=1, xmin=min_ma, xmax=max_ma, group=interval),
		  fill=rep(c("grey60", "grey80"), 4), alpha=0.5) +
	geom_vline(aes(xintercept=307), linetype="dotted")+
	geom_text(aes(x=midpoint, y=0.5, label=interval), size=2)+
	theme_classic() +
	theme(axis.line.y = element_blank(),
	      aspect.ratio=0.05,
	      axis.title.y=element_blank(),
	      axis.text.y=element_blank(),
	      axis.ticks.y=element_blank()) +

	scale_x_reverse("Time (Ma)")
	# scale_colour_discrete(values=)
	# geom_hline(aes(yintercept=0), color="black")

p2 <- interval_data %>% ggplot() +
	geom_rect(aes(ymin=0, ymax=1, xmax=max(max_ma), xmin=min(min_ma)), fill=colour_1) +
	geom_vline(aes(xintercept=307), linetype="dotted")+
	geom_text(aes(x=295, y=0.5), label="MODEL A") +
	theme_void() +
	theme(aspect.ratio=0.1)+
	scale_x_reverse()

p3 <- interval_data %>% ggplot() +
	geom_rect(aes(ymin=0, ymax=1, xmin=min(min_ma), xmax=307), fill=colour_2) +
	geom_rect(aes(ymin=0, ymax=1, xmax=max(max_ma), xmin=307), fill=colour_1) +
	geom_vline(aes(xintercept=307), linetype="dotted")+
	geom_text(aes(x=315, y=0.5), label="MODEL A") +
	geom_text(aes(x=290, y=0.5), label="MODEL B") +
	theme_void() +
	theme(aspect.ratio=0.1)+
	scale_x_reverse()

gga_out <- ggarrange(p1, p2, p3, nrow = 3, ncol=1, labels=c("", "Ignoring impact of CRC (Scenario A)", "Including impact of CRC (Scenarios B and C)"), hjust = 0, align = "v")

pdf(file.path("figures", "conceptual_figure.pdf"), 7, 4)
print(gga_out)
dev.off()

kable(interval_data)
