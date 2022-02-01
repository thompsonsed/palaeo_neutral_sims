## Testing the upscaling in biodiversity for contiguous landscapes using the Preston function

library(tidyverse)
library(viridis)

source("preston.R")
fig_dir <- "../../../Figures/Paleo/"
df <- expand.grid(a_local=10^seq(2, 9, 1),
		  a_region=10^seq(3, 9, 1),
		  sigma=2^seq(1, 5, 1),
		  nu=10^seq(-4, -9, -0.1),
		  mult=c(10, 100, 1000))
df <- df %>% filter(a_region > a_local, nu*mult < 1) %>%
	mutate(richness_local_init=S_contig(a_local, nu, sigma^2),
	       richness_local_end=S_contig(a_local, nu*mult, sigma^2),
	       richness_region_init=S_contig(a_region, nu, sigma^2),
	       richness_region_end=S_contig(a_region, nu*mult, sigma^2),
	       prop_local = richness_local_end/richness_local_init,
	       prop_region = richness_region_end/richness_region_init,
	       se = prop_local/prop_region)

p <- df %>%
	mutate(region_mult = a_region/a_local) %>%
	filter(region_mult %in% c(10, 100, 1000),
	       a_local %in% 10^seq(2, 6, 1),
	       sigma==8) %>%
	ggplot()+
	geom_line(aes(x=nu, y=se*100, colour=a_local, group=a_local))+
	scale_colour_viridis(expression(paste("Size of local\ncommunity (", A[local], ")")),
			     trans="log10",
			     labels=scales::trans_format("log10", scales::math_format()))+
	scale_linetype_discrete(expression(paste(sigma)))+
	scale_x_log10("Speciation rate",
		      labels=scales::trans_format("log10", scales::math_format()))+
	scale_y_continuous("Sampling effectiveness (%)")+
	facet_grid(region_mult ~mult, labeller = labeller(region_mult=as_labeller(c("10"="A[region]==10*A[local]",
										    "100"="A[region]==100*A[local]",
										    "1000"="A[region]==1000*A[local]"), default = label_parsed),
							  mult=as_labeller(c("10"="nu[t[1]]==10*nu[t[0]]",
							  		   "100"="nu[t[1]]==100*nu[t[0]]",
							  		   "1000"="nu[t[1]]==1000*nu[t[0]]"), default = label_parsed)))+
	theme_classic()
pdf(file.path(fig_dir, "paleo_contiguous_detection.pdf"), 6.5, 5)
print(p)
dev.off()
