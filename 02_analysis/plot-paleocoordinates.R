library(tidyverse)


fixedspreadmaps <- c("Tr1","J6","K7","Pg1")
fixedspreadmaps <- as.character(fixed.intervals$bin)
pdf(paste(folder.name, "/figures/figure-global-occurrence-maps.pdf", sep = ""), width = 16, height = 23)
par(mfrow = c(5,2))
for (i in 1:length(fixedspreadmaps)) {
	pm_plot(interval = paste("scotese_", round_any(fixed.intervals[fixed.intervals$bin == fixedspreadmaps[i], "midpoint"], 10),"m", sep = ""), #sample(strsplit(fixed.intervals$stages[fixedspreadmaps[i], ], split = ", ")[[1]], 1)
		dat = int_data[[fixedspreadmaps[i]]],
		point.col = "total_occurrences",
		point.size = "fixed", fixed.cex = 0.5,
		mst.tree = global.mst.trees[[fixedspreadmaps[i]]],
		add.title = FALSE,
		grid.cell.size = 1,
		lo.col = "darkblue", mid.col = "pink", hi.col = "red", line.col = "black")
	text(0,85, fixedspreadmaps[i], cex = 2)
	par(xpd = T); text(-185,100,letters[i],cex = 1.5, font = 2,family = 'sans'); par(xpd = F)
}
graphics.off()





