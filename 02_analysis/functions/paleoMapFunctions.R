pm_getmap <- function(interval, colsea = "#00509010",
		      colland = "#66666660",
		      do.plot = TRUE){
	## we might hack this with "with" or "null" for avoiding NOTE on check: 'no visible binding for global variable'
	## see: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when

	load(paste("./data/", interval, ".rda", sep = ""))
	# if user does not set plot=FALSE plot the shape file
	# get the shape file with help function getShape to open lazyload data
	assign("shape", get(interval))

	if (do.plot == TRUE) {
		par(mar = c(0,0,0,0))
		plot(shape, col = "white", border = FALSE)
		rect(xleft = -180, xright = 180, ybottom = -90,
		     ytop = 90, col = colsea,
		     border = FALSE)
		plot(shape, col = colland, border = F, add = T)
	}
	# return the shape file
	shape
}

#@TODO make the plot point cex sizes scale more sensibly
#@TODO prevent points from overplotting when grid-cell sizes are increased
pm_plot <- function(interval, dat, mst.tree = NULL, grid.cell.size = NULL,
		    point.col = c('total_collections','total_occurrences','total_species','prop_singletons'),
		    point.size = c('fixed','total_species','total_occurrences','total_collections','prop_singletons'),
		    add.title = TRUE,
		    add.grid = TRUE,
		    point.pch = 16,
		    fixed.cex = 2,
		    colsea = "#00509010",
		    colland = "#66666660",
		    colpoints = "#99000020",
		    lo.col = "dodgerblue4", mid.col = "green", hi.col = "yellow", line.col = "darkred",
		    floor.cex = 0.5,
		    add.legend = TRUE,
		    cex=1){

	color.gradient <- function(x, colors = c(lo.col,mid.col,hi.col), colsteps = 20) {
		return( colorRampPalette(colors) (colsteps)[findInterval(x, seq(min(x),max(x), length.out = colsteps))] )
	}

	point.col <- match.arg(point.col, c('total_occurrences','total_collections','total_species','prop_singletons'))
	col.title <- c("No. Occurrences","No. Collections","No. Species","% Singletons"); names(col.title) <- c('total_occurrences','total_collections','total_species','prop_singletons'); col.title <- col.title[point.col]
	point.size <- match.arg(point.size, c('fixed','total_species','total_occurrences','total_collections','prop_singletons'))
	point.title <- c("", "No. Species", "No. Occurrences", "No. Collections", "% Singletons"); names(point.title) <- c('fixed','total_species','total_occurrences','total_collections','prop_singletons'); point.title <- point.title[point.size]
	#getting the shape file for the map and the data for plotting it on the map
	shape <- pm_getmap(interval = interval, do.plot = FALSE)
	#plotting the map and the data
	if (class(dat) == "data.frame") {
		if (!is.null(grid.cell.size)) {
			dat$paleolatdec <- round_any(dat$paleolatdec, grid.cell.size)
			dat$paleolngdec <- round_any(dat$paleolngdec, grid.cell.size)
			dat$paleocoord <- paste(dat$paleolatdec,dat$paleolngdec, sep = ",")
		}
		#defines size and axes of the plot
		dat$occurrence.binomial <- paste(dat$occurrence.genus_name,dat$occurrence.species_name)

		total_occurrences <- as.data.frame(table(dat[, c('paleocoord','occurrence_no')]$paleocoord))
		rownames(total_occurrences) <- total_occurrences$Var1; total_occurrences$Var1 <- NULL; colnames(total_occurrences) <- "total_occurrences"
		total_collections <- as.data.frame(table(dat[!duplicated(dat[c("collection_no","paleocoord")]), c('paleocoord','collection_no')]$paleocoord))
		rownames(total_collections) <- total_collections$Var1; total_collections$Var1 <- NULL; colnames(total_collections) <- "total_collections"
		total_species <- as.data.frame(table(dat[!duplicated(dat[c("occurrence.binomial","paleocoord")]), c('paleocoord','occurrence.binomial')]$paleocoord))
		rownames(total_species) <- total_species$Var1; total_species$Var1 <- NULL; colnames(total_species) <- "total_species"
		unique.pcoords <- rownames(total_species)
		prop_singletons <- sapply(unique.pcoords, function(x) sum(table(dat[dat$paleocoord == x, "occurrence.binomial"]) == 1) / nrow(dat[dat$paleocoord == x, ]))
		prop_singletons <- as.data.frame(prop_singletons)
		cell.data <- cbind(total_occurrences,total_collections,total_species,prop_singletons)

		unique.dat <- dat[!duplicated(dat$paleocoord),c("paleocoord","paleolatdec","paleolngdec")]; rownames(unique.dat) <- unique.dat$paleocoord; unique.dat$paleocoord <- NULL
		plotting.data <- merge(unique.dat,cell.data, by = "row.names")

		if (point.col == "total_occurrences") {
			cols <- color.gradient(log(plotting.data$total_occurrences))
			col.min <- log(min(plotting.data$total_occurrences))
			col.max <- log(max(plotting.data$total_occurrences))
			col.min.legend <- min(plotting.data$total_occurrences)
			col.max.legend <- max(plotting.data$total_occurrences)
		}
		if (point.col == "total_collections") {
			cols <- color.gradient(log(plotting.data$total_collections))
			col.min <- log(min(plotting.data$total_collections))
			col.max <- log(max(plotting.data$total_collections))
			col.min.legend <- min(plotting.data$total_collections)
			col.max.legend <- max(plotting.data$total_collections)
		}
		if (point.col == "total_species") {
			cols <- color.gradient(log(plotting.data$total_species))
			col.min <- log(min(plotting.data$total_species))
			col.max <- log(max(plotting.data$total_species))
			col.min.legend <- min(plotting.data$total_species)
			col.max.legend <- max(plotting.data$total_species)
		}
		if (point.col == "prop_singletons") {
			cols <- colorRampPalette(colors = c(lo.col,mid.col,hi.col)) (20)[findInterval(plotting.data$prop_singletons, seq(0,1, length.out = 20))]
			col.min <- 0
			col.max <- 1
			col.min.legend <- 0
			col.max.legend <- 1
		}


		if (point.size == "fixed") {
			cex.sizes <- rep(fixed.cex, nrow(plotting.data))
		}
		if (point.size == "total_occurrences") {
			cex.sizes <- log10(plotting.data$total_occurrences)*2
			cex.sizes.legend <- plotting.data$total_occurrences
		}
		if (point.size == "total_collections") {
			cex.sizes <- log10(plotting.data$total_collections)*2
			cex.sizes.legend <- plotting.data$total_collections
		}
		if (point.size == "total_species") {
			cex.sizes <- log10(plotting.data$total_species)*2
			cex.sizes.legend <- plotting.data$total_species
		}
		if (point.size == "prop_singletons") {
			cex.sizes <- ((1 - plotting.data$prop_singletons) * 3) + 0.1
			cex.sizes.legend <- 1 - plotting.data$prop_singletons
		}

		par(mar = c(1,1,1,1))
		plot(1, xlim = c(-180,180), ylim = c(-90,90), type = "n", asp = 1, bty = "n", axes = FALSE)
		plot(shape, col = "white", border = FALSE, add = T)
		rect(xleft = -180, xright = 180, ybottom = -90,
		     ytop = 90, col = colsea,
		     border = FALSE)
		plot(shape, col = colland, border = FALSE, add = T)

		if (add.grid == TRUE) {
			# if (is.null(grid.cell.size)) {grid.overlay <- 10} else {grid.overlay <- grid.cell.size}
			xpd <- par()$xpd
			par(xpd = FALSE)
			grid.overlay <- 10
			lng_x1 <- lng_x0 <- seq(from = -180, to = 180, by = grid.overlay)
			lng_y0 <- rep(-90, length(lng_x0)); lng_y1 <- -(lng_y0)
			segments(x0 = lng_x0, y0 = lng_y0, x1 = lng_x1, y1 = lng_y1, lwd = 0.5, lty = 3, col = "darkgrey")
			lat_y1 <- lat_y0 <- seq(from = -90, to = 90, by = grid.overlay)
			lat_x0 <- rep(-180, length(lat_y0)); lat_x1 <- -(lat_x0)
			segments(x0 = lat_x0, y0 = lat_y0, x1 = lat_x1, y1 = lat_y1, lwd = 0.5, lty = 3, col = "darkgrey")
			axis(side = 2, pos = -180, at = seq(from = -90, to = 90, by = 45))
			axis(side = 1, pos = -90, at = seq(from = -180, to = 180, by = 45))
			par(xpd = xpd)
		}

		if (!is.null(mst.tree)) {
			pairwise.pcoords <- mst.tree[c("from","to")]
			for (i in 1:nrow(mst.tree)) {
				x0 <- mean(subset(dat, collection_no == pairwise.pcoords[i,1])$paleolngdec)
				y0 <- mean(subset(dat, collection_no == pairwise.pcoords[i,1])$paleolatdec)
				x1 <- mean(subset(dat, collection_no == pairwise.pcoords[i,2])$paleolngdec)
				y1 <- mean(subset(dat, collection_no == pairwise.pcoords[i,2])$paleolatdec)
				segments(x0, y0, x1, y1, col = line.col)
			}
		}
		points(plotting.data$paleolngdec,
		       plotting.data$paleolatdec,
		       pch = point.pch, col = adjustcolor(cols, alpha.f = 0.5),
		       cex = cex.sizes + floor.cex
		       # cex = cex*2
		)

		if (add.legend == TRUE) {
			SDMTools::legend.gradient(pnts = cbind(x = c(-175,-165,-165,-175), y = c(-85,-85,-45,-45)),
						  color.gradient(seq(from = col.min, to = col.max, length.out = 20)),
						  limits = round(range(c(col.min.legend, col.max.legend)),1),
						  title = col.title
			)

			if (point.size != "fixed") {
				legend.seq <- seq(from = min(cex.sizes), to = max(cex.sizes), length.out = 5)
				legend.seq.legend <- seq(from = min(cex.sizes.legend), to = max(cex.sizes.legend), length.out = 5)
				legend(-175, 0, pch = point.pch, pt.cex = legend.seq + floor.cex, legend = paste("   ", round_any(legend.seq.legend,1), sep = ""), bty = "n", col = adjustcolor("blue", alpha.f = 0.5), title = point.title)
			}

		}
		if (add.title == TRUE) {
			text(0,85, interval, cex = 2)
		}

	}
}


pm_plotCoords <- function(interval, dat, mst.tree,
			  colsea = "#00509010",
			  colland = "#66666660",
			  colpoints = "#99000020",
			  point.pch = 16,
			  cex = 1){

	#getting the shape file for the map and the data for plotting it on the map
	shape <- pm_getmap(interval = interval, do.plot = FALSE)
	#plotting the map and the data
	if (class(dat) == "data.frame") {
		#defines size and axes of the plot
		par(mar = c(0,0,0,0))
		plot(shape, col = "white", border = FALSE, xlim = c(-180,180), ylim = c(-90,90))
		rect(xleft = -180, xright = 180, ybottom = -90,
		     ytop = 90, col = colsea,
		     border = TRUE)
		plot(shape, col = colland, border = FALSE, add = T)
		points(dat$lngdec,
		       dat$latdec,
		       pch = point.pch, col = colpoints,
		       cex = cex)
		if (!is.null(mst.tree)) {
			pairwise.coords <- mst.tree[c("from","to")]
			for (i in 1:nrow(mst.tree)) {
				x0 <- mean(subset(dat, collection_no == pairwise.coords[i,1])$lngdec)
				y0 <- mean(subset(dat, collection_no == pairwise.coords[i,1])$latdec)
				x1 <- mean(subset(dat, collection_no == pairwise.coords[i,2])$lngdec)
				y1 <- mean(subset(dat, collection_no == pairwise.coords[i,2])$latdec)
				segments(x0, y0, x1, y1, col = "blue")
			}
		}
	}
}
