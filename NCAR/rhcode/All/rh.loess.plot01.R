###################################
##Plotting the time series for each vertex of kd tree
##
##initialize the rhipe and setup the directories
source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$loess <- "loess01"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

rst <- rhread(file.path(rh.datadir, par$dataset, "spatial", par$loess))
result <- do.call("rbind", lapply(rst, "[[", 2))
result$factor <- factor(
	rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times=length(unique(result$fac))), 
	levels = paste("Period", c(9:1))
)
if(par$dataset == "precip") {
	ylab <- "Precipitation (millimeters)"
}else if(par$dataset == "tmax") {
	ylab <- "Maximum Temperature (degree centigrade)"
}else {
	ylab <- "Minimum Temperature (degree centigrade)"
}
result$time2 <- c(rep(0:143,8), 0:83)

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

order <- ddply(
	.data = result,
    .variables = "fac",
    .fun = summarise,
    lon = lon[1],
    lat = lat[1]
)
order.st <- as.character(
	order[order(order$lon, order$lat, decreasing=TRUE), ]$fac
)
result$fac <- factor(
	result$fac, 
	levels = order.st
)

trellis.device(
	device = postscript, 
	file = paste(local.output, "/scatterplot_of_vertices_", par$dataset, ".ps", sep = ""), 
	color = TRUE, 
	paper = "legal"
)
  for(i in levels(result$fac)) {
	b <- xyplot( fitted ~ time2 | factor,
    	data = subset(result, fac == i),
        xlab = list(label = "Month"),
        ylab = list(label = ylab),
		main = list(label = paste("Vertex of ", "(", 
			unique(subset(result, fac == i)$lat), ", ",
			unique(subset(result, fac == i)$lon), ")", sep="")
		),
		type = "b",
		pch = 16,
		cex = 0.5,
		layout = c(1,9),
		strip = FALSE,
		grib = TRUE,
		xlim = c(0, 143),
		scales = list(
			y = list(relation = 'same', alternating = TRUE), 
			x = list(at = seq(0, 143, by = 12), relation = 'same')
		),
		panel = function(...) {
			panel.abline(v = seq(0, 145, by = 12), color = "lightgrey", lty = 3, lwd =0.5)
			panel.xyplot(...)
		}
	)
	print(b)
  }
dev.off()

library(maps)
us.map <- map('state', plot = FALSE, fill = TRUE)
vertices <- result[seq(1, nrow(result), by=1236), c("lon", "lat", "fac")]
vertices$fac <- factor(
	vertices$fac, 
	levels = 1:length(unique(vertices$fac))
)
vertices <- vertices[with(vertices, order(fac)),]
trellis.device(
	device = postscript, 
	file = paste(local.output, "/map_vertices_", par$dataset, ".ps", sep = ""), 
	color=TRUE, 
	paper="legal"
)
a <- xyplot( lat ~ lon,
	data  = vertices,
	xlab  = list(label = "Longitude"),
	ylab  = list(label = "Latitude"),
	main  = list(label = "Vertices of K-D Tree"),
	type  = "p",
	cex   = 0.5,
	col   = "red",
	pch   = 16,
	panel = function(x, y, ...) {
		panel.xyplot(x, y, ...)
		for(i in seq(1, nrow(vertices), by = 2)){
			panel.segments(x[i], y[i], x[i+1], y[i+1], col = "red")
		}
		panel.segments(
			c(x[1], x[2]), c(y[1], y[2]), 
			c(x[3], x[4]), c(y[3], y[4]), 
			col="red"
		)
	    panel.polygon(us.map$x,us.map$y)	
	}
)
print(a)
dev.off()

result$fac <- factor(
	result$fac, 
	levels = order.st
)
trellis.device(
	device = postscript, 
	file = paste(local.output, "/QQ_plot_of_month_", par$dataset, ".ps", sep = ""), 
	color = TRUE, 
	paper = "legal"
)
  for(i in levels(result$fac)){
    a <- qqmath(~ fitted | month,
			data = subset(result, fac == i),
			distribution = qnorm,
			aspect = "xy",
			layout = c(12,1),
			pch  = 16,
			cex  = 0.5,
			main = list(label = paste("Vertex of ", "(",
				unique(subset(result, fac == i)$lat), ", ",
				unique(subset(result, fac == i)$lon), ")", sep="")
			),
			xlab = list(label = "Unit normal quantile"),
			ylab = list(label = ylab, cex=1.2),
#     scales = list(x = list(cex=1.5), y = list(cex=1.5)),
			prepanel = prepanel.qqmathline,
			panel = function(x, y,...) {
				panel.grid()
				panel.qqmathline(x, y=x)
				panel.qqmath(x, y, ...)
			}
    )
    print(a)
  }
dev.off()


dd <- ddply(
	.data = result,
	.variables = c("fac","month"),
	.fun = summarise,
	mean = mean(fitted)
)
mm <- dd[rep(row.names(dd), each=103),]
result <- result[with(result, order(fac, month, year)), ]
result$central <- result$fitted - mm$mean
trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/", par$dataset, 
		"_vertices_conditional_month.ps", 
		sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
  for(i in levels(result$fac)) {
    b <- xyplot( central ~ year | month,
		data = subset(result, fac == i),
		xlab = list(label = "Year"),
		ylab = list(label = ylab),
        main = list(label = paste("Vertex of ", "(",
            unique(subset(result, fac == i)$lat), ", ",
            unique(subset(result, fac == i)$lon), ")", sep="")
        ),
		type = "p",
		pch = 16,
		cex = 0.5,
		layout = c(4,3),
		strip = TRUE,
		key = list(
			text = list(
				label = c("centralized fitted","loess smoothing")
			), 
			lines = list(
				pch = c(".", ""), 
				cex = 4, 
				lwd = 1.5, 
				type = c("p","l"), 
				col = col[1:2]
			),
			columns = 2
		),		
		scales = list(
			y = list(relation = 'same', alternating=TRUE), 
			x = list(tick.number=10, relation='same')
		),
		panel = function(x, y, ...) {
			panel.xyplot(x, y, ...)
			panel.loess(x, y, span=2/3, degree=1, col=col[2],...)
		}
	)
	print(b)
  }
dev.off()
