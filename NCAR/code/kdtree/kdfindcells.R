#######################################
## Sample stations from a1950 stations
#######################################

## build a kdt-tree using the 7,738 stations after 1950
## randomly sample one station from each cell

library(plyr)
library(maps)
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$family <- "symmetric"
par$span <- 0.05
par$loess <- "loess04.bystation.all.10pc"
source("~/Projects/Spatial/NCAR/myloess/kdtree.R")
source("~/Rhipe/rhinitial.R")
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
load("~/Projects/Spatial/NCAR/RData/info.RData")
dyn.load("~/Projects/Spatial/NCAR/myloess/src/cppkdtree.so")
source("~/Projects/Spatial/NCAR/code/spatial/mykdtree.R")

## to get total stations after 1950
## there are 7,738 stations
rst <- rhread(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), par$loess
	)
)
result <- do.call("rbind", lapply(rst, "[[", 2))
stations.a1950.tmax <- as.character(unique(result$station.id))
## use c++ kdtree code to get the index of each obs in each cell
## locs is all stations in a1950
locs <- UStinfo[UStinfo$station.id %in% stations.a1950.tmax, ]
row.names(locs) <- 1:nrow(locs)
rst <- cppkdtree(as.matrix(locs[,c(3,4)]), 100)
## idx is a data.frame which sample one station for each cell
idx <- ddply(
	.data = rst,
	.variable = "leaf",
	.fun = function(r) {
		set.seed(100)
		sample(r$idx, 1)
	}
)
## places is stations with lon and lat in each cell
places <- locs[row.names(locs) %in% idx[,2],]
tmax.sample.a1950 <- data.frame(
	station.id = as.character(locs[row.names(locs) %in% idx[,2],1]),
	leaf = idx$leaf[order(idx$V1)]
)
save(
	list = "tmax.sample.a1950", 
	file = file.path(local.datadir, "samplestation.a1950.RData")
)

## vert is using my own kdtree function to get the vertices information
vert <- kdtree(locs[, c(3,4)], 0.0125)

us.map <- map('state', plot = FALSE, fill = TRUE)

trellis.device(
	device = postscript, 
	file = paste(local.output, "/map_vertices_", par$dataset, ".ps", sep = ""), 
	color=TRUE, 
	paper="legal"
)
a <- xyplot( lat ~ lon,
	data  = vert[[1]],
	xlab  = list(label = "Longitude"),
	ylab  = list(label = "Latitude"),
	main  = list(label = "Stations from kdtree"),
  key=list(
    text = list(label=c(
      "station location",
      "kdtree cell boundry"
    )),
    lines = list(
      pch = 16, 
      cex = 0.7, 
      lwd = 1,
      lty = 2, 
      type = c("p","l"), 
      col = c("blue", "red")
    ), 
    columns = 2
  ),
	panel = function(x, y, ...) {
		panel.xyplot(
			x = places$lon, 
			y = places$lat, 
			cex = 0.7, col = "blue", 
			pch = 16, type = "p", ...
		)
		for(k in 1:nrow(places)){
			panel.text(
				places$lon[k], places$lat[k], tmax.sample.a1950$leaf[k], 
				adj = c(0,0), col = col[1]
			)
		}
		for(i in seq(1, nrow(vert[[1]]), by = 2)){
			panel.segments(
				x[i], y[i], x[i+1], y[i+1], col = "red", lty = 2, lwd = 0.7
			)
		}
		panel.segments(
			c(x[1], x[2]), c(y[1], y[2]), 
			c(x[3], x[4]), c(y[3], y[4]), 
			col="red", lty=2, lwd = 0.7
		)
	  panel.polygon(us.map$x,us.map$y)	
	  panel.xyplot(
	  	x = x, 
	  	y = y, 
	  	col = "red", type = "p", 
	  	pch = 16, cex = 0.3, ...
	  )
	}
)
print(a)
dev.off()

## quantile plot of elevation in each cell of kdt-tree
locs$leaf <- rst$leaf[order(rst$idx)]
trellis.device(
	device = postscript, 
	file = paste(local.output, "/elev.dist.bycell.", par$dataset, ".ps", sep = ""), 
	color=TRUE, 
	paper="legal"
)
a <- qqmath(~ log2(elev+128) | factor(leaf),
	data = locs, 
	distribution = qunif,
	layout = c(5,3),
#  par.settings = list(layout.heights = list(strip = 1.5)),
	pch  = 16,
	aspect = 1,
	cex  = 0.5,
	xlab = list(label = "f-value"),
	ylab = list(label = "Log base 2 elevation"),
	prepanel = prepanel.qqmathline,
	panel = function(x, y,...) {
		panel.grid()
		panel.qqmath(x, y, ...)
	}
)
print(a)
dev.off()