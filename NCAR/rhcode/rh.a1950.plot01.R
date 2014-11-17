###################################################################
##spatial loess fit at all stations after 1950 with interpolation##
###################################################################

# loess01 is set to be span=0.05, degree=2
library(maps)
library(lattice)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col
us.map <- map('state', plot = FALSE, fill = TRUE)

source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$loess <- "loess01"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

rst <- rhread(file.path(rh.datadir, par$dataset, "spatial", "a1950", par$loess))
time <- data.frame(
	do.call("rbind",lapply(rst, "[[", 1)), 
	stringsAsFactors=FALSE
)
t <- time[c(rep(1:nrow(time), each=nrow(rst[[1]][[2]]))), ]
names(t) <- c("year", "month")
result <- cbind(t, do.call("rbind", lapply(rst, "[[", 2)))
result$month <- factor(
	result$month, 
	levels = c(
		"Jan","Feb","Mar","Apr","May","June",
		"July","Aug", "Sep", "Oct", "Nov", "Dec"
	)
)
result <- result[order(result$lon, result$lat),]
result$location <- factor(rep(1:nrow(rst[[1]][[2]]), each = 576))

## contour plot of smoothing residuals from spatial loess
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.contour.loess.resid.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
contourplot(
	resid.fit ~ lon * lat | month * year,
	data = result,
  region = TRUE,
  at = seq(-3, 3, 0.5),
  cuts = 10,
  layout = c(2,2),
  xlab = "Longitude",
  ylab = "Latitude",
  main = "Smoothing Loess Residuals",
 	panel = function(x, y, z, ...) {
 		panel.contourplot(x,y,z,...)
    panel.polygon(us.map$x,us.map$y, border = "red") 
 	}
 )
dev.off()

## QQ plot of smoothing residuals for each location
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.dist.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
qqmath(~ resid.fit | location,
	data = result,
	distribution = qnorm,
	layout = c(6,3),
	strip = strip.custom(
		factor.levels = paste(
			result$lon[seq(1,nrow(result),by=576)], 
			result$lat[seq(1,nrow(result),by=576)], 
			sep=", "
		)
	),
#  par.settings = list(layout.heights = list(strip = 1.5)),
	pch  = 16,
	aspect = 1,
	cex  = 0.5,
	xlab = list(label = "Unit normal quantile"),
	ylab = list(label = "Smoothing Loess Residuals"),
	prepanel = prepanel.qqmathline,
	panel = function(x, y,...) {
			panel.grid()
			panel.qqmathline(x, y=x)
			panel.qqmath(x, y, ...)
	}
)
dev.off()

## calculate the mean and std for each location and 
## plot the mean against lon and lat
result.mean <- ddply(
	.data = result,
	.variable = "location",
	.fun = summarize,
	mean = mean(resid.fit),
	std = sd(resid.fit),
	lon = unique(lon),
	lat = unique(lat)
)

trellis.device(
	postscript, 
	file = paste(
		local.output, 
		"/a1950.loess.residmean.vs.lon.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
xyplot(
	mean ~ lon | as.factor(lat),
	data = result.mean,
	aspect = "xy",
	pch = 1,
	cex = 0.7,
	xlab = "Longitude",
	ylab = "Mean of Smoothing Residual",
	layout = c(4,2),
	panel = function(x,y,...){
		panel.xyplot(x,y,...)
		panel.abline(h=0, col = "red", lty = 6)
	}
)
dev.off()

trellis.device(
	postscript, 
	file = paste(
		local.output, 
		"/a1950.loess.residmean.vs.lat.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
xyplot(
	mean ~ lat | as.factor(lon),
	data = result.mean,
	xlab = "Latitude",
	ylab = "Mean of Smoothing Residual",
	aspect = "xy",
	pch = 1,
	cex = 0.7,
	layout = c(5, 2),
	panel = function(x,y,...){
		panel.xyplot(x,y,...)
		panel.abline(h=0, col = "red", lty = 6)
	}
)
dev.off()

trellis.device(
	postscript, 
	file = paste(
		local.output, 
		"/a1950.loess.residstd.vs.lon.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
xyplot(
	std ~ lon | as.factor(lat),
	data = result.mean,
	aspect = "xy",
	pch = 1,
	cex = 0.7,
	xlab = "Longitude",
	ylab = "Standard Deviation of Smoothing Residual",
	layout = c(6, 2),
	panel = function(x,y,...){
		panel.xyplot(x,y,...)
	}
)
dev.off()

trellis.device(
	postscript, 
	file = paste(
		local.output, 
		"/a1950.loess.residstd.vs.lat.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
xyplot(
	std ~ lat | as.factor(lon),
	data = result.mean,
	xlab = "Latitude",
	ylab = "Standard Deviation of Smoothing Residual",
	aspect = "xy",
	pch = 1,
	cex = 0.7,
	layout = c(7, 1),
	panel = function(x,y,...){
		panel.xyplot(x,y,...)
	}
)
dev.off()

## preparing for contourplot for elevation
load("stations.a1950.RData")
load("info.RData")
dyn.load("~/Projects/Spatial/NCAR/myloess/shareLib/myloess2.so")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")

loc <- subset(UStinfo, station.id %in% stations.a1950.tmax)[,-5]
new.grid <- expand.grid(
	lon = seq(-124, -67, by = 0.25),
	lat = seq(25, 49, by = 0.25)
)
instate <- !is.na(map.where("state", new.grid$lon, new.grid$lat))
new.grid <- new.grid[instate, ]
elev.fit <- my.loess2( elev ~ lon + lat,
	data = loc,
	degree = 2, 
	span = 0.05
)
grid.fit <- my.predict.loess(
	object = elev.fit,
	newdata = data.frame(
		lon = new.grid$lon,
		lat = new.grid$lat
	)
)
new.grid$elev <- grid.fit

##plotting
trellis.device(
	postscript, 
	file = paste(
		local.output, "/contour.loess.elev.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
contourplot(
	elev ~ lon * lat,
	data = new.grid,
  region = TRUE,
#  at = seq(-3, 3, 0.5),
  cuts = 30,
  xlab = "Longitude",
  ylab = "Latitude",
  main = "Contourplot for Elevation",
 	panel = function(x, y, z, ...) {
 		panel.contourplot(x,y,z,...)
    panel.polygon(us.map$x,us.map$y, border = "blue") 
 	}
 )
dev.off()

##coast distance
library(sp)
library(maps)

## single point for a simple test
pts <- matrix(c(-88, 36), ncol = 2)
## simple map data
mp <- map("usa", plot = FALSE)
## convert coords to matrix and dump NA
xy.coast <- cbind(mp$x, mp$y)[!is.na(mp$x), ]
sea.coast <- xy.coast[c(195:4238,6095:6713),]

## container for all the nearest points matching the input
closest.points <- matrix(0, ncol = 2, nrow = nrow(pts))
for (i in 1:nrow(pts)) {
   closest.points[i, 1:2] <- sea.coast[which.min(spDistsN1(sea.coast,
pts, longlat = TRUE)), ]
}