#################################################
##spatial loess fit at all stations after 1950 ##
#################################################

# loess02 is set to be span=0.05, degree=2
# plot the residuals only at stations
library(maps)
library(lattice)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$loess <- "loess02.bystation.10pc"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

rst <- rhread(file.path(rh.datadir, par$dataset, "spatial", "a1950", par$loess))
result <- do.call("rbind", lapply(rst, "[[", 2))
result$residual <- result$tmax - result$fitted
result$month <- factor(
	result$month, 
	levels = c(
		"Jan","Feb","Mar","Apr","May","June",
		"July","Aug", "Sep", "Oct", "Nov", "Dec"
	)
)
## QQ plot of residuals overall
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.dist.overall.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- qqmath(~ residual,
	data = result,
	distribution = qnorm,
# par.settings = list(layout.heights = list(strip = 1.5)),
	pch  = 16,
	aspect = 1,
	cex  = 0.5,
	xlab = list(label = "Unit normal quantile"),
	ylab = list(label = "Loess Residuals"),
	prepanel = prepanel.qqmathline,
	panel = function(x, y,...) {
			panel.grid()
			panel.qqmathline(x, y=x)
			panel.qqmath(x, y, ...)
	}
)
print(a)
dev.off()

## QQ plot of residuals for each location
# order stations by the diviation from the normal distribution
myfun <- function(data){
  yy <- quantile(data$residual, c(0.25, 0.75))
  xx <- qnorm(c(0.25, 0.75))
  r <- diff(yy)/diff(xx)
  x <- qnorm(ppoints(length(data$residual)))
  y <- r*x + yy[1] - xx[1]*r
  div <- sum(abs(sort(data$residual) - y))
}
#order stations by the mean of the residual for each station
result$station.id <- reorder(
	result$station.id, 
	result$residual, 
	function(r){abs(mean(r))}
)
#or order stations by the lat and lon
od <- as.character(
	unique(result[order(result$lon, result$lat, decreasing = TRUE), ]$station.id)
)
result$station.id <- factor(result$station.id, levels=od)

trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.dist.bystation2.", #bystation is ordered by mean, bystation2 is ordered by lat and lon
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- qqmath(~ residual | station.id,
	data = result,
	distribution = qnorm,
	layout = c(5,3),
#  par.settings = list(layout.heights = list(strip = 1.5)),
	pch  = 16,
	aspect = 1,
	cex  = 0.5,
	scale = list(y=list(relation="free")),
	xlab = list(label = "Unit normal quantile"),
	ylab = list(label = "Loess Residuals"),
#	prepanel = prepanel.qqmathline,
	panel = function(x, y,...) {
			panel.grid()
			panel.qqmathline(x, y=x)
			panel.qqmath(x, y, ...)
	}
)
print(a)
dev.off()

## QQ plot of residual for each time point
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.dist.bytime2.",  #bytime is year*month
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- qqmath(~ residual | month*year,
	data = result,
	distribution = qnorm,
	layout = c(6,2),
#  par.settings = list(layout.heights = list(strip = 1.5)),
	pch  = 16,
#	aspect = 1,
	cex  = 0.5,
	xlab = list(label = "Unit normal quantile"),
	ylab = list(label = "Loess Residuals"),
	prepanel = prepanel.qqmathline,
	panel = function(x, y,...) {
			panel.grid()
			panel.qqmathline(x, y=x)
			panel.qqmath(x, y, ...)
	}
)
print(a)
dev.off()


## scatter plot of residual against fited value conditional on month
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.vs.fit.bytime.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- xyplot(residual ~ fitted | month*year,
	data = result,
	layout = c(4,3),
#  par.settings = list(layout.heights = list(strip = 1.5)),
	pch  = 16,
#  aspect = "xy",
	cex  = 0.2,
	ylab = list(label = "Loess Residuals"),
	xlab = list(label = "LOess Fitted Value"),
	panel = function(x, y,...) {
			panel.xyplot(x, y, ...)
	}
)
print(a)
dev.off()

result$station.id <- reorder(
	result$station.id, 
	result$residual, 
	function(r){abs(mean(r))}
)
## scatter plot residual against fitted value for each station
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.vs.fit.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- xyplot( residual ~ fitted | station.id,
	data = result, 
	layout = c(4,3),
#  par.settings = list(layout.heights = list(strip = 1.5)),
	pch  = 16,
#  aspect = "xy",
	cex  = 0.2,
#	scale = list(y=list(relation="free")),
	ylab = list(label = "Loess Residuals"),
	xlab = list(label = "Loess Fitted Value"),
	panel = function(x, y,...) {
			panel.xyplot(x, y, ...)
	}
)
print(a)
dev.off()

## scatter plot of residual/fitted over time for each station
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.fit.vs.time.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in unique(result$station.id)) {
  a <- xyplot(fitted ~ time | cut(time, seq(0, 600,by=120)),
  	data = subset(result, station.id == i),
		layout = c(1,5),
		strip = FALSE,
		pch  = 16,
  	aspect = "xy",
		cex  = 0.4,
		scale = list(
			y = list(
				relation = "same", 
				alternating = TRUE
			),
			x = list(
				at = seq(0, 576, by = 12),
				relation = "free",
				axs = "r"
			)
		),
  	main = list(label = paste("Station", i)),
		ylab = list(label = "Loess Fitted Value"),
		xlab = list(label = "Month"),
		prepanel = function(x ,y){
			index <- findInterval(mean(x), c(0,121,241,361,481,600))
			lim <- list(c(0, 120), c(121, 240), c(241, 360), c(361, 480), c(481, 600))
			list(xlim = lim[[index]])
		},
		panel = function(x, y,...) {
			panel.xyplot(x, y,...)
		}
	)
  print(a)
}
dev.off()

## scatter plot of residual against lon/lat conditional on interval of lat/lon
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.vs.lat.bytime.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in sort(unique(result$year))) {
	for(j in levels(result$month)) {
		a <- xyplot(residual ~ lat | equal.count(lon, 20, overlap=0),
			data = subset(result, year == i & month == j),
			layout = c(10,2), #for vs.lon layout is c(10,2)
			strip=strip.custom(var.name = "Longitude", strip.levels=rep(FALSE, 2)),
			pch  = 16,
			cex  = 0.4,
			scale = list(
				y = list(
					relation = "same", 
					alternating = TRUE
				)
			),
    	key=list(
        text = list(label=c(
        	"residuals",
          "loess smoothing: span=0.75, degree=2 "
        )),
        lines = list(
            pch=c(".",""), 
            cex=4, 
            lwd=1.5, 
            type=c("p","l"), 
            col=col[1:2]
        ), 
        columns = 2
    	),
  		main = list(label = paste(i, j)),
			ylab = list(label = "Loess Residuals"),
			xlab = list(label = "Latitude"),
			prepanel = function(x,y,...) {
				prepanel.loess(x, y, span = 3/4, degree=2)
			},
			panel = function(x, y,...) {
				panel.xyplot(x, y,...)
				panel.loess(x,y,span = 3/4, degree=2, col=col[2],...)
			}
		)
		print(a)
	}
}
dev.off()

## compare the stl fit with spatial loess fit 
rst <- rhread(file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess02.bystation.stl"))
result <- do.call("rbind", lapply(rst, "[[", 2))
result$spatial.r <- result$tmax - result$fitted
result$stl.r <- result$tmax - result$stlfit
result$station.id <- reorder(
	result$station.id, 
	result$stl.r, 
	function(r){-length(r)}
)
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.vs.remainder.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
	a <- xyplot(spatial.r ~ stl.r | station.id,
		data = result,
		layout = c(5,2),
		pch  = 16,
		cex  = 0.4,
		aspect = 1,
		scale = list(
			y = list(
				relation = "same", 
				alternating = TRUE
			)
		),
		ylab = list(label = "Spatial loess Residuals"),
		xlab = list(label = "STL+ Remainders"),
		panel = function(x, y,...) {
			panel.xyplot(x, y,...)
			panel.abline(a=0, b=1, color="black", lty=1)
		}
	)
	print(a)
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