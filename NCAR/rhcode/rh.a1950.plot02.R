#################################################
## Plotting at all stations after 1950
## 
## - plot021:
##     dignostic plots conditional on year*month
##     all stations 7,738 are included
## - plot02:
##     dignostic plots condtional on stations
##     sampled stations from kd-tree built based on 4,978 stations
##     have over 300 obs
##
## - loess02: is set to be without elevation
## - loess04: is set to be with elevation
#################################################


# plot the residuals only at stations
library(maps)
library(lattice)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
#par$family <- "gaussian"
par$span <- 0.05
par$family <- "symmetric"
#par$span <- 0.025
par$loess <- "loess02.bystation.10pc"
#par$loess <- "loess04.bystation.10pc"
if(pmatch(par$loess, "loess02.bystation.10pc", nomatch = 0)){
	title <- paste("from Robust Fit without Elev span=", par$span, sep="")
}else{
	title <- paste("from Robust Fit with Elev span=", par$span, sep="")
}
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

rst <- rhread(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), par$loess
	)
)
result <- do.call("rbind", lapply(rst, "[[", 2))
result$residual <- result$tmax - result$fitted
result$station.id <- as.character(result$station.id)
result$month <- factor(
	result$month, 
	levels = c(
		"Jan","Feb","Mar","Apr","May","June",
		"July","Aug", "Sep", "Oct", "Nov", "Dec"
	)
)


################################
#Plots for sampled 128 stations#
################################
#order stations by the mean of the residual for sampled stations
load(file.path(local.datadir, "samplestation.a1950.RData"))
sub <- result[with(result, station.id %in% tmax.sample.a1950$station.id),]
sub$station.id <- reorder(
	sub$station.id, 
	sub$residual, 
	function(r){mean(r, na.rm=TRUE)}
)
## change the labels of stations on the strip
tmax.sample.a1950$station.id <- factor(
	tmax.sample.a1950$station.id, 
	levels=levels(sub$station.id)
)
striplabel <- tmax.sample.a1950[order(tmax.sample.a1950$station.id),]$leaf

## QQ plot of residuals for sample locations from kd-tree
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.dist.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- qqmath(~ residual | factor(station.id, labels = striplabel),
	data = sub,
	distribution = qnorm,
	layout = c(5,3),
	pch  = 16,
	aspect = 1,
	cex  = 0.3,
	main = list(
		label = paste("Normal Quantiles of Residuals", title)
	),
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

## scatter plot residual against fitted value for sampled stations
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.vs.fit.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- xyplot( residual ~ fitted | factor(station.id, labels=striplabel),
	data = sub, 
	layout = c(4,3),
#  par.settings = list(layout.heights = list(strip = 1.5)),
	pch  = 16,
	cex  = 0.3,
  key=list(
    text = list(label=c(
      "residuals",
      "loess smoothing: span=0.75, degree=2"
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
	scale = list(
		y = list(relation = "free"),
		x = list(relation = "free")
	),
	main = list(
		label = paste("Residuals vs. Fitted", title)
	),
	ylab = list(label = "Loess Residuals"),
	xlab = list(label = "Loess Fitted Value"),
	panel = function(x, y,...) {
		panel.xyplot(x, y, ...)
		panel.loess(x,y,span=0.75, degree=2, col=col[2],...)
		panel.abline(h=0, col = "red", lwd = 0.5, lty = 2)
	}
)
print(a)
dev.off()

## scatter plot of residual over time for sampled stations
sub <- sub[order(sub$station.id, sub$time),]
sub$time1 <- rep(c(rep(0:119, 4), 0:95), times = 128)
sub$factor <- factor(
    rep(rep(paste("Period", 1:5), c(rep(120,4),96)), times=128), 
    levels=paste("Period", c(5:1))
)
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.vs.time.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in levels(sub$station.id)) {
  a <- xyplot(residual ~ time1 | factor,
  	data = subset(sub, station.id == i),
		layout = c(1,5),
		strip = FALSE,
		pch  = 16,
#  	aspect = "xy",
		cex  = 0.4,
		main = list(
			label = paste("Residuals vs. Month", title)
		),
		ylab = list(
			label = paste("Loess residuals")
		),
		sub = paste("Station from cell",tmax.sample.a1950$leaf[with(tmax.sample.a1950, which(station.id == i))]),
		xlab = list(label = "Month"),
#		prepanel = function(x ,y){
#			index <- findInterval(mean(x), c(0,121,241,361,481,600))
#			lim <- list(c(0, 120), c(121, 240), c(241, 360), c(361, 480), c(481, 600))
#			list(xlim = lim[[index]])
#		},
		panel = function(x, y,...) {
			panel.xyplot(x, y,...)
			panel.abline(h=0, col="red", lty=2, lwd=0.8)
		}
	)
  print(a)
}
dev.off()
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.vs.time2.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in levels(sub$station.id)) {
  a <- xyplot(residual ~ time,
  	data = subset(sub, station.id == i&!is.na(residual)),
		layout = c(1,1),
		strip = FALSE,
		pch  = 16,
		cex  = 0.4,
		scale = list(
			x = list(
				relation = "free",
				tick.number = 10
			)
		),
  	key = list(
    	text = list(label=c(
    		"residuals",
      	"span=0.35, degree=1",
      	"span=0.15, degree=2"
    	)),
    	lines = list(
      	lwd=1.5, 
      	cex=0.7,
      	pch=16,
      	type=c("p","l","l"), 
      	col=col[1:3]
    	), 
    	columns = 3
  	),
		main = list(
			label = paste("Residuals vs. Month", title)
		),
		sub = paste("Station from cell",tmax.sample.a1950$leaf[with(tmax.sample.a1950, which(station.id == i))]),
		ylab = list(label = "Loess residuals"),
		xlab = list(label = "Month"),
		panel = function(x, y,...) {
			panel.xyplot(x, y,...)
			panel.abline(h=0, col="red", lty=2, lwd=0.8)
			panel.loess(
				x, y, span = 0.35, degree = 1, 
				col=col[2], evaluation = 100
			)
			panel.loess(
				x, y, span = 0.15, degree = 2, 
				col=col[3], evaluation = 100
			)
		}
	)
  print(a)
}
dev.off()

#scatter plot of fit over time for sampled stations
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.fit.vs.time.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in levels(sub$station.id)) {
	data <- subset(sub, station.id == i)
  a <- xyplot(fitted ~ time1 | factor,
  	data = data,
		layout = c(1,5),
		strip = FALSE,
  	aspect = "xy",
		scale = list(
      y = list(
      	relation = 'same', 
      	tick.number=4, 
      	alternating=TRUE
      ), 
      x = list(
      	at = seq(0, 143, by=12), 
      	relation = 'same'
      )
    ),
  	key = list(
    	text = list(label=c(
      	"Raw obs",
      	"Fitted value"
    	)),
    	lines = list(
      	pch=1, 
      	cex=0.5, 
      	lwd=1.5, 
      	type=c("p","l"), 
      	col=col[1:2]
    	), 
    	columns = 2
  	),
		main = list(
			label = paste("Residuals vs. Month", title)
		),
		sub = paste(
			"Station from cell",
			tmax.sample.a1950$leaf[with(tmax.sample.a1950, which(station.id == i))]
		),
  	xlim = c(0,119),
  	ylim = c(
  		min(min(data$tmax, na.rm = TRUE), min(data$fitted)) - 1,
  		max(max(data$tmax, na.rm = TRUE), max(data$fitted)) + 1
  	),
  	ylab = list(label = "Loess fitted value"),
		xlab = list(label = "Month"),
		panel = function(x, y, subscripts ,...) {
      panel.xyplot(
      	x = sort(x), 
      	y = y[order(x)], 
      	type="b", col=col[2], pch=16,cex =0.5,lwd=1, ...
      )
			panel.xyplot(
				x = data[subscripts,]$time1, 
				y = data[subscripts,]$tmax,
				type="p", col=col[1], cex=0.5, ...
			)
      panel.abline(
      	v=seq(0,119, by=12), color="lightgrey", lty=3, lwd=0.5
      )
		}
	)
  print(a)
}
dev.off()

##scatter plot of residuals over year conditional on month in a year
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.vs.year.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in levels(sub$station.id)){
  b <- xyplot( residual ~ (as.numeric(year)-1950) | month,
    data = subset(sub, station.id==i),
    xlab = list(label = "Year"),
		ylab = list(label = "Loess residuals"),
		main = list(
			label = paste("Residuals vs. Year", title)
		),
		sub = paste(
			"Station from cell",
			tmax.sample.a1950$leaf[with(tmax.sample.a1950, which(station.id == i))]
		),
  	key = list(
    	text = list(label=c(
      	"residuals",
      	"loess smoothing: span=0.75, degree=1"
    	)),
    	lines = list(
      	pch=16, 
      	cex=0.5, 
      	lwd=1.5, 
      	type=c("p","l"), 
      	col=col[1:2]
    	), 
    	columns = 2
  	),
    pch = 16,
    cex = 0.5,
    layout = c(12,1),
    strip = TRUE,
    scales = list(
    	y = list(relation = 'same', alternating=TRUE), 
    	x = list(tick.number = 5, relation='same')
    ),
    panel = function(x,y,...){
      panel.abline(h = 0, color = "black", lty = 1)
      panel.xyplot(x,y,...)
		  panel.loess(x,y,span = 0.75, degree = 1, col = col[2],...)
    }
  )
  print(b)
}
dev.off()
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.vs.year2.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in levels(sub$station.id)){
  b <- xyplot( residual ~ (as.numeric(year)-1950) | month,
    data = subset(sub, station.id == i),
    xlab = list(label = "Year"),
		ylab = list(label = "Loess residuals"),
		main = list(
			label = paste("Residuals vs. Year", title)
		),
		sub = paste(
			"Station from cell",
			tmax.sample.a1950$leaf[with(tmax.sample.a1950, which(station.id == i))]
		),
	  type = "b",	     
    pch  = 16,
    cex  =0.5,
    layout = c(2,6),
    strip = TRUE,
    scales = list(
    	y = list(relation = 'same', alternating=TRUE), 
    	x = list(tick.number=10, relation='same')
    ),
    panel = function(x,y,...) {
      panel.abline(h = 0, color = "black", lty = 1)
      panel.xyplot(x,y,...)
    }
  )
  print(b)
}
dev.off()

##QQ plot of residuals condtional on month in a year for sampled stations
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.dist.month.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in levels(sub$station.id)){
  a <- qqmath(~ residual | month,
    data = subset(sub, station.id == i),
    distribution = qnorm,
    aspect = 1,
    pch = 16,
    cex = 0.5,
    layout = c(4,3),
    xlab = list(label="Unit normal quantile", cex=1.2),
		ylab = list(label = "Loess residuals"),
		main = list(
			label = paste("Residuals vs. Year", title)
		),
		sub = paste(
			"Station from cell",
			tmax.sample.a1950$leaf[with(tmax.sample.a1950, which(station.id == i))]
		),
    prepanel = prepanel.qqmathline,
    panel = function(x, y,...) {
      panel.grid(lty=3, lwd=0.5, col="black",...)
      panel.qqmathline(x, y=x)
      panel.qqmath(x, y,...)
    }
  )
  print(a)
}
dev.off()

############
#Other code#
############
## compare the stl fit with spatial loess fit 
rst <- rhread(
	file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess02.bystation.stl")
)
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