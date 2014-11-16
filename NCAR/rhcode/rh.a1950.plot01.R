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