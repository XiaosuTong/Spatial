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
## 24lat + 10lon = -400 is the line to cut off the left lower corner
## 13lat - 20lon = 1925 is the line to cut off the right lower corner
result <- subset(
	x = result, 
	subset = resid.fit >= -3 & resid.fit <= 3 & 
	(24*lat + 10*lon) >= -380 & (13*lat - 20*lon) >= 1925
)
result$month <- factor(
	result$month, 
	levels = c(
		"Jan","Feb","Mar","Apr","May","June",
		"July","Aug", "Sep", "Oct", "Nov", "Dec"
	)
)
trellis.device(
	postscript, 
	file = paste(local.output, "/a1950.contour.loess.resid", par$dataset, ".ps", sep = ""), 
	color = TRUE,
	paper = "legal"
)
contourplot(
	resid.fit ~ lon * lat | month*year,
	data = result,
  cuts = 12, 
  region = TRUE,
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

trellis.device(
	postscript, 
	file = paste(local.output, "/a1950.loess.resid.vs.lon", par$dataset, ".ps", sep = ""), 
	color = TRUE,
	paper = "legal"
)
xyplot(
	resid.fit ~ lon | lat,
	data = result,
	xlab = "Longitude",
	ylab = "Smoothing Residual",
	layout = c(3, 3)
)
dev.off()
trellis.device(
	postscript, 
	file = paste(local.output, "/a1950.loess.resid.vs.lat", par$dataset, ".ps", sep = ""), 
	color = TRUE,
	paper = "legal"
)
xyplot(
	resid.fit ~ lat | lon,
	data = result,
	xlab = "Latitude",
	ylab = "Smoothing Residual",
	layout = c(3, 3)
)
dev.off()