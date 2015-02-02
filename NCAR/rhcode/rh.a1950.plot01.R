###################################################################
##spatial loess fit at all stations after 1950 with interpolation##
###################################################################

# loess01 is set to be span=0.05, degree=2 
# plot the residual smoothing surface using grid points
library(maps)
library(lattice)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col
us.map <- map('state', plot = FALSE, fill = TRUE)

source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$loess <- "loess03.elev"
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

## contour plot of smoothing residuals from spatial loess
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.contour.loess.resid.elev2.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in c(sort(unique(result$year)))){
	for(j in levels(result$month)){
		a <- contourplot(
			resid.fit ~ lon * lat,
			data = subset(result, year == i & month == j),
  		region = TRUE,
 		# at = seq(-3, 3, 0.5),
  		cuts = 20,
  		xlab = "Longitude",
  		ylab = "Latitude",
  		main = paste("Smoothing Loess Residuals", j, i),
 			panel = function(x, y, z, ...) {
 				panel.contourplot(x,y,z,...)
    		panel.polygon(us.map$x,us.map$y, border = "red") 
 			}
 		)
 		print(a)
 	}
}
dev.off()