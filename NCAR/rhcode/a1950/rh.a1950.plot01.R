###################################################################
##spatial loess fit at all stations after 1950 with interpolation##
###################################################################

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
par$loess <- "loess03.rob"
#par$family <- "gaussian"
#par$span <- 0.05
par$span <- 0.025
par$family <- "symmetric"
title <- paste("Smoothing Residuals of Robust Fit with Elev span=", par$span, sep="")
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")


## contour plot of smoothing residuals from spatial loess
rst <- rhread(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), par$loess
	)
)
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

trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.contour.loess.resid.0.02503", #with elev at the end of the plot name 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in c(sort(unique(result$year)))){
	for(j in levels(result$month)){
		b <- levelplot(
			resid.fit ~ lon * lat,
			data = result,
			subset = year == i & month == j, #& resid.fit < 5 & resid.fit > -5,
  		region = TRUE,
 		  at = c(seq(-4,-1,0.2),seq(-0.8,0.8,0.1), seq(1,4,0.2)),
 # 		cuts = 10,
  		col.regions = colorRampPalette(c("blue", "yellow","red")),
  	#	colorkey = list(at = seq(-1, 1, 0.1)),
  		xlab = "Longitude",
  		ylab = "Latitude",
  		main = paste("Smoothing Loess Residuals", j, i),
 			panel = function(x, y, z, ...) {
 				panel.levelplot(x,y,z,...)
    		panel.polygon(us.map$x,us.map$y, border = "red") 
 			}
 		)
 		print(b)
 	}
}
dev.off()


