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


job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		year <- as.numeric(map.keys[[r]][1])
		month <- map.keys[[r]][2]
#		b <- contourplot(
#			resid.fit ~ lon * lat,
#			data = rst[[1]][[2]],
#  		region = TRUE,
# 		  at = c(seq(-4,-1,0.5),seq(-0.8,0.8,0.2), seq(1,4,0.5)),
#  		cuts = 20,
#  	#	colorkey = list(at = seq(-3, 3, 0.2)),
#  		xlab = "Longitude",
#  		ylab = "Latitude",
#  		#main = paste(year, month, title),
# 			panel = function(x, y, z, ...) {
# 				panel.contourplot(x,y,z,...)
#    		panel.polygon(us.map$x,us.map$y, border = "red") 
# 			}
# 		)
# 		od <- (year - 1950)*12 + which(month.abb==substr(month,1,3))
		b <- xyplot(1:10~1:10)
 		rhcollect(1, b)
	})
})
job$reduce <- expression(
	pre = {
		combine <- list()
	},
	reduce = {
		combine <- reduce.values
	},
	post = {
		trellis.device(
			postscript, 
			file = paste(
				"./tmp", "/a1950.contour.loess.resid.", #with elev at the end of the plot name 
				par$dataset, ".ps", sep = ""
			), 
			color = TRUE,
			paper = "legal"
		)
   	for(i in 1:length(combine)){print(combine[[i]])}
		dev.off()
	}
)
job$setup <- expression(
	map = {
		library(lattice)
	  library(maps, lib.loc = lib.loc)
	},
	reduce = {
		library(lattice)
	}
)
job$parameters <- list(
	par = par,
	title = title,
	us.map = us.map,
	col = col
)
job$input <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), par$loess
	), 
	type = "sequence"
)
job$output <- rhfmt(
	file.path(
		rh.output, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), par$loess
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 1,
	rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 5
job$copyFiles = TRUE
job$jobname <- file.path(
	rh.output, par$dataset, "spatial", "a1950",
	par$family, paste("sp", par$span, sep=""), par$loess
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)