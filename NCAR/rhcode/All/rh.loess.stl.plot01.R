###################################
##Plotting the stl components for each vertex of kd tree
##
##initialize the rhipe and setup the directories
source("~/Rhipe/rhinitial.R")
par <- list()
par$dataset <- "tmax"
par$loess <- "loess01"
par$machine <- "gacrux"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

rst <- rhread(
	file.path(
		rh.datadir, 
		par$dataset, 
		"spatial", 
		paste(par$loess, "stl", sep="." )
	)
)
result <- do.call("rbind", lapply(rst, "[[", 2))
result$factor <- factor(
	rep(
		rep(paste("Period", 1:9), c(rep(144,8),84)), 
		times = length(unique(result$fac))
	), 
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
spatial <- ddply(
	.data = result,
	.variables = "fac",
	.fun = summarise,
	 lon = unique(lon),
	 lat = unique(lat)
)
#######################################
##trend + seasonal component vs time
#######################################
trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/scatterplot_of_trend+seasonal_vertices_", 
		par$dataset, 
		".ps", 
		sep = ""
	), 
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
			layout = c(1,9),
			strip = FALSE,
			xlim = c(0, 143),
			scales = list(
				y = list(relation = 'same', alternating = TRUE), 
				x = list(at = seq(0, 143, by = 12), relation = 'same')
			),
			panel = function(x, y, subscripts, ...) {
				panel.abline(v = seq(0,145, by=12), color="lightgrey", lty = 3, lwd = 0.5)
        panel.xyplot(x, y, type="p", col=col[1], pch=16, cex=0.5, ...)
        sub <- subset(result, fac == i)
        panel.xyplot(
        	x = sub[subscripts,]$time2, 
        	y = (sub[subscripts,]$fitted - sub[subscripts,]$fc.remainder), 
        	type = "l", 
        	col = col[2], 
        	lwd = 1, ...
        )
			}
		)
		print(b)
  }
dev.off()


###########################################
##trend component and yearly mean vs time
###########################################
tmp.trend <- result[, c(
	!(names(result) %in% c("seasonal", "remainder", "data.weights", "factor"))
	)
]
dr <- ddply(
	.data = tmp.trend,
	.variables = c("fac","year"),
	.fun = summarise,
	mean = mean(fitted)
)
mm <- dr[rep(row.names(dr), each=12),]
tmp.trend <- tmp.trend[with(tmp.trend, order(fac, year, month)),]
tmp.trend <- cbind(tmp.trend, mean= mm$mean)

trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/scatterplot_of_lowfreq_vertex_", 
		par$dataset, 
		".ps", 
		sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
	b <- xyplot( mean ~ time | fac,
		data = tmp.trend,
		xlab = list(label = "Month", cex = 1.2),
		ylab = list(label = ylab, cex = 1.2),
		strip = strip.custom(
			par.strip.text = list(cex = 1), 
			factor.levels = paste(spatial$lat, spatial$lon, sep = ", ")
		),
		par.settings = list(
			layout.heights = list(strip = 1)
		),
		xlim = c(0, 1235),
		pch = 16,
		layout = c(10,1),
		aspect = "xy",
		key = list(
			text = list(label=c("low frequency component","yearly mean")), 
			lines = list(
				pch = c("","."), 
				cex = 4, 
				lwd = 1.5, 
				type = c("l","p"), 
				col = col[1:2]
			), 
			columns = 2
		),
		cex = 0.3,
		scales = list(
			y = list(relation = 'free'), 
			x = list(at=seq(0, 1236, by=600), relation = 'same')
		),
		prepanel = function(x, y, subscripts,...){
			v <- tmp.trend[subscripts,]
			ylim <- range(v$mean)
 			ans <- prepanel.default.xyplot(v$time, v$fc.trend, ...)
			ans$ylim <- range(ans$ylim, ylim)
			ans
		},
		panel = function(x, y, subscripts, ...) {
			panel.xyplot(x, y, type="p", col=col[2], ...)
			v <- tmp.trend[subscripts,]	
			panel.xyplot(v$time, v$fc.trend, type="l", col=col[1], ...)
		}
	)
    print(b)
dev.off()

trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/scatterplot_of_middlefreq_vertex_", 
		par$dataset, 
		".ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
	b <- xyplot( (mean-fc.trend) ~ time | fac,
		data = tmp.trend,
		xlab = list(label = "Month", cex = 1.2),
		ylab = list(label = ylab, cex = 1.2),
		strip = strip.custom(
			par.strip.text = list(cex = 1), 
			factor.levels = paste(spatial$lat, spatial$lon, sep = ", ")
		),
		par.settings = list(
			layout.heights = list(strip = 1)
		),
		xlim = c(0, 1235),
		layout = c(2,1),
		aspect = "xy",
		key = list(
			text = list(
				label = c(
					"middle frequency component",
					"yearly mean-low frequency component"
				)
			),
			lines = list(
				pch = c("","."), 
				cex = 4, 
				lwd = 1.5, 
				type =c ("l","p"), 
				col = col[1:2]
			), 
			columns = 2
		),
		scales = list(
			y = list(relation = 'free'), 
			x = list(at=seq(0, 1236, by=120), relation = 'same')
		),
		prepanel = function(x,y,subscripts,...){
			v <- tmp.trend[subscripts,]
			ylim <- range(v$mean-v$fc.trend)
			ans <- prepanel.default.xyplot(v$time, v$fc.second, ...)
			ans$ylim <- range(ans$ylim, ylim)
			ans
		},
		panel = function(x, y, subscripts, ...) {
			panel.xyplot(x,y,type="p", col=col[2], pch=16, cex=0.3,...)
			v <- tmp.trend[subscripts,]
			panel.xyplot(v$time, v$fc.second, type="l", col=col[1], ...)
			panel.abline(v=seq(0,1236, by=60), color="lightgrey", lty=3, lwd=0.5) 
		}
	)
	print(b)
dev.off()

rm(tmp.trend)

######################################################################
#QQ plot and time series plot of Pooled remainder of max temperature 
######################################################################
#Create the QQ plot of temperature for one station
trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/QQ_plot_of_remainder_vertex_", 
		par$dataset, 
		".ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
	a <- qqmath(~ fc.remainder | fac,
		data = result,
		distribution = qnorm,
		strip = strip.custom(
			par.strip.text = list(cex = 1),
			factor.levels = paste(spatial$lat, spatial$lon, sep = ", ")
		),
		aspect = "xy",
		pch = 16,
		cex = 0.3,
		layout = c(5,1),
		xlab = list(label = "Unit normal quantile", cex = 1.2),
		ylab = list(label = ylab, cex=1.2),
		prepanel = prepanel.qqmathline,
		panel = function(x, y,...) {
			panel.grid(lty = 3, lwd = 0.5, col = "black",...)
			panel.qqmathline(x, y = x)
			panel.qqmath(x, y,...)
		}
	)
	print(a)
dev.off()

tmp.remainder <- result
tmp.remainder$factor <- factor(
	tmp.remainder$factor, 
	levels = paste("Period", c(3,2,1,6,5,4,9,8,7))
)
trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/scatterplot_of_remainder_vertex_", 
		par$dataset, 
		".ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
	for(i in levels(tmp.remainder$fac)){
		b <- xyplot( fc.remainder ~ time2 | factor,
			data = subset(tmp.remainder, fac == i),
			xlab = list(label = "Month", cex = 1.2),
			ylab = list(label = ylab, cex = 1.2),
			main = list(label = paste("Vertex of ", "(",
				unique(subset(result, fac == i)$lat), ", ",
				unique(subset(result, fac == i)$lon), ")", sep="")
			),
			type = "p",
			layout = c(1,3),
			pch = 16,
			cex = 0.5,
			xlim = c(0, 143),
			strip = FALSE,
			scales = list(
				y = list(relation = 'same', alternating=TRUE), 
				x = list(at=seq(0, 143, by=12), relation='same')
			),
			panel = function(x,y,...) {
				panel.abline(h=0, v=seq(0,143, by=12), color="lightgrey", lty=3, lwd=0.5)
				panel.xyplot(x,y,...)
				panel.loess(x,y,degree=2,span=1/4, col=col[2], ...)
			}
		)
		print(b)
  }
dev.off()

trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/scatterplot_of_remainder2_vertex_", 
		par$dataset, 
		".ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
  for(i in levels(tmp.remainder$fac)){
		b <- xyplot( fc.remainder ~ time,
			data = subset(tmp.remainder, fac == i),
			xlab = list(label = "Month", cex = 1.2),
			ylab = list(label = paste(ylab), cex = 1.2),
			main = list(label = paste("Vertex of ", "(",
				unique(subset(result, fac == i)$lat), ", ",
				unique(subset(result, fac == i)$lon), ")", sep="")
			),
			type = "p",
			pch = 16,
			cex = 0.5,
			xlim = c(0, 1235),
			key = list(
				text = list(
					label=c(
						"remainder",
						"degree=2,span=0.15",
						"degree=1,span=0.35"
					)
				), 
				lines = list(
					pch = c(".", "", ""), 
					cex = 4, 
					lwd = 1.5, 
					type = c("p","l", "l"), 
					col = col[1:3]
				), 
				columns = 3
			),
			scales = list(
				y = list(relation = 'same', alternating=TRUE), 
				x = list(at=seq(0, 1235, by=120), relation='same')
			),
			panel = function(x, y,...) {
				panel.abline(h=0)
				panel.xyplot(x,y,...)
				panel.loess(
					x = x, 
					y = y, 
					degree = 2, 
					span = 0.15, 
					col = col[2], 
					evaluation = 100,...
				)
				panel.loess(
					x = x, 
					y = y, 
					degree = 1, 
					span = 0.35, 
					col = col[3], 
					evaluation = 100,...
				)
			}
		)
		print(b)
  }
dev.off()

#################################################
##Auto correlation ACF for the remainder
#################################################
ACF <- ddply(
	.data=tmp.remainder,
	.variables="fac",
	.fun= summarise,
	correlation = c(acf(fc.remainder, plot=FALSE)$acf),
	lag = c(acf(fc.remainder, plot=FALSE)$lag)
)

trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/acf_of_remainder_vertex_", 
		par$dataset, 
		".ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
  for(i in levels(ACF$fac)){
		b <- xyplot( correlation ~ lag,
			data = subset(ACF, fac == i & lag!=0),
			xlab = list(label = "Lag", cex = 1.2),
			ylab = list(label = "ACF", cex = 1.2),
			main = list(label = paste("Vertex of ", "(",
				unique(subset(result, fac == i)$lat), ", ",
				unique(subset(result, fac == i)$lon), ")", sep="")
			),
			type = "h",
			panel = function(x,y,...) {
				panel.abline(h=0)
				panel.xyplot(x,y,...)
			}
		)
		print(b)
  }
dev.off()

###############################################
##remainder of vertex conditional on month
###############################################
trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/scatterplot_of_remainder_vertex_conditional_month_", 
		par$dataset, 
		".ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
  for(i in levels(tmp.remainder$fac)){
		b <- xyplot( fc.remainder ~ year | month,
			data = subset(tmp.remainder, fac == i),
			xlab = list(label = "Year", cex = 1.2),
			ylab = list(label = ylab, cex = 1.2),
			main = list(label = paste("Vertex of ", "(",
				unique(subset(result, fac == i)$lat), ", ",
				unique(subset(result, fac == i)$lon), ")", sep="")
			),
			pch = 16,
			cex = 0.5,
		layout = c(12,1),
		strip = TRUE,
		scales = list(
			y = list(relation = 'same', alternating = TRUE), 
			x = list(tick.number = 10, relation = 'same')
		),
		panel = function(x,y,...){
			panel.abline(h=0, color="black", lty=1)
			panel.xyplot(x,y,...)
			panel.loess(x,y,span=3/4, degree=1, col=col[2],...)
		}
	)
	print(b)
  }
dev.off()

trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/scatterplot_of_remainder2_vertex_conditional_month_", 
		par$dataset, 
		".ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
  for(i in levels(tmp.remainder$fac)){
	b <- xyplot( fc.remainder ~ year | month,
			data = subset(tmp.remainder,fac == i),
			xlab = list(label = "Year", cex = 1.2),
			ylab = list(label = ylab, cex = 1.2),
			main = list(label = paste("Vertex of ", "(",
				unique(subset(result, fac == i)$lat), ", ",
				unique(subset(result, fac == i)$lon), ")", sep="")
			),
			type= "b",
			pch=16,
			cex=0.5,
			layout = c(2,6),
			strip = TRUE,
			scales = list(
				y = list(relation = 'same', alternating = TRUE), 
				x = list(tick.number = 10, relation = 'same')
			),
			panel = function(x,y,...) {
 				panel.abline(h=0, color="black", lty=1)
				panel.xyplot(x,y,...)
			}
		)
		print(b)
  }
dev.off()

trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/QQ_plot_of_remainder_vertex_conditional_month_", 
		par$dataset, 
		".ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
  for(i in levels(tmp.remainder$fac)){
		a <- qqmath(~ fc.remainder | month,
			data = subset(tmp.remainder, fac == i),
			distribution = qnorm,
			aspect = 1,
			pch = 16,
			cex = 0.3,
	 		layout = c(4,3),
			xlab = list(label="Unit normal quantile", cex=1.2),
			ylab = list(label = ylab, cex = 1.2),
			main = list(label = paste("Vertex of ", "(",
				unique(subset(result, fac == i)$lat), ", ",
				unique(subset(result, fac == i)$lon), ")", sep="")
			),
			prepanel = prepanel.qqmathline,
			panel = function(x, y,...) {
				panel.grid(lty=3, lwd=0.5, col="black",...)
				panel.qqmathline(x, y = x)
				panel.qqmath(x, y,...)
			}
		)
		print(a)
  }
dev.off()

rm(tmp.remainder)

###########################################
##Seasonal component conditional on month
###########################################
tmp.seasonal <- result[, c(
	!(names(result) %in% c("weights","trend","remainder"))
	)
]

trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/scatterplot_of_seasonal_vertex_conditional_month_", 
		par$dataset,
		".ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
  for(i in order.st){
		b <- xyplot( data.seasonal ~ year | month,
			data = subset(tmp.seasonal, fac == i),
			xlab = list(label = "Year", cex = 1.2),
			ylab = list(label = ylab, cex = 1.2),
			main = list(label = paste("Vertex of ", "(",
				unique(subset(result, fac == i)$lat), ", ",
				unique(subset(result, fac == i)$lon), ")", sep="")
			),
			type = c("p"),
			pch = 16,
			cex = 0.5,
			layout = c(12,1),
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

trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/scatterplot_of_seasonal2_vertex_conditional_month_", 
		par$dataset, 
		".ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
	b <- xyplot( data.seasonal ~ month | fac,
		data = subset(tmp.seasonal, year == 1895),
		xlab = list(label = "Month", cex = 1.2),
		ylab = list(label = ylab, cex = 1.2),
		type = "b",
		strip = strip.custom(
			par.strip.text = list(cex = 1),
			factor.levels = paste(spatial$lat, spatial$lon, sep = ", ")
		),
		pch = 16,
		aspect = "xy",
		cex = 0.5,
		layout = c(5,2),
		scales = list(
			y = list(relation = 'same', alternating=TRUE), 
			x = list(at=c(1, 3, 5, 7, 9, 11), relation='same')
		),
		panel = function(x,y,...) {
			panel.xyplot(x,y,...)
		}
	)
	print(b)
dev.off()

rm(tmp.seasonal)