###################################
##Plotting the stl components for each vertex of kd tree
##
##initialize the rhipe and setup the directories
source("~/Rhipe/rhinitial.R")
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
library(lattice)
library(plyr)
par <- list()
par$dataset <- "tmax"
par$loess <- "loess01"

rst <- rhread(file.path(rh.datadir, par$dataset, "spatial", paste(par$loess, "stl", sep="." )))
result <- do.call("rbind", lapply(rst, "[[", 2))
result$factor <- factor(rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times=66), levels=paste("Period", c(9:1)))
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

#######################################
##trend + seasonal component vs time
#######################################
trellis.device(postscript, file = paste(local.output, "/scatterplot_of_trend+seasonal_vertices_", par$dataset, ".ps", sep = ""), color=TRUE, paper="legal")
  for(i in unique(result$fac)) {
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
        	panel.xyplot(sub[subscripts,]$time2, 
				(sub[subscripts,]$trend + sub[subscripts,]$seasonal), 
				type = "l", 
				col = col[2], 
				lwd=1, ...)
        }
	)
	print(b)
  }
dev.off()


###########################################
##trend component and yearly mean vs time
###########################################
tmp.trend <- result[, c(!(names(result) %in% c("seasonal", "remainder", "data.weights", "factor")))]

dr <- ddply(.data = tmp.trend,
        .variables = c("fac","year"),
        .fun = summarise,
         mean = mean(fitted)
)
mm <- dr[rep(row.names(dr), each=12),]
tmp.trend <- tmp.trend[with(tmp.trend, order(fac, year, month)),]
tmp.trend <- cbind(tmp.trend, mean= mm$mean)

order <- ddply(.data = tmp.trend,
                .variables = "fac",
                .fun = summarise,
                 mean = mean(fitted)
)
order.st <- as.character(order[order(order$mean, decreasing=TRUE), ]$fac)
tmp.trend$fac <- factor(tmp.trend$fac, levels=order.st)

trellis.device(postscript, file = paste(local.output, "/scatterplot_of_trend_vertex_", par$dataset, ".ps", sep = ""), color=TRUE, paper="legal")
	b <- xyplot( mean ~ time | fac,
		data = tmp.trend,
		xlab = list(label = "Month", cex = 1.2),
		ylab = list(label = ylab, cex = 1.2),
		strip = strip.custom(
			par.strip.text = list(cex = 1), 
			factor.levels = paste(tmp.trend$lat, tmp.trend$lon, sep = "")
		),
		par.settings = list(layout.heights = list(strip = 1)),
		xlim = c(0, 1235),
		pch = 16,
		layout = c(5,2),
		aspect = "xy",
		key = list(
			text = list(label=c("low frequency component","yearly mean")), 
			lines = list(pch=c("","."), cex=4, lwd=1.5, type=c("l","p"), col=col[1:2]), 
			columns = 2
		),
		cex = 0.3,
		scales = list(y = list(relation = 'free'), x=list(at=seq(0, 1236, by=600), relation = 'same')),
		prepanel = function(x, y, subscripts,...){
			v <- tmp.trend[subscripts,]
			ylim <- range(v$mean)
 			ans <- prepanel.default.xyplot(v$time, v$trend, ...)
			ans$ylim <- range(ans$ylim, ylim)
			ans
		},
		panel = function(x, y, subscripts, ...) {
			panel.xyplot(x, y, type="p", col=col[2], ...)
			v <- tmp.trend[subscripts,]	
			panel.xyplot(v$time, v$trend, type="l", col=col[1], ...)
		}
	)
    print(b)
dev.off()

