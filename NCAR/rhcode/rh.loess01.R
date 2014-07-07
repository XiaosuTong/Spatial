###############################################
##
##get the loess fit at each vertix of the kd tree using the whole 
##dataset. The output key is the vertix id, and value is the 1236 
##fitted value for that vertix.
#


##initialize the rhipe and setup the directories
source("~/Rhipe/rhinitial.R")
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/rhcode/my.loess.R")
##set up the parameters
par <- list()
par$dataset <- "tmax"
par$N <- 1236
par$span <- 0.2
par$degree <- 2
##construct the kd tree for given parameters of span and degree
##station information is in USpinfo and UStinfo object.
if(par$dataset == "precip"){
	load(file.path(local.datadir, "USpinfo.RData"))
    info <- USpinfo
}else{
	load(file.path(local.datadir, "UStinfo.RData"))
    info <- UStinfo
}
rm(list=grep("US", ls(), value=T))
source("~/Projects/Spatial/NCAR/code/spatial/kdtree.R")
kd <- kdtree(info[, c("station.id","lon","lat")], alpha = par$span/5)
rhsave(list=("kd"), file=file.path(rh.datadir, par$dataset, "Rdata", paste("kd", ".RData", sep="")))

##loess fit at each vertices
job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
	i <- ceiling(map.keys[[r]]/12)
	j <- map.keys[[r]] - (i - 1)*12
	month <- c("Jan","Feb","Mar","Apr","May","June","July","Aug", "Sep", "Oct", "Nov", "Dec")
	m <- month[j]
	y <- i + 1894
	v <- subset(get(par$dataset), year == y & month == m)[, c("station.id", "elev", "lon", "lat", par$dataset)]
	lo.fit <- my.loess( get(par$dataset) ~ lon + lat, 
		data    = v, 
		degree  = par$degree, 
		span    = par$span, 
		control = loess.control(surface = "direct")
	)
	value <- predict(lo.fit, data.frame(lon = kd$vertix$lon, lat = kd$vertix$lat))
	value <- cbind(kd$vertix, value)
	names(value) <- c("lon","lat","fitted")
	value$year <- rep(y, nrow(value))
	value$month <- rep(m, nrow(value))
	value$fac <- seq_len(nrow(value))
	lapply(1:nrow(value), function(k){
		key <- value$fac[k]
		rhcollect(key, value[k, ])
	})
  })
})
job$reduce <- expression(
	pre={
		combined <- data.frame()
	},
	reduce={
		combined <- rbind(combined, do.call(rbind, reduce.values))
	},
	post={
		if(nrow(combined) == 1236) {
			rhcounter("BUG", "1236", 1)
		}
		stopifnot(nrow(combined) == 1236)
		combined$month <- factor(combined$month, 
			levels = c("Jan", "Feb", "Mar", "Apr", "May", "June", 
    		"July", "Aug", "Sep", "Oct", "Nov", "Dec")
		)
		combined <- combined[with(combined, order(year, month)), ]
		combined$time <- 0:1235
		rhcollect(reduce.key, combined)
	}	
)
job$setup <- expression(
	map = {
		load("kd.RData")
	    load(paste(par$dataset, "RData", sep="."))
	}
)
job$shared <- c(
	file.path(rh.datadir, par$dataset, "Rdata", "kd.RData"),
	file.path(rh.datadir, par$dataset, "Rdata", paste(par$dataset, "RData", sep="."))
)
job$parameters <- list(par = par, my.loess = my.loess, my.simple = my.simple)
job$input <- c(par$N, 242) 
job$output <- rhfmt(file.path(rh.datadir, par$dataset, "spatial", "loess01"), type="sequence")
job$mapred <- list(mapred.reduce.tasks = 66, rhipe_reduce_buff_size=10000)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "loess01")
job$readback <- FALSE
job$noeval <- TRUE
job.mr <- do.call("rhwatch", job)
z <- rhex(job.mr, async = TRUE)


###################################
##Plotting the time series for each vertex of kd tree
##
rst <- rhread("/wsc/tongx/Spatial/tmp/tmax/spatial/loess01")
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
trellis.device(postscript, file = paste(local.output, "/scatterplot_of_vertices_", par$dataset, ".ps", sep = ""), color=TRUE, paper="legal")
  for(i in unique(result$fac)) {
	b <- xyplot( fitted ~ time2 | factor,
    	data = subset(result, fac == i),
        xlab = list(label = "Month"),
        ylab = list(label = ylab),
		main = list(label = paste("Vertex of ", "(", 
			unique(subset(result, fac == i)$lat), ", ",
			unique(subset(result, fac == i)$lon), ")", sep="")
		),
		type = "b",
		pch=16,
		cex=0.5,
		layout = c(1,9),
#       strip.left = TRUE,
		strip = FALSE,
#       aspect= 0.1,
		grib = TRUE,
		xlim = c(0, 143),
#       scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
		scales = list(y = list(relation = 'same', alternating=TRUE), x=list(at=seq(0,143,by=12), relation='same')),
		panel = function(...) {
			panel.abline( v=seq(0,145,by=12), color="lightgrey", lty=3, lwd=0.5)
			panel.xyplot(...)
		}
	)
	print(b)
  }
dev.off()

library(maps)
us.map <- map('state', plot = FALSE, fill = TRUE)
trellis.device(postscript, file = paste(local.output, "/map_vertices_", par$dataset, ".ps", sep = ""), color=TRUE, paper="legal")
a <- xyplot( lat ~ lon,
	data  = result,
	xlab  = list(label = "Longitude"),
	ylab  = list(label = "Latitude"),
	main  = list(label = "Vertices of K-D Tree"),
	type  = "p",
	cex   = 0.5,
	col   = "red",
	pch   = 16,
	panel = function(x, y, ...) {
		panel.xyplot(x, y, ...)
	    panel.polygon(us.map$x,us.map$y)	
	}
)
print(a)
dev.off()

trellis.device(postscript, file = paste(local.output, "/QQ_plot_of_month_", par$dataset, ".ps", sep = ""), color=TRUE, paper="legal")
  for(i in unique(result$fac)){
    a <- qqmath(~ fitted | month,
        data = subset(result, fac == i),
        distribution = qnorm,
        aspect = "xy",
        layout = c(12,1),
        pch  = 16,
        cex  = 0.5,
        main = list(label = paste("Vertex of ", "(",
            unique(subset(result, fac == i)$lat), ", ",
            unique(subset(result, fac == i)$lon), ")", sep="")
        ),
		xlab = list(label = "Unit normal quantile"),
		ylab = list(label = ylab, cex=1.2),
#       scales = list(x = list(cex=1.5), y = list(cex=1.5)),
        prepanel = prepanel.qqmathline,
        panel = function(x, y,...) {
                panel.grid()
                panel.qqmathline(x, y=x)
                panel.qqmath(x, y, ...)
        }
    )
    print(a)
  }
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_precipitation_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
    b <- xyplot( precip ~ year | month,
             data = subset(tmp,station.id==i),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.2),
         main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = c("p"),
             pch=16,
             cex=0.5,
         layout = c(4,3),
         strip = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
         scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 40,by=10), v=seq(1900,2000,by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
          panel.loess(x,y,span=2/3,degree=1, col=col[2],...)
             }
        )
        print(b)
     }
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_precipitation_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
    b <- xyplot( precip ~ year | month,
             data = subset(tmp,station.id==i),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.2),
         main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = c("p"),
             pch=16,
             cex=0.5,
         layout = c(4,3),
         strip = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
         scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 40,by=10), v=seq(1900,2000,by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
          panel.loess(x,y,span=2/3,degree=1, col=col[2],...)
             }
        )
        print(b)
     }
dev.off()

