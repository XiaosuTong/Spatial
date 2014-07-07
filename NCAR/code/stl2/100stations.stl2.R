############################################################################
##Max temperature, 100 stations, seasonal+trend with raw data by stl2
##########################################################################

#load the data
library(lattice)
library(plyr)
library(stl2)
library(maps)
outputdir <- "~/Projects/Spatial/NCAR/output/"
datadir <- "~/Projects/Spatial/NCAR/RData/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))
load(paste(datadir,"stations.RData", sep=""))

dataset <- "tmax"
#parameter <- list(sw="periodic", sd=1, tw=1141, td=1, inner=10, outer=0, flag=FALSE)
#parameter <- list(sw="periodic", sd=1, tw=1855, td=1, fcw=c(1855,121), fcd=c(1,2), inner=10, outer=0, flag=TRUE)
parameter <- list(sw="periodic", sd=1, tw=1855, td=1, fcw=c(1855,241), fcd=c(1,1), inner=10, outer=0, flag=TRUE)

if(dataset == "tmax"){
ylab <- "Maximum Temperature (degrees centigrade)"
}else if(dataset == "tmin"){
        ylab <- "Minimum Temperature (degrees centigrade)"
}else {
        ylab <- "Precipitation (millimeters)"
}
if(dataset %in% c("tmax", "tmin")){
        data <- UStemp
        datainfo <- UStinfo
}else{
        data <- USppt
        ddatainfo <- USpinfo
}
rm(list=grep("US", ls(), value=TRUE))
#find the 100 stations with largest observation number for max temperature
stations <- get(paste("stations", dataset, sep="."))
tmp <- data[with(data, station.id %in% stations),]
if(any(with(tmp, is.na(get(dataset))) == TRUE)) stop("The first 100 stations have NA")

month <- tmp$month
levels(month) <- c(1:12)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
column <- c("station.id", "year","month","date", dataset)
tmp <- tmp[with(tmp, order(station.id, date)), column]
names(tmp)[dim(tmp)[2]] <- "response"
tmp$factor <- factor(rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times=100), levels=paste("Period", c(9:1)))
tmp$time <- c(rep(0:143,8), 0:83)

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

if(parameter$flag){
dr <- ddply(tmp, "station.id", function(r){
	do.call("cbind", stl2(r$response, r$date,
	n.p=12,
	s.window = parameter$sw,
	s.degree = parameter$sd,
	t.window = parameter$tw,
	t.degree = parameter$td,
	fc.window = parameter$fcw,
	fc.degree = parameter$fcd,
	inner = parameter$inner,
	outer = parameter$outer)[c("data","fc")])
})
tmp <- cbind(tmp, dr[, c(!(names(dr) %in% c("station.id", "data.raw", "data.sub.labels")))])
names(tmp)[grep("fc.fc", names(tmp))] <- c("fc.trend", "fc.second")
}else{
dr <- ddply(tmp, "station.id", function(r){
	stl2(r$response, r$date, 
	n.p = 12, 
	s.window = parameter$sw, 
	s.degree = parameter$sd, 
	t.window = parameter$tw, 
	t.degree = parameter$td, 
	inner = parameter$inner, 
	outer = parameter$outer)$data
})
tmp <- cbind(tmp, dr[, c(!(names(dr) %in% c("station.id", "raw", "sub.labels")))])
}

##################################
##trend+seasonal time series plot
##################################
trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_with_stl2_trend+seasonal_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
	for(i in stations){
    	b <- xyplot( response ~ time | factor,
        	data = subset(tmp, station.id==i),
            xlab = list(label = "Month", cex = 1.2),
            ylab = list(label = paste("Station", i, ylab, sep=" "), cex = 1.2),
            layout = c(1,9),
			aspect= 0.06,
			strip = FALSE,
            grib = TRUE,
            xlim = c(0, 143),
#           scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
            scales = list(y = list(relation = 'same', tick.number=4, alternating=TRUE), x=list(at=seq(0, 143, by=12), relation='same')),
            panel = function(x,y,subscripts,...) {
                 panel.abline(v=seq(0,145, by=12), color="lightgrey", lty=3, lwd=0.5)
                 panel.xyplot(x, y, type="p", col=col[1], pch=16, cex=0.5, ...)
				 sub <- subset(tmp, station.id==i)
				 panel.xyplot(sub[subscripts,]$time, (sub[subscripts,]$trend+sub[subscripts,]$seasonal), type="l", col=col[2], lwd=1, ...)            
            }
        )
        print(b)
   }
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_with_stl2_seasonal+fc_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
    for(i in stations){
        b <- xyplot( response ~ time | factor,
            data = subset(tmp, station.id==i),
            xlab = list(label = "Month", cex = 1.2),
            ylab = list(label = paste("Station", i, ylab, sep=" "), cex = 1.2),
            layout = c(1,9),
            aspect= 0.06,
            strip = FALSE,
            grib = TRUE,
            xlim = c(0, 143),
#           scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
            scales = list(y = list(relation = 'same', tick.number=4, alternating=TRUE), x=list(at=seq(0, 143, by=12), relation='same')),
            panel = function(x,y,subscripts,...) {
                 panel.abline(v=seq(0,145, by=12), color="lightgrey", lty=3, lwd=0.5)
                 panel.xyplot(x, y, type="p", col=col[1], pch=16, cex=0.5, ...)
                 sub <- subset(tmp, station.id==i)
                 panel.xyplot(sub[subscripts,]$time, (sub[subscripts,]$data.seasonal+sub[subscripts,]$fc.trend+sub[subscripts,]$fc.second), type="l", col=col[2], lwd=1, ...)
            }
        )
        print(b)
   }
dev.off()

##################################################################################
#QQ plot and time series plot of Pooled remainder of max temperature for 100 stations
##################################################################################
#Create the QQ plot of temperature for one station
trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_stl2_remainder_", dataset ,"_of_100_stations", ".ps", sep = ""), color=TRUE, paper="legal")
	a <- qqmath(~ fc.remainder | unique(station.id),
		data = tmp,
		distribution = qnorm,
		aspect = 1,
		pch = 16,
		cex = 0.3,
		layout = c(6,3),
        xlab = list(label="Unit normal quantile", cex=1.2),
        ylab = list(label=ylab, cex=1.2),
#       scales = list(x = list(cex=1.5), y = list(cex=1.5)),
        prepanel = prepanel.qqmathline,
        panel = function(x, y,...) {
        	panel.grid(lty=3, lwd=0.5, col="black",...)
            panel.qqmathline(x, y=x)
            panel.qqmath(x, y,...)
        }
	)
    print(a)
dev.off()

tmp.remainder <- tmp
tmp.remainder$time1 <- rep(0:1235,100)
tmp.remainder$factor <- factor(tmp$factor, levels=paste("Period", c(3,2,1,6,5,4,9,8,7)))
trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_with_stl2_remainder_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        b <- xyplot( fc.remainder ~ time | factor,
             data = subset(tmp.remainder, station.id==i),
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = paste("Station", i, ylab), cex = 1.2),
             type = "p",
             layout = c(1,3),
             pch = 16,
             cex = 0.5,
	    	 xlim = c(0, 143),
             strip = FALSE,
             grib = TRUE,
#            scales = list(y = list(relation='free', cex=1.5), x=list(relation='free',format="%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(at=seq(0, 143, by=12), relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=0, v=seq(0,143, by=12), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
		  	 	  panel.loess(x,y,degree=2,span=1/4, col=col[2], ...)
             }
        )
        print(b)
   }
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_with_stl2_remainder2_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        b <- xyplot( fc.remainder ~ time1,
             data = subset(tmp.remainder, station.id==i),
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = paste("Station", i, ylab), cex = 1.2),
             type = "p",
             pch = 16,
             cex = 0.5,
             xlim = c(0, 1235),
			 key = list(text=list(label=c("remainder","degree=2,span=0.15","degree=1,span=0.35")), lines=list(pch=c(".", "", ""), cex=4, lwd=1.5, type=c("p","l", "l"), col=col[1:3]), columns=3),
#           scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(at=seq(0, 1235, by=120), relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=0)
                  panel.xyplot(x,y,...)
                  panel.loess(x,y,degree=2,span=0.15, col=col[2], evaluation=100,...)
		  		  panel.loess(x,y,degree=1,span=0.35, col=col[3], evaluation=100,...)
             }
        )
        print(b)
   }
dev.off()


#################################################
##Auto correlation ACF for the remainder
#################################################
ACF <- ddply(.data=tmp.remainder,
	.variables="station.id",
	.fun= summarise,
	 correlation = c(acf(fc.remainder, plot=FALSE)$acf),
	 lag = c(acf(fc.remainder, plot=FALSE)$lag) 
)

trellis.device(postscript, file = paste(outputdir, "acf_of_", dataset, "_with_stl2_remainder_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        b <- xyplot( correlation ~ lag,
             data = subset(ACF, station.id==i & lag!=0),
             xlab = list(label = "Lag", cex = 1.2),
             ylab = list(label = paste("Station", i, "ACF"), cex = 1.2),
#             main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = "h",
             panel = function(x,y,...) {
                  panel.abline(h=0)
                  panel.xyplot(x,y,...)
             }
        )
        print(b)
   }
dev.off()

rm(tmp.remainder)
##########################################################################
##time series plot of trend component and yearly average of raw from stl2
##########################################################################
tmp.trend <- tmp[, c(!(names(tmp) %in% c("data.seasonal", "data.remainder", "data.trend", "data.weights", "factor")))]

dr <- ddply(.data = tmp.trend,
		.variables = c("station.id","year"),
		.fun = summarise,
		 mean = mean(response)
)
mm <- dr[rep(row.names(dr), each=12),]
tmp.trend <- cbind(tmp.trend, mean= mm$mean)

order <- ddply(.data = tmp.trend,
                .variables = "station.id",
                .fun = summarise,
                 mean = mean(response)
)
order.st <- as.character(order[order(order$mean, decreasing=TRUE), ]$station.id)
tmp.trend$station.id <- factor(tmp.trend$station.id, levels=order.st)
tmp.trend$time <- rep(0:1235,100)

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_with_stl2_fc.trend_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( mean ~ time | station.id, 
             data = tmp.trend,
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = ylab, cex = 1.2),
#            strip = strip.custom(par.strip.text= list(cex = 1.5)),
#            par.settings = list(layout.heights = list(strip = 1.5)),
	     	 xlim = c(0, 1235),
             pch = 16,
			 layout = c(10,1),
	     	 aspect="xy",
             key=list(text=list(label=c("low frequency component","yearly mean")), lines=list(pch=c("","."), cex=4, lwd=1.5, type=c("l","p"), col=col[1:2]), columns=2),
             cex = 0.3,
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'same',format = "%b %Y", tick.number=8), cex=1.2),
             scales = list(y = list(relation = 'free'), x=list(at=seq(0, 1236, by=600), relation = 'same')),
			 prepanel = function(x,y,subscripts,...){
					v <- tmp.trend[subscripts,]
					ylim <- range(v$mean)
					ans <- prepanel.default.xyplot(v$time, v$fc.trend, ...)
					ans$ylim <- range(ans$ylim, ylim)
					ans
			 },
             panel = function(x, y, subscripts, ...) {
                  panel.xyplot(x,y,type="p", col=col[2], ...)
		  		  v <- tmp.trend[subscripts,]
		  		  panel.xyplot(v$time, v$fc.trend, type="l", col=col[1], ...)
             }
        )
	print(b)
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_with_stl2_fc.second_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( (mean-fc.trend) ~ time | station.id,
             data = tmp.trend,
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = ylab, cex = 1.2),
#            strip = strip.custom(par.strip.text= list(cex = 1.5)),
#            par.settings = list(layout.heights = list(strip = 1.5)),
             xlim = c(0, 1235),
             layout = c(1,2),
             aspect="xy",
             key=list(text=list(label=c("middle frequency component","yearly mean-low frequency component")), lines=list(pch=c("","."), cex=4, lwd=1.5, type=c("l","p"), col=col[1:2]), columns=2),
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'same',format = "%b %Y", tick.number=8), cex=1.2),
             scales = list(y = list(relation = 'free'), x=list(at=seq(0, 1236, by=60), relation = 'same')),
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

######################################################
##Trend + remainder component and loess from stl2
######################################################
tmp.fc <- tmp[, c(!(names(tmp) %in% c("weights","trend", "remainder")))]
tmp.fc$tr <- tmp.fc$response - tmp.fc$seasonal
tmp.fc$time <- rep(0:1235,100)

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_after_seasonal_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tr ~ time | station.id,
             data = tmp.fc,
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = ylab, cex = 1.2),
			 layout = c(1,1),
#             aspect="xy",
             xlim = c(0, 1235),	
             key=list(type="l", text=list(label=c("after seasonal","loess smoothing")),  lines=list(lty = c(1,6), lwd=1.5, col=col[1:2]), columns=2),
             scales = list(y = list(relation = 'free'), x=list(at=seq(0, 1236, by=300), relation = 'same')),
             prepanel = function(x,y,...) prepanel.loess(x,y,span=900/1236, degree=2, ...),
             panel = function(x,y,...) {
                  panel.loess(x,y,degree=2, span=900/1236, col=col[2], type = "l", lty=1, ...)
                  panel.xyplot(x,y,col=col[1], type= "p", cex=0.3, ...)
             }
        )
        print(b)
dev.off()

rm(tmp.fc)
#####################################################
##time series plot of trend component loess from stl2
#####################################################
trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_with_stl2_fc.trend_loess_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( fc.trend ~ time | station.id,
             data = tmp.trend[order(tmp.trend$time),],
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = ylab, cex = 1.2),
#            strip = strip.custom(par.strip.text= list(cex = 1.5)),
#            par.settings = list(layout.heights = list(strip = 1.5)),
             aspect="xy",
	     	 layout = c(5,2),
	     	 xlim = c(0, 1235),
             key=list(type="l", text=list(label=c("loess smoothing","trend component")),  lines=list(lty = c(1,6), lwd=1.5, col=col[1:2]), columns=2),
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'same',format = "%b %Y", tick.number=8), cex=1.2),
             scales = list(y = list(relation = 'free'), x=list(at=seq(0, 1236, by=300), relation = 'same')),
	     	 prepanel = function(x,y,...) prepanel.loess(x,y,span=3/4, degree=1, ...),
             panel = function(x,y,...) {
#                 panel.abline(v=seq(0,1236, by=300), color="lightgrey", lty=3, lwd=0.5)
                  panel.loess(x,y,degree=1, span=3/4, col=col[1], type = "l", lty=1, ...)
		  		  panel.xyplot(x,y,col=col[2], type= "l", lty = 6, ...)
             }
        )
	print(b)
dev.off()


##############################
##
##############################
dr <- ddply(.data = tmp,
        .variables = c("station.id", "year"), 
        .fun = summarise,
         mean = mean(response),
		 range = max(response) - min(response)
)
order <- ddply(.data = tmp,
        .variables = "station.id",
        .fun = summarise,
         mean = mean(response)
)
order.st <- as.character(order[order(order$mean, decreasing=TRUE), ]$station.id)
dr$station.id <- factor(dr$station.id, levels=order.st)

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_range_vs_mean_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
	a <- xyplot(mean ~ range | station.id, 
		data = dr,
		pch = 16,
        xlab = list(label = "Yearly Range", cex = 1.2),
        ylab = list(label = "Yearly Mean", cex = 1.2),
		cex = 0.5,
		layout =c (10,1),
		scales = list(y = list(relation = 'same'), x=list(relation = 'same')),

	)
	print(a)
dev.off()


##########################################################################################
##QQ plot of the slopes of loess of trend component, ordered stations by average of slopes
#		 lat = lat[1:1235],
#		 slope = diff(loess)
#)
#
##Calculate the mean of slopes for each station, and then order the stations by mean slope.
#mean.slope <- ddply(.data = slope.trend,
#                .variables = "station.id",
#                .fun = summarise,
#		 elev = unique(elev),
#		 lon = unique(lon),
#		 lat = unique(lat),
#                 mean = mean(slope)
#)
#mean.slope <- mean.slope[order(mean.slope$mean),]
#station.or <- as.character(mean.slope$station.id)
#
##Attend to add a small map in the corner of the trend loess plot
#us.map <- map('state', plot = FALSE, fill = TRUE)
#
#trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_tmax_slope_of_trend_loess",".ps", sep = ""), color=TRUE, paper="legal")
#    for(i in station.or){
#        a <- qqmath(~ slope,
#                data = subset(slope.trend, station.id == i),
#                distribution = qunif,
#                aspect = 1,
#                pch = 16,
#                cex = 0.5,
##                main = list(label= paste("Station ", i, sep=""), cex=2),
#                xlab = list(label="f-value", cex=1.2),
#                ylab = list(label="Maximum temperature(degrees centigrade)", cex=1.2),
##                scales = list(x = list(cex=1.5), y = list(cex=1.5)),
#                prepanel = prepanel.qqmathline,
#                panel = function(x, y,...) {
#                        #panel.grid(lty=3, lwd=0.5, col="black",...)
#                        panel.qqmath(x, y,...)
#                }
##        )
#        b <- xyplot(lat ~ lon,
#                data = subset(tmp.trend, station.id==i),
#                xlab = NULL,
#                ylab = NULL,
#                pch = 16,
#                cex = 0.4,
#                xlim = c(-127,-65),
#                ylim = c(23,50),
#                col = "red",
##               scales = list(x = list(draw=FALSE), y = list(draw=FALSE)),
##               strip = strip.custom(par.strip.text= list(cex = 1.5)),
##               par.settings = list(layout.heights = list(strip = 1.5)),
#                panel = function(x,y,...) {
#                        panel.polygon(us.map$x,us.map$y,lwd=0.2)
#                        panel.xyplot(x,y,...)
#                }
#        )
#        plot.new()
#        title(paste("Station ", i, sep=""), cex=1.5)
#        print(a, pos=c(0.4,0,1,1), newpage=FALSE, more=TRUE)
#        print(b, pos=c(0,0.5,0.4,1), more=FALSE)
#    }
#dev.off()
#
##scatter plot of average slops vs spatial factor
#trellis.device(postscript, file = paste(outputdir, "scatter_plot_of_average_slopes_tmax",".ps", sep = ""), color=TRUE, paper="legal")
#        c <- xyplot(mean ~ log2(elev),
#                data = mean.slope,
#                xlab = list(label = "Elevation (meter)", cex=1.2),
#                ylab = list(label = "Average slope of trend", cex=1.2),
#                pch = 16,
#                cex = 0.7,
##               scales = list(x = list(draw=FALSE), y = list(draw=FALSE)),
##               strip = strip.custom(par.strip.text= list(cex = 1.5)),
##               par.settings = list(layout.heights = list(strip = 1.5)),
#                panel = function(x,y,...) {
#                        panel.xyplot(x,y,...)
#                }
#        )
#	print(c)
#dev.off()
