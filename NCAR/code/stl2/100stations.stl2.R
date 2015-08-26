############################################################################
##Max temperature, 100 stations, seasonal+trend with raw data by stl2
##########################################################################

#load the data
library(lattice)
library(plyr)
library(maps)

##################################
##trend+seasonal time series plot
##################################
fitRaw <- function(data=rst, outputdir, target="tmax", size = "letter", test = TRUE){

  data$factor <- factor(
    x = rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times=100),
    levels = paste("Period", c(9:1))
  )
  data$time <- c(rep(0:143,8), 0:83) 

  stations <- unique(data$station.id)

  if(test) {
    stations <- stations[1]
  }

  if (target == "tmax") {
    ylab <- "Maximum Temperature (degrees centigrade)"
  } else if (target == "tmin") {
    ylab <- "Minimum Temperature (degrees centigrade)"
  } else {
    ylab <- "Precipitation (millimeters)"
  }

  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("fitted", "100stations", target, "ps", sep=".")),
    color = TRUE, 
    paper = size
  )
    for(i in stations){
      sub <- subset(data, station.id==i)
      b <- xyplot( resp ~ time | factor,
        , data = sub
        , xlab = list(label = "Month", cex = 1.5)
        , ylab = list(label = ylab, cex = 1.5)
        , main = list(label=paste("Station ", i, sep=""), cex=1)
        , layout = c(1,9)
        , aspect= "xy"
        , strip = FALSE,
        , xlim = c(0, 143)
        , key=list(
          text = list(label=c("raw","fitted")), 
          lines = list(pch=16, cex=0.7, lwd=1.5, type=c("p","l"), col=col[1:2]),
          columns=2
        )
        , scales = list(
            y = list(tick.number=4), 
            x = list(at=seq(0, 143, by=12), relation='same')
          )
        , panel = function(x,y,subscripts,...) {
            panel.abline(v=seq(0,145, by=12), color="lightgrey", lty=3, lwd=0.5)
            panel.xyplot(x, y, type="p", col=col[1], pch=16, cex=0.5, ...)
            if (!any(grepl("fc", names(rst)))) {
              panel.xyplot(sub[subscripts,]$time, (sub[subscripts,]$trend+sub[subscripts,]$seasonal), type="l", col=col[2], lwd=1, ...)            
            } else {
              panel.xyplot(sub[subscripts,]$time, (sub[subscripts,]$data.seasonal+sub[subscripts,]$fc.trend+sub[subscripts,]$fc.second), type="l", col=col[2], lwd=1, ...)
            }
          }
      )
      print(b)
    }
  dev.off()

}

##################################################################################
#QQ plot and time series plot of Pooled remainder of max temperature for 100 stations
##################################################################################
#Create the QQ plot of temperature for one station
remainder <- function(data=rst, outputdir, target="tmax", size = "letter", test = TRUE) {

  stations <- unique(data$station.id)

  if(test) {
    stations <- stations[1]
  }

  if (target == "tmax") {
    ylab <- "Maximum Temperature (degrees centigrade)"
  } else if (target == "tmin") {
    ylab <- "Minimum Temperature (degrees centigrade)"
  } else {
    ylab <- "Precipitation (millimeters)"
  }
  
  idx <- grep("remainder", names(data))
  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("remainder", "100stations", target, "ps", sep=".")), 
    color = TRUE, 
    paper = size
  )
    a <- qqmath( ~ data[, idx] | unique(station.id),
      , data = data, 
      , distribution = qnorm,
      , aspect = 1,
      , pch = 16,
      , cex = 0.3,
      , layout = c(4,3),
      , xlab = list(label="Unit normal quantile", cex=1.5),
      , ylab = list(label=ylab, cex=1.5),
      , scales = list(x = list(cex=1.2), y = list(cex=1.2)),
      , prepanel = prepanel.qqmathline,
      , panel = function(x, y,...) {
          panel.grid(lty=3, lwd=0.5, col="black",...)
          panel.qqmathline(x, y=x)
          panel.qqmath(x, y,...)
        }
    )
    print(a)
  dev.off()

}

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
ACF <- ddply(
    .data = tmp.remainder,
	.variables = "station.id",
	.fun = summarise,
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

dr <- ddply(
    .data = tmp.trend,
	.variables = c("station.id","year"),
	.fun = summarise,
	mean = mean(response)
)
mm <- dr[rep(row.names(dr), each=12),]
tmp.trend <- cbind(tmp.trend, mean= mm$mean)

order <- ddply(
    .data = tmp.trend,
    .variables = "station.id",
    .fun = summarise,
    mean = mean(response)
)
order.st <- as.character(order[order(order$mean, decreasing=TRUE), ]$station.id)
tmp.trend$station.id <- factor(tmp.trend$station.id, levels=order.st)
tmp.trend$time <- rep(0:1235,100)

trellis.device(
    device = postscript, 
    file = paste(
        outputdir, 
        "scatterplot_of_", 
        dataset, 
        "_with_stl2_fc.trend_for_100_stations",
        ".ps", sep = ""
    ), 
    color = TRUE, 
    paper = "legal"
)
b <- xyplot( mean ~ time | station.id, 
    data = tmp.trend,
    xlab = list(label = "Month", cex = 1.2),
    ylab = list(label = ylab, cex = 1.2),
#   strip = strip.custom(par.strip.text= list(cex = 1.5)),
#   par.settings = list(layout.heights = list(strip = 1.5)),
	xlim = c(0, 1235),
    pch = 16,
	layout = c(10,1),
	aspect="xy",
    key=list(
        text = list(label=c("low frequency component","yearly mean")), 
        lines = list(pch=c("","."), cex=4, lwd=1.5, type=c("l","p"), col=col[1:2]),
        columns=2
    ),
    cex = 0.3,
    scales = list(
        y = list(relation = 'free'), 
        x=list(at=seq(0, 1236, by=600), relation = 'same')
    ),
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

trellis.device(
    device = postscript, 
    file = paste(
        outputdir, "scatterplot_of_", dataset, 
        "_with_stl2_fc.second_for_100_stations",
        ".ps", sep = ""
    ), 
    color = TRUE, 
    paper = "legal"
)
b <- xyplot( (mean-fc.trend) ~ time | station.id,
    data = tmp.trend,
    xlab = list(label = "Month", cex = 1.2),
    ylab = list(label = ylab, cex = 1.2),
#   strip = strip.custom(par.strip.text= list(cex = 1.5)),
#   par.settings = list(layout.heights = list(strip = 1.5)),
    xlim = c(0, 1235),
    layout = c(2,1),
    aspect="xy",
    key=list(
        text = list(label=c(
            "middle frequency component",
            "yearly mean-low frequency component"
        )),
        lines = list(
            pch=c("","."), 
            cex=4, 
            lwd=1.5, 
            type=c("l","p"), 
            col=col[1:2]
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

######################################################
##Trend + remainder component and loess from stl2
######################################################
tmp.fc <- tmp[, c(!(names(tmp) %in% c("weights","trend", "remainder")))]
tmp.fc$tr <- tmp.fc$response - tmp.fc$seasonal
tmp.fc$time <- rep(0:1235,100)

trellis.device(
    device = postscript, 
    file = paste(
        outputdir, "scatterplot_of_", dataset, 
        "_after_seasonal_for_100_stations",".ps", sep = ""
    ), 
    color = TRUE, 
    paper = "legal"
)
b <- xyplot( tr ~ time | station.id,
    data = tmp.fc,
    xlab = list(label = "Month", cex = 1.2),
    ylab = list(label = ylab, cex = 1.2),
	layout = c(1,1),
    xlim = c(0, 1235),	
    key=list(
        type = "l", 
        text = list(
            label = c("after seasonal","loess smoothing")
        ),  
        lines = list(
            lty = c(1,6), 
            lwd = 1.5, 
            col = col[1:2]
        ), 
        columns = 2
    ),
    scales = list(
        y = list(relation = 'free'), 
        x = list(at=seq(0, 1236, by=300), relation = 'same')
    ),
    prepanel = function(x,y,...) {
        prepanel.loess(x,y,span=900/1236, degree=2, ...)
    },
    panel = function(x,y,...) {
        panel.loess(
            x = x,
            y = y,
            degree = 2, 
            span = 900/1236, 
            col = col[2], 
            type = "l", lty=1, ...
        )
        panel.xyplot(
            x = x,
            y = y,
            col = col[1], 
            type= "p", cex=0.3, ...
        )
    }
)
print(b)
dev.off()

rm(tmp.fc)
#####################################################
##time series plot of trend component loess from stl2
#####################################################
trellis.device(
    device = postscript, 
    file = paste(
        outputdir, "scatterplot_of_", dataset, 
        "_with_stl2_fc.trend_loess_for_100_stations",
        ".ps", sep = ""
    ), 
    color = TRUE, 
    paper = "legal"
)
b <- xyplot( fc.trend ~ time | station.id,
    data = tmp.trend[order(tmp.trend$time),],
    xlab = list(label = "Month", cex = 1.2),
    ylab = list(label = ylab, cex = 1.2),
#   strip = strip.custom(par.strip.text= list(cex = 1.5)),
#   par.settings = list(layout.heights = list(strip = 1.5)),
    aspect = "xy",
	layout = c(5,2),
	xlim = c(0, 1235),
    key = list(
        type = "l", 
        text = list(
            label = c("loess smoothing","trend component")
        ),  
        lines = list(
            lty = c(1,6), 
            lwd = 1.5, 
            col = col[1:2]
        ), 
        columns=2
    ),
    scales = list(
        y = list(relation = 'free'), 
        x = list(at=seq(0, 1236, by=300), relation = 'same')
    ),
	prepanel = function(x,y,...) {
        prepanel.loess(x, y, span = 3/4, degree = 1, ...)
    },
    panel = function(x,y,...) {
        panel.loess(
            x = x,
            y = y,
            degree = 1, 
            span = 3/4, 
            col = col[1], 
            type = "l", lty=1, ...
        )
		panel.xyplot(
            x = x, 
            y = y,
            col = col[2], 
            type= "l", lty = 6, ...
        )
    }
)
print(b)
dev.off()


##################################
##scatter plot pf range vs mean
##################################
dr <- ddply(
    .data = tmp,
    .variables = c("station.id", "year"), 
    .fun = summarise,
    mean = mean(response),
	range = max(response) - min(response)
)
order <- ddply(
    .data = tmp,
    .variables = "station.id",
    .fun = summarise,
    mean = mean(response)
)
order.st <- as.character(order[order(order$mean, decreasing=TRUE), ]$station.id)
dr$station.id <- factor(dr$station.id, levels=order.st)

trellis.device(
    device = postscript, 
    file = paste(
        outputdir, 
        "scatterplot_of_", 
        dataset, 
        "_range_vs_mean_for_100_stations",
        ".ps", 
        sep = ""
    ), 
    color=TRUE, 
    paper="legal"
)
	a <- xyplot(mean ~ range | station.id, 
		data = dr,
		pch = 16,
        xlab = list(label = "Yearly Range", cex = 1.2),
        ylab = list(label = "Yearly Mean", cex = 1.2),
		cex = 0.5,
		layout =c(10,1),
		scales = list(
            y = list(relation = 'same'), 
            x = list(relation = 'same')
        )
	)
	print(a)
dev.off()


#temporal residual against spatial covariates conditional on month
trellis.device(
    device = postscript, 
    file = paste(
        "temporal remainder vs. lat",
        ".ps", 
        sep = ""
    ), 
    color=TRUE, 
    paper="legal"
)
    a <- xyplot(fc.remainder ~ lat | date,
        data = tmp,
        pch = 16,
        cex = 0.6,
        layout = c(4,3)
    )
    print(a)
dev.off()


##########################################################################################
##QQ plot of the slopes of loess of trend component, ordered stations by average of slopes
##########################################################################################

#Get the loess estimate of trend component
#tmp.trend <- tmp[,c(1:9,15)]
#tmp.trend$time <- rep(0:1235,100)
#loess.trend <- ddply(.data = tmp.trend,
#           .variables = "station.id",
#           .fun = summarise,
#            elev = elev,
#            lon = lon,
#            lat = lat,
#            time = time,
#            loess = loess(trend ~ time, span=0.4, degree=2)$fitted
#)
#
##Since the x-axis are time with 1 unit, 
##the difference between two loess estimate is the approximate of the slope.
##Each station will have 1235 slope approximations.
#slope.trend <- ddply(.data = loess.trend,
#       .variables = "station.id",
#       .fun = summarise,
#        elev = elev[1:1235],
#        lon = lon[1:1235],
#        lat = lat[1:1235],
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
