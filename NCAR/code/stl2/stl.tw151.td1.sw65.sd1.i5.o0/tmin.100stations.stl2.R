############################################################################
##Min temperature, 100 stations, seasonal+trend with raw data by stl2
##########################################################################

#load the data
library(lattice)
library(plyr)
library(stl2)
library(maps)
datadir <- "~/Projects/Spatial/NCAR/RData/"
outputdir <- "~/Projects/Spatial/NCAR/output/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)


#find the 100 stations with largest observation number for min temperature
TMincount <- sort(TMincount, decreasing=TRUE)
stations <- head(names(TMincount),100)

UStemp$month <- factor(UStemp$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))
tmp <- UStemp[UStemp[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$tmin)),]

month <- tmp$month
levels(month) <- c(1:12)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$date),]
tmp <- tmp[order(tmp[,1]),-9]
tmp$factor <- factor(rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times=100), levels=paste("Period", c(9:1)))
tmp$time <- c(rep(0:143,8), 0:83)

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

dr <- ddply(tmp, "station.id", function(r){stl2(r$tmin, r$date, n.p=12, s.window=65, t.window=151, inner = 5, outer = 0)$data})
tmp <- cbind(tmp, dr)
tmp <- tmp[order(tmp$station.id,tmp$date),]

##################################
##trend+seasonal time series plot
##################################
tmp1 <- tmp[,1:11]
names(tmp1)[8] <- "response"
tmp1$group <- rep("raw", 123600)
tmp1 <- tmp1[,c(1:7,9,10,11,8,12)]
tmp2 <- tmp[,c(1:7,9,10,11,14,15)]
tmp2$response <- tmp2$seasonal + tmp2$trend
tmp2 <- tmp2[,-c(11,12)]
tmp2$group <- rep("seasonal+trend", 123600)

rst <- rbind(tmp1, tmp2)

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmin_with_stl2_trend+seasonal_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        b <- xyplot( response ~ time | factor,
             data = rst[rst$station.id==i,],
	     groups = group,
             xlab = list(label = "Time", cex = 1.2),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.2),
	     main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = c("p","l"),
	     distribute.type= TRUE,
             layout = c(1,9),
             aspect= 0.06,
	     pch = 16,
	     cex = 0.5,
	     strip = FALSE,
#	     strip.left = TRUE,
             grib = TRUE,
             xlim = c(0, 143),
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(at=seq(0,143,by=12), relation='same')),
             panel = function(...) {
                  panel.abline(h=seq(min(rst$response), max(rst$response),by=5), v=seq(0,145,by=12), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
                 
             }
        )
        print(b)
   }
dev.off()

##################################################################################
#QQ plot and time series of Pooled remainder of min temperature for 100 stations
##################################################################################

#Create the QQ plot of temperature for one station
trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_stl2_remainder_tmin_of_100_stations", ".ps", sep = ""), color=TRUE, paper="legal")
        a <- qqmath(~ remainder | station.id,
                data = tmp,
                distribution = qnorm,
		aspect = 1,
		pch=16,
		cex=0.5,
		layout= c(6,3),
#                main = list(label= paste("Station ", i, sep=""), cex=2),
                xlab = list(label="Unit normal quantile", cex=1.2),
                ylab = list(label="Minimum temperature(degrees centigrade)", cex=1.2),
#                scales = list(x = list(cex=1.5), y = list(cex=1.5)),
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
tmp.remainder$factor <- factor(tmp$factor, levels=paste("Period", c(3,2,1,6,5,4,9,8,7)))
trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmin_with_stl2_remainder_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        b <- xyplot( remainder ~ time | factor,
             data = subset(tmp.remainder, station.id==i),
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.2),
             main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = "p",
             layout = c(1,3),
             pch = 16,
             cex = 0.5,
#             aspect= "xy",
	     xlim = c(0,143),
             strip = FALSE,
#             strip.left = TRUE,
             grib = TRUE,
#            xlim = c(start, end),
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(at = seq(0, 143, by=12), relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=0, v=seq(0,145, by=12), color="lightgrey", lty=3, lwd=0.5)
		  panel.loess(x,y,degree=2,span=1/4, col=col[2])
                  panel.xyplot(x,y,...)
             }
        )
        print(b)
   }
dev.off()


##################################################
##time series plot of seaosnal+remainder from stl2
##################################################
tmp$sr <- tmp$seasonal + tmp$remainder

trellis.device(postscript, file = paste("scatterplot_of_tmin_with_stl2_seasonal+remainder_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        b <- xyplot( sr ~ time | factor,
             data = subset(tmp, station.id==i),
             xlab = list(label = "Month", cex = 1.5),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.5),
             main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = "p",
             layout = c(1,9),
             pch = 16,
             cex = 0.5,
             aspect= 0.06,
             strip = FALSE,
             strip.left = TRUE,
             grib = TRUE,
#            xlim = c(start, end),
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='same')),
             panel = function(...) {
                  panel.abline(h=seq(-10,10,by=5), v=seq(0,140, by=12), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)

             }
        )
        print(b)
   }
dev.off()


###########################################################################
##time series plot of trend component and yearly average of raw from stl2
###########################################################################
tmp.trend <- tmp[,c(1:9,15)]

dr <- ddply(.data = tmp.trend,
                .variables = c("station.id","year"),
                .fun = summarise,
                 mean = mean(tmin)
)
mm <- dr[rep(row.names(dr), each=12),]
tmp.trend <- cbind(tmp.trend, mean= mm$mean)

order <- ddply(.data = tmp.trend,
                .variables = "station.id",
                .fun = summarise,
                 mean = mean(tmin)
)
order.st <- as.character(order[order(order$mean, decreasing=TRUE), ]$station.id)
tmp.trend$station.id <- factor(tmp.trend$station.id, levels=order.st)

tmp1 <- tmp.trend[,1:10]
tmp1$time <- rep(0:1235,100)
tmp2 <- tmp.trend[,c(1:9,11)]
tmp2$time <- rep(0:1235,100)
names(tmp1)[10] <- names(tmp2)[10] <- "response"
tmp.trend <- rbind(tmp1, tmp2)
tmp.trend$group <- rep(c("trend","yearly mean"), each = 123600)

#Attend to add a small map in the corner of the trend loess plot
us.map <- map('state', plot = FALSE, fill = TRUE)

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmin_with_stl2_trend_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in levels(tmp.trend$station.id)){
        b <- xyplot( response ~ time,
             data = subset(tmp.trend, station.id==i),
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.2),
             type = c("l","p"),
             groups = group,
             distribute.type= TRUE,
             pch = 16,
	     xlim = c(0,1235),
             cex = 0.3,
             key=list(type="l", text=list(label=c("Trend component","Yearly average")),  lines=list(lwd=1.5, col=col[1:2]), columns=2),
             aspect = 1,
#            strip = strip.custom(par.strip.text= list(cex = 1.5)),
#            par.settings = list(layout.heights = list(strip = 1.5)),
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'same',format = "%b %Y", tick.number=8), cex=1.2),
             scales = list(y = list(relation = 'free'), x=list(at=seq(0, 1235, by=120), relation = 'same')),
             panel = function(...) {
                  panel.abline(h=seq(-10,30,by=1), v=seq(0, 1236, by=120), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
        )
        a <- xyplot(lat ~ lon,
                data = subset(tmp.trend, station.id==i),
                xlab = NULL,
                ylab = NULL,
                pch = 16,
                cex = 0.4,
                xlim = c(-127,-65),
                ylim = c(23,50),
                col = "red",
#               scales = list(x = list(draw=FALSE), y = list(draw=FALSE)),
#               strip = strip.custom(par.strip.text= list(cex = 1.5)),
#               par.settings = list(layout.heights = list(strip = 1.5)),
                panel = function(x,y,...) {
                        panel.polygon(us.map$x,us.map$y,lwd=0.2)
                        panel.xyplot(x,y,...)
                }
        )
        plot.new()
        title(paste("Station ", i, sep=""), cex=1.5)
        print(b, pos=c(0.4,0,1,0.95), newpage=FALSE, more=TRUE)
        print(a, pos=c(0,0.55,0.4,1), more=FALSE)
   }
dev.off()


#####################################################
##Time series plot of trend+remainder against time
#####################################################
tmp1 <- tmp[,c(1:9,15, 16)]
tmp1$tr <- tmp1$trend + tmp1$remainder
tmp2 <- tmp[, c(1:9, 15)]
names(tmp2)[10] <- names(tmp1)[12] <- "response"
tmp1$group <- rep("trend+remainder", 123600)
tmp2$group <- rep("trend", 123600)
tmp2$time <- rep(0:1235,100)
tmp1$time <- rep(0:1235,100)
tmp.tr <- do.call("rbind", list(tmp1[,-c(10:11)], tmp2))

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmin_with_stl2_trend+remainder_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        b <- xyplot( response ~ time,
             data = subset(tmp.tr, station.id==i),
             main = list(label = paste("Station ", i, sep=""), cex=1.5),
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.2),
             groups = group,
             distribute.type = TRUE,
             type = c("l","p"),
#            strip = strip.custom(par.strip.text= list(cex = 1.5)),
#            par.settings = list(layout.heights = list(strip = 1.5)),
             xlim = c(0, 1235),
             pch = 16,
             aspect="xy",
             key=list(type="l", text=list(label=c("Trend","Trend+remainder", "loess")),  lines=list(lwd=1.5, col=col[1:3]), columns=3),
             cex = 0.3,
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'same',format = "%b %Y", tick.number=8), cex=1.2),
             scales = list(y = list(relation = 'free'), x=list(at=seq(0, 1235, by=120), relation = 'same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(0,40,by=2), v=seq(0,1236, by=120), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x,y,span=1/3, degree=2, col = col[3], ...)

             }
        )
        print(b)
    }
dev.off()

################################################
##time series plot of trend component loess from stl2
################################################
tmp.trend <- tmp[,c(1:9,15)]
tmp.trend$time <- rep(0:1235,100)

order <- ddply(.data = tmp.trend,
                .variables = "station.id",
                .fun = summarise,
                 mean = mean(tmin)
)
order.st <- as.character(order[order(order$mean, decreasing=TRUE), ]$station.id)
tmp.trend$station.id <- factor(tmp.trend$station.id, levels=order.st)

#Attend to add a small map in the corner of the trend loess plot
us.map <- map('state', plot = FALSE, fill = TRUE)

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmin_with_stl2_trend_loess_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in levels(tmp.trend$station.id)){
        b <- xyplot( trend ~ time,
             data = subset(tmp.trend[order(tmp.trend$time),], station.id==i),
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.2),
#            layout = c(1,1),
#            strip = strip.custom(par.strip.text= list(cex = 1.5)),
#            par.settings = list(layout.heights = list(strip = 1.5)),
             pch = 16,
	     lwd = 1,
#            aspect="xy",
	     xlim = c(0, 1235),
#            grib = TRUE,
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'same',format = "%b %Y", tick.number=8), cex=1.2),
             scales = list(y = list(relation = 'free'), x=list(at=seq(0, 1235, by=120), relation = 'same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-10,40,by=0.5), v=seq(0,1236, by=120), color="lightgrey", lty=3, lwd=0.5)
                  panel.loess(x,y,degree=1, span=1/3,...)
#                 panel.xyplot(x,y,...)
             }
        )
        a <- xyplot(lat ~ lon,
                data = subset(tmp.trend, station.id==i),
                xlab = NULL,
                ylab = NULL,
#                layout = c(1,1),
                pch = 16,
                cex = 0.4,
                xlim = c(-127,-65),
                ylim = c(23,50),
                col = "red",
#               scales = list(x = list(draw=FALSE), y = list(draw=FALSE)),
#                strip = strip.custom(par.strip.text= list(cex = 1.5)),
#                par.settings = list(layout.heights = list(strip = 1.5)),
                panel = function(x,y,...) {
                        panel.polygon(us.map$x,us.map$y,lwd=0.2)
                        panel.xyplot(x,y,...)
                }
        )
        plot.new()
        title(paste("Station ", i, sep=""), cex=1.5)
        print(b, pos=c(0.4,0,1,0.95), newpage=FALSE, more=TRUE)
        print(a, pos=c(0,0.55,0.4,1), more=FALSE)
   }
dev.off()


##########################################################################################
##QQ plot of the slopes of loess of trend component, ordered stations by average of slopes
##########################################################################################

#Get the loess estimate of trend component
tmp.trend <- tmp[,c(1:9,15)]
tmp.trend$time <- rep(0:1235,100)
loess.trend <- ddply(.data = tmp.trend,
                        .variables = "station.id",
                        .fun = summarise,
                         elev = elev,
                         lon = lon,
                         lat = lat,
                         time = time,
                         loess = loess(trend ~ time, span=0.4, degree=2)$fitted
)

#Since the x-axis are time with 1 unit, 
#the difference between two loess estimate is the approximate of the slope.
#Each station will have 1235 slope approximations.
slope.trend <- ddply(.data = loess.trend,
                .variables = "station.id",
                .fun = summarise,
                 elev = elev[1:1235],
                 lon = lon[1:1235],
                 lat = lat[1:1235],
                 slope = diff(loess)
)

#Calculate the mean of slopes for each station, and then order the stations by mean slope.
mean.slope <- ddply(.data = slope.trend,
                .variables = "station.id",
                .fun = summarise,
                 elev = unique(elev),
                 lon = unique(lon),
                 lat = unique(lat),
                 mean = mean(slope)
)
mean.slope <- mean.slope[order(mean.slope$mean),]
station.or <- as.character(mean.slope$station.id)

#Attend to add a small map in the corner of the trend loess plot
us.map <- map('state', plot = FALSE, fill = TRUE)

trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_tmin_slope_of_trend_loess",".ps", sep = ""), color=TRUE, paper="legal")
    for(i in station.or){
        a <- qqmath(~ slope,
                data = subset(slope.trend, station.id == i),
                distribution = qunif,
                aspect = 1,
                pch = 16,
                cex = 0.5,
#                main = list(label= paste("Station ", i, sep=""), cex=2),
                xlab = list(label="f-value", cex=1.2),
                ylab = list(label="Minimum temperature(degrees centigrade)", cex=1.2),
#                scales = list(x = list(cex=1.5), y = list(cex=1.5)),
                prepanel = prepanel.qqmathline,
                panel = function(x, y,...) {
                        #panel.grid(lty=3, lwd=0.5, col="black",...)
                        panel.qqmath(x, y,...)
                }
        )
        b <- xyplot(lat ~ lon,
                data = subset(tmp.trend, station.id==i),
                xlab = NULL,
                ylab = NULL,
                pch = 16,
                cex = 0.4,
                xlim = c(-127,-65),
                ylim = c(23,50),
                col = "red",
#               scales = list(x = list(draw=FALSE), y = list(draw=FALSE)),
#               strip = strip.custom(par.strip.text= list(cex = 1.5)),
#               par.settings = list(layout.heights = list(strip = 1.5)),
                panel = function(x,y,...) {
                        panel.polygon(us.map$x,us.map$y,lwd=0.2)
                        panel.xyplot(x,y,...)
                }
        )
        plot.new()
        title(paste("Station ", i, sep=""), cex=1.5)
        print(a, pos=c(0.4,0,1,1), newpage=FALSE, more=TRUE)
        print(b, pos=c(0,0.5,0.4,1), more=FALSE)
    }
dev.off()

#scatter plot of average slops vs spatial factor
trellis.device(postscript, file = paste(outputdir, "scatter_plot_of_average_slopes_tmin",".ps", sep = ""), color=TRUE, paper="legal")
        c <- xyplot(mean ~ lat,
                data = mean.slope,
                xlab = list(label = "Elevation (meter)", cex=1.2),
                ylab = list(label = "Average slope of trend", cex=1.2),
                pch = 16,
                cex = 0.7,
#               scales = list(x = list(draw=FALSE), y = list(draw=FALSE)),
#               strip = strip.custom(par.strip.text= list(cex = 1.5)),
#               par.settings = list(layout.heights = list(strip = 1.5)),
                panel = function(x,y,...) {
                        panel.xyplot(x,y,...)
                }
        )
        print(c)
dev.off()


#####################################################
##trend for a given level of longitude and latitude
#####################################################
#tmp.trend <- tmp.trend[order(tmp.trend$lon, tmp.trend$station.id),]
#tmp.trend$lon.fac <- rep(1:5,each=1236*2*20)
#tmp.trend <- tmp.trend[order(tmp.trend$lat, tmp.trend$station.id),]
#tmp.trend$lat.fac <- rep(1:5, each=1236*2*20)

#trellis.device(postscript, file = paste("scatterplot_of_tmin_with_stl2_trend_for_100_stations_ordered",".ps", sep = ""), color=TRUE, paper="legal")
#   for(j in 1:5){
#      for(i in 1:5){
#	if(dim(subset(tmp.trend, lon.fac==i&lat.fac==j))[1]!=0){
#        b <- xyplot( response ~ date | station.id,
#             data = subset(tmp.trend, lon.fac==i&lat.fac==j),
#             xlab = list(label = "Month", cex = 1.5),
#             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.5),
#	     main = list(label = paste("Longitude Level ", i, " and Latitude Level ", j, sep=""), cex = 1.5),
#             type = c("l","p"),
#	     groups = group,
#             distribute.type= TRUE,
##	     strip.left = TRUE,
#             layout = c(1,1),
#             pch = 16,
#             cex = 0.3,
##	     aspect = "xy",
#             strip = strip.custom(par.strip.text= list(cex = 1.5)),
#             par.settings = list(layout.heights = list(strip = 1.5)),
##             grib = TRUE,
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'same',format = "%b %Y", tick.number=8), cex=1.2),
##             scales = list(y = list(relation = 'same', cex=1), x=list(tick.number=10, cex=1.2, relation = 'same')),
#            panel = function(...) {
#                  panel.abline(h=seq(-10,30,by=1), v=seq(start, end, by="60 month"), color="lightgrey", lty=3, lwd=0.5)
#                  panel.xyplot(...)
#             }
#        )
#        print(b)
#	}
#     }
#  }
#dev.off()



