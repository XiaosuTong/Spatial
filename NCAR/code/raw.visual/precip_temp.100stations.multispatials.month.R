###################################################################################################
##Precipitation and temperature for the 100 most observations stations
##Row data against one of spatial factor conditional on month*another spatial factor
###################################################################################################

#Set up the directory and load the data.
library(lattice)
library(plyr)
datadir <- "~/Projects/Spatial/NCAR/RData/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))
USppt$month <- factor(USppt$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))
UStemp$month <- factor(UStemp$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))

#Count the observation number for each station.
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)
TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

##############################################################################################################
##scatter plot for precipitation vs latitude/longitude/elevation conditional on month*another spatial factor
###############################################################################################################
Pcount <- sort(Pcount, decreasing=TRUE)
stations <- head(names(Pcount),100)

tmp <- USppt[USppt[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$precip)),]

month <- tmp$month
levels(month) <- c(1:12)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$date),]
#start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
#end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$station.id, tmp$month, tmp$year),]


##precipitation vs latitude condtional on month*elevation
tmp <- tmp[order(tmp$elev, tmp$station.id),]
tmp$elev.fc <- rep(paste("Elev Level", 1:4, sep=" "), each=1236*25)

trellis.device(postscript, file = paste("scatterplot_of_precipitation_vs_lat_for_100_stations_conditional_month*elev",".ps", sep = ""), color=TRUE, paper="legal")
	b <- xyplot( precip ~ lat | month*as.factor(elev.fc),
             data = tmp[order(tmp$lat),],
#	     groups = group,
             xlab = list(label = "Latitude (degrees)", cex = 1.5),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.5),
             type = c("p"),
#	     col = col[2:1],
#             distribute.type= TRUE,
             pch=16,
	     lwd=2, 
             cex=0.5,
#             strip = strip.custom(par.strip.text= list(cex = 1.5)),
#             par.settings = list(layout.heights = list(strip = 1.5)),
	     layout = c(1,1),
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
	     scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 80,by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)		  
             }
        )
        print(b)
dev.off()

##precipitation vs latitude conditional on month*longitude
tmp <- tmp[order(tmp$lon, tmp$station.id),]
tmp$lon.fc <- rep(paste("Lon Level", 1:4, sep=" "), each=1236*25)

trellis.device(postscript, file = paste("scatterplot_of_precipitation_vs_lat_for_100_stations_conditional_month*lon",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( precip ~ lat | month*as.factor(lon.fc),
             data = tmp[order(tmp$lat),],
#            groups = group,
             xlab = list(label = "Latitude (degrees)", cex = 1.5),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.5),
             type = c("p"),
#            col = col[2:1],
#             distribute.type= TRUE,
             pch=16,
             lwd=2,
             cex=0.5,
	     strip = TRUE,
#             strip = strip.custom(par.strip.text= list(cex = 1.5)),
#             par.settings = list(layout.heights = list(strip = 1.5)),
             layout = c(1,1),
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 100,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()


########################################################################

##precipitation vs longitude condtional on month*elevation
trellis.device(postscript, file = paste("scatterplot_of_precipitation_vs_lon_for_100_stations_conditional_month*elev",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( precip ~ lon | month*as.factor(elev.fc),
             data = tmp[order(tmp$lon),],
#             groups = group,
             xlab = list(label = "Longitude (degrees)", cex = 1.5),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type= TRUE,
             pch=16,
             lwd=2,
             cex=0.5,
             layout = c(1,1),
             strip = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 80,by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()


##precipitation vs longitude condtional on month*latitude
tmp <- tmp[order(tmp$lat, tmp$station.id),]
tmp$lat.fc <- rep(paste("Lat Level", 1:4, sep=" "), each=1236*25)

trellis.device(postscript, file = paste("scatterplot_of_precipitation_vs_lon_for_100_stations_conditional_month*lat",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( precip ~ lon | month*as.factor(lat.fc),
             data = tmp[order(tmp$lon),],
#             groups = group,
             xlab = list(label = "Longitude (degrees)", cex = 1.5),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type= TRUE,
             pch=16,
             lwd=2,
             cex=0.5,
             layout = c(1,1),
             strip = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 80,by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

########################################################################

##precipitation vs elevation conditional on month*longitude
trellis.device(postscript, file = paste("scatterplot_of_precipitation_vs_elev_for_100_stations_conditional_month*lon",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( precip ~ log2(elev) | month*as.factor(lon.fc),
             data = tmp[order(tmp$elev),],
#             groups = group,
             xlab = list(label = "Log of Elevation (log base 2 meters)", cex = 1.5),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type= TRUE,
             pch=16,
             lwd=2,
             cex=0.5,
             layout = c(1,1),
             strip = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 80,by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()


##precipitation vs elevation conditional on month*latitude
trellis.device(postscript, file = paste("scatterplot_of_precipitation_vs_elev_for_100_stations_conditional_month*lat",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( precip ~ log2(elev) | month*as.factor(lat.fc),
             data = tmp[order(tmp$elev),],
#             groups = group,
             xlab = list(label = "Log of Elevation (log base 2 meters)", cex = 1.5),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type= TRUE,
             pch=16,
             lwd=2,
             cex=0.5,
             layout = c(1,1),
             strip = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 80,by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

####################################################################################################################
##scatter plot for maximum temperature vs latitude/longitude/elevation conditional on month*another spatial factor
####################################################################################################################
TMaxcount <- sort(TMaxcount, decreasing=TRUE)
stations <- head(names(TMaxcount),100)

tmp <- UStemp[UStemp[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$tmax)),]


month <- tmp$month
levels(month) <- c(1:12)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$date),]
#start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
#end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$station.id, tmp$month, tmp$year),]

##tmax vs latitude condtional on month*elevation
tmp <- tmp[order(tmp$elev, tmp$station.id),]
tmp$elev.fc <- rep(paste("Elev Level", 1:4, sep=" "), each=1236*25)

trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_lat_for_100_stations_conditional_month*elev",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ lat | month*as.factor(elev.fc),
             data = tmp[order(tmp$lat),],
#	     groups = group,
             xlab = list(label = "Latitude (degrees)", cex = 1.5),
             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             pch=16,
             cex=0.5,	 
	     lwd = 2,    
             layout = c(1,1),
	     strip = TRUE,
#             aspect= 1,
             grib = TRUE,
#             xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,50,by=2), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()


##tmax vs latitude condtional on month*longitude
tmp <- tmp[order(tmp$lon, tmp$station.id),]
tmp$lon.fc <- rep(paste("Lon Level", 1:4, sep=" "), each=1236*25)

trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_lat_for_100_stations_conditional_month*lon",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ lat | as.factor(lon.fc)*as.factor(elev.fc),
             data = tmp[order(tmp$lat),],
#            groups = group,
             xlab = list(label = "Latitude (degrees)", cex = 1.5),
             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             pch=16,
             cex=0.5,
             lwd = 2,
             layout = c(1,1),
             strip = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,50,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

##tmax vs latitude conditional on longitude*elevation for each month
trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_lat_for_100_stations_conditional_month*lon*elev",".ps", sep = ""), color=TRUE, paper="legal")
      for(i in levels(tmp$month)){
        b <- xyplot( tmax ~ lat | as.factor(lon.fc)*as.factor(elev.fc),
             data = subset(tmp[order(tmp$lat),], month==i),
#            groups = group,
             xlab = list(label = "Latitude (degrees)", cex = 1.5),
             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.5),
             main = list(label = i, cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             pch=16,
             cex=0.5,
             lwd = 2,
             layout = c(4,4),
             strip = FALSE,
             strip.left = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,50,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                 # panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
      }
dev.off()


#############################################################################

##tmax vs longitude conditional on month*latitude
tmp <- tmp[order(tmp$lat, tmp$station.id),]
tmp$lat.fc <- rep(paste("Lat Level", 1:4, sep=" "), each=1236*25)

trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_lon_for_100_stations_conditional_month*lat",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ lon | month*as.factor(lat.fc),
             data = tmp[order(tmp$lon),],
#             groups = group,
             xlab = list(label = "Longitude (degrees)", cex = 1.5),
             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             pch=16,
             cex=0.5,
             lwd = 2,
             layout = c(1,1),
             strip = TRUE,
#             aspect= 1,
             grib = TRUE,
#             xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,50,by=2), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

##tmax vs longitude conditional on month*elevation
trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_lon_for_100_stations_conditional_month*elev",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ lon | month*as.factor(elev.fc),
             data = tmp[order(tmp$lon),],
#             groups = group,
             xlab = list(label = "Longitude (degrees)", cex = 1.5),
             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             pch=16,
             cex=0.5,
             lwd = 2,
             layout = c(1,1),
             strip = TRUE,
#             aspect= 1,
             grib = TRUE,
#             xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,50,by=2), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

##tmax vs longitude conditional on latitude*elevation for each month
trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_lon_for_100_stations_conditional_month*lat*elev",".ps", sep = ""), color=TRUE, paper="legal")
      for(i in levels(tmp$month)){
        b <- xyplot( tmax ~ lon | as.factor(lat.fc)*as.factor(elev.fc),
             data = subset(tmp[order(tmp$lon),], month==i),
#            groups = group,
             xlab = list(label = "Longitude (degrees)", cex = 1.5),
             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.5),
             main = list(label = i, cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             pch=16,
             cex=0.5,
             lwd = 2,
             layout = c(4,4),
             strip = FALSE,
             strip.left = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,50,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                 # panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
      }
dev.off()


##############################################################

##tmax vs elevation condtional on month*latitude
trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_elev_for_100_stations_conditional_month*lat",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ log2(elev) | month*as.factor(lat.fc),
             data = tmp[order(tmp$elev),],
#             groups = group,
             xlab = list(label = "Log of Elevation (log base 2 meters)", cex = 1.5),
             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type= TRUE,
             pch=16,
             lwd=2,
             cex=0.5,
             layout = c(1,1),
             strip = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 80,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()


##tmax vs elevation condtional on month*longitude
trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_elev_for_100_stations_conditional_month*lon",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ log2(elev) | month*as.factor(lon.fc),
             data = tmp[order(tmp$elev),],
#             groups = group,
             xlab = list(label = "Log of Elevation (log base 2 meters)", cex = 1.5),
             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type= TRUE,
             pch=16,
             lwd=2,
             cex=0.5,
             layout = c(1,1),
             strip = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 80,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

##tmax vs elevation conditional on latitude*longitude for each month
trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_elev_for_100_stations_conditional_month*lat*lon",".ps", sep = ""), color=TRUE, paper="legal")
      for(i in levels(tmp$month)){
        b <- xyplot( tmax ~ log2(elev) | as.factor(lat.fc)*as.factor(lon.fc),
             data = subset(tmp[order(tmp$elev),], month==i),
#            groups = group,
             xlab = list(label = "Log of Elevation (log base 2 meters)", cex = 1.5),
             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.5),
             main = list(label = i, cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             pch=16,
             cex=0.5,
             lwd = 2,
             layout = c(4,4),
             strip = FALSE,
             strip.left = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,50,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                 # panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
      }
dev.off()

###########################################################################################
##scatter plot for minimum temperature vs latitude/longitude/elevation conditional on month
###########################################################################################
TMincount <- sort(TMincount, decreasing=TRUE)
stations <- head(names(TMincount),100)

tmp <- UStemp[UStemp[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$tmin)),]

month <- tmp$month
levels(month) <- c(1:12)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$date),]
#start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
#end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$station.id, tmp$month, tmp$year),]


##tmin vs latitude condtional on month*elevation
tmp <- tmp[order(tmp$elev, tmp$station.id),]
tmp$elev.fc <- rep(paste("Elev Level", 1:4, sep=" "), each=1236*25)

trellis.device(postscript, file = paste("scatterplot_of_tmin_vs_lat_for_100_stations_conditional_month*elev",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmin ~ lat | month*as.factor(elev.fc),
             data = tmp[order(tmp$lat),],
#	     groups = group,
             xlab = list(label = "Latitude (degrees)", cex = 1.5),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
	     strip = TRUE,
             pch=16,
             cex=0.5,
	     lwd =2,
             layout = c(1,1),
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,40,by=2), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loessax, y, degree=1, span=2/3, col = col[2], ...)
f latitude of the 100 stations for minimum temperature is from 24.55 to 46.90. Minimum temperature of 100 stations over 103 years for a given month is plotted against the latitude of stations on one page. For one station, 103 observations of minimum temperature are all plotted at which the latitude of the station is. A purple smooth line is drawn on the plot which is a loess curve of raw observation with degree=1, span=2/3. 
             }
        )
        print(b)
dev.off()


##tmin vs latitude condtional on month*longitude
tmp <- tmp[order(tmp$lon, tmp$station.id),]
tmp$lon.fc <- rep(paste("Lon Level", 1:4, sep=" "), each=1236*25)

trellis.device(postscript, file = paste("scatterplot_of_tmin_vs_lat_for_100_stations_conditional_month*lon",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmin ~ lat | month*as.factor(lon.fc),
             data = tmp[order(tmp$lat),],
#            groups = group,
             xlab = list(label = "Latitude (degrees)", cex = 1.5),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             strip = TRUE,
             pch=16,
             cex=0.5,
             lwd =2,
             layout = c(1,1),
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,40,by=2), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()



##############################################################################################

##tmin vs longitude condtional on month*latitude
tmp <- tmp[order(tmp$lat, tmp$station.id),]
tmp$lat.fc <- rep(paste("Lat Level", 1:4, sep=" "), each=1236*25)

trellis.device(postscript, file = paste("scatterplot_of_tmin_vs_lon_for_100_stations_conditional_month*lat",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmin ~ lon | month*as.factor(lat.fc),
             data = tmp[order(tmp$lon),],
#             groups = group,
             xlab = list(label = "Longitude (degrees)", cex = 1.5),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             strip = TRUE,
             pch=16,
             cex=0.5,
             lwd =2,
             layout = c(1,1),
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,40,by=2), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

#tmin vs longitude conditional on month*elevation
trellis.device(postscript, file = paste("scatterplot_of_tmin_vs_lon_for_100_stations_conditional_month*elev",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmin ~ lon | month*as.factor(elev.fc),
             data = tmp[order(tmp$lon),],
#             groups = group,
             xlab = list(label = "Longitude (degrees)", cex = 1.5),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             strip = TRUE,
             pch=16,
             cex=0.5,
             lwd =2,
             layout = c(1,1),
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,40,by=2), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

##############################################################

##tmin vs elevation conditional on month*latitude
trellis.device(postscript, file = paste("scatterplot_of_tmin_vs_elev_for_100_stations_conditional_month*lat",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmin ~ log2(elev) | month*as.factor(lat.fc),
             data = tmp[order(tmp$elev),],
#             groups = group,
             xlab = list(label = "Log of Elevation (log base 2 meters)", cex = 1.5),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type= TRUE,
             pch=16,
             lwd=2,
             cex=0.5,
             layout = c(1,1),
             strip = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-50, 80,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()


##tmin vs elevation conditional on month*longitude
trellis.device(postscript, file = paste("scatterplot_of_tmin_vs_elev_for_100_stations_conditional_month*lon",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmin ~ log2(elev) | month*as.factor(lon.fc),
             data = tmp[order(tmp$elev),],
#             groups = group,
             xlab = list(label = "Log of Elevation (log base 2 meters)", cex = 1.5),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type= TRUE,
             pch=16,
             lwd=2,
             cex=0.5,
             layout = c(1,1),
             strip = TRUE,
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='free')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-50, 80,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

