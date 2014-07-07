###################################################################################################
##Precipitation and temperature, 100 most observations stations 
##Raw data AGAINST one of spatial factor conditional on month
###################################################################################################

#Set up the directory and load the data.
library(lattice)
library(plyr)
datadir <- "~/Projects/Spatial/NCAR/RData/"
outputdir <- "~/Projects/Spatial/NCAR/output/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))
USppt$month <- factor(USppt$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))
UStemp$month <- factor(UStemp$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))

#Count the observation number for each station.
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)
TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

#####################################################################################
##scatter plot for precipitation vs latitude/longitude/elevation conditional on month
#####################################################################################
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
start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$station.id, tmp$month, tmp$year),]

########################################################################################
#A very naive way to calculate the loess estimate of each observation. If the number 
#of observation is relative large, then there is not too much difference between manully 
#calculating loess estimate and using panel.loess function.
#
#lat.dr <- ddply(.data=tmp,
#	    .variables= c("month"), 
#	    .fun = summarise,
#	     station.id = station.id,
#	     fitted = loess(precip~lat, degree=1, span=1/3)$fitted)
#
#lat.dr <- lat.dr[order(lat.dr$station.id, lat.dr$month),]
#tmp.lat <- cbind(tmp, loess=lat.dr$fitted)
#tmp1 <- tmp.lat[,c(1:7, 9,8)]
#tmp2 <- tmp.lat[,c(1:7, 9,10)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp.lat <- rbind(tmp2, tmp1)
#tmp.lat$group <- rep(c("loess","raw"), each=123600)
########################################################################################


##precipitation vs latitude conditional on month
trellis.device(postscript, file = paste("scatterplot_of_precipitation_vs_lat_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
	b <- xyplot( precip ~ lat | month,
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
	     layout = c(1,1),
             strip = strip.custom(par.strip.text= list(cex = 1.5)),
             par.settings = list(layout.heights = list(strip = 1.5)),

#             aspect= "xy",
             grib = TRUE,
#             xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
	     scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 80,by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
		  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()


#########################################################################
#lon.dr <- ddply(.data=tmp,
#		.variables = "month",
#		.fun = summarise,
#		 station.id = station.id,
#		 fitted = loess(precip~lon, degree=1, span=1/3)$fitted)
#
#lon.dr <- lon.dr[order(lon.dr$station.id, lon.dr$month),]
#tmp.lon <- cbind(tmp, loess=lon.dr$fitted)
#tmp1 <- tmp.lon[, c(1:7, 9,8)]
#tmp2 <- tmp.lon[, c(1:7, 9,10)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp.lon <- rbind(tmp2, tmp1)
#tmp.lon$group <- rep(c("loess","raw"), each=123600)
##########################################################################

##precipitation vs longitude conditional month
trellis.device(postscript, file = paste("scatterplot_of_precipitation_vs_lon_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( precip ~ lon | month,
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
             strip = strip.custom(par.strip.text= list(cex = 1.5)),
             par.settings = list(layout.heights = list(strip = 1.5)),
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 80,by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

###############################################################################
#elev.dr <- ddply(.data=tmp,
#                .variables = "month",
#                .fun = summarise,
#                 station.id = station.id,
#                 fitted = loess(precip~log2(elev), degree=1, span=1/3)$fitted)
#
#elev.dr <- elev.dr[order(elev.dr$station.id, elev.dr$month),]
#tmp.elev <- cbind(tmp, loess=elev.dr$fitted)
#tmp1 <- tmp.elev[, c(1:7, 9,8)]
#tmp2 <- tmp.elev[, c(1:7, 9,10)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp.elev <- rbind(tmp2, tmp1)
#tmp.elev$group <- rep(c("loess","raw"), each=123600)
#################################################################################

##Precipitation vs elevation conditional on month
trellis.device(postscript, file = paste("scatterplot_of_precipitation_vs_elev_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( precip ~ log2(elev) | month,
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
             strip = strip.custom(par.strip.text= list(cex = 1.5)),
             par.settings = list(layout.heights = list(strip = 1.5)),
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 80,by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()


###########################################################################################
##scatter plot for maximum temperature vs latitude/longitude/elevation conditional on month
###########################################################################################
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


########################################################################
#lat.dr <- ddply(.data=tmp,
#            .variables= "month",
#            .fun = summarise,	    
#             station.id = station.id,
#             fitted = loess(tmax ~ lat, degree = 1, span = 1/3)$fitted)
#
#lat.dr <- lat.dr[order(lat.dr$station.id, lat.dr$month),]
#tmp.lat <- cbind(tmp, loess=lat.dr$fitted)
#tmp1 <- tmp.lat[,c(1:7, 10,9)]
#tmp2 <- tmp.lat[,c(1:7,10,11)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp.lat <- rbind(tmp2, tmp1)
#tmp.lat$group <- rep(c("loess","raw"), each=123600)
#########################################################################

##maximum temperature vs latitude conditional on month
trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_lat_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ lat | month,
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
             strip = strip.custom(par.strip.text= list(cex = 1.5)),
             par.settings = list(layout.heights = list(strip = 1.5)),
#             aspect= 1,
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,50,by=2), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

#############################################################################
#lon.dr <- ddply(.data=tmp,
#            .variables= "month",
#            .fun = summarise,
#             station.id = station.id,
#             fitted = loess(tmax ~ lon, degree = 1, span = 1/3)$fitted)
#
#lon.dr <- lon.dr[order(lon.dr$station.id, lon.dr$month),]
#tmp.lon <- cbind(tmp, loess=lon.dr$fitted)
#tmp1 <- tmp.lon[,c(1:7, 10,9)]
#tmp2 <- tmp.lon[,c(1:7,10,11)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp.lon <- rbind(tmp2, tmp1)
#tmp.lon$group <- rep(c("loess","raw"), each=123600)
##############################################################################

##maximum temperature vs longitude conditional on month
trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_lon_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ lon | month,
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
             strip = strip.custom(par.strip.text= list(cex = 1.5)),
             par.settings = list(layout.heights = list(strip = 1.5)),
#             aspect= 1,
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,50,by=2), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()


#############################################################################
#elev.dr <- ddply(.data=tmp,
#                .variables = "month",
#                .fun = summarise,
#                 station.id = station.id,
#                 fitted = loess(tmax~log2(elev), degree=1, span=1/3)$fitted)
#
#elev.dr <- elev.dr[order(elev.dr$station.id, elev.dr$month),]
#tmp.elev <- cbind(tmp, loess=elev.dr$fitted)
#tmp1 <- tmp.elev[, c(1:7, 10,9)]
#tmp2 <- tmp.elev[, c(1:7, 10,11)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp.elev <- rbind(tmp2, tmp1)
#tmp.elev$group <- rep(c("loess","raw"), each=123600)
##############################################################################

##maximum temperature vs elevation conditional on month
trellis.device(postscript, file = paste("scatterplot_of_tmax_vs_elev_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ log2(elev) | month,
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
#             aspect= "xy",
             grib = TRUE,
             strip = strip.custom(par.strip.text= list(cex = 1.5)),
             par.settings = list(layout.heights = list(strip = 1.5)),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-20, 80,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
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

#####################################################################
#lat.dr <- ddply(.data=tmp,
#            .variables= "month",
#            .fun = summarise,
#	     station.id=station.id,
#             fitted = loess(tmin ~ lat, degree=1, span=1/3)$fitted)
#
#lat.dr <- lat.dr[order(lat.dr$station.id, lat.dr$month),]
#tmp.lat <- cbind(tmp, loess=lat.dr$fitted)
#tmp1 <- tmp.lat[,c(1:7,10,8)]
#tmp2 <- tmp.lat[,c(1:7,10,11)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp.lat <- rbind(tmp2, tmp1)
#tmp.lat$group <- rep(c("loess","raw"), each=123600)
######################################################################

##minimum temperature vs latitude conditional on month
trellis.device(postscript, file = paste("scatterplot_of_tmin_vs_lat_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmin ~ lat | month,
             data = tmp[order(tmp$lat),],
#	     groups = group,
             xlab = list(label = "Latitude (degrees)", cex = 1.5),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             pch=16,
             cex=0.5,
	     lwd =2,
             layout = c(1,1),
#             aspect= "xy",
             grib = TRUE,
             strip = strip.custom(par.strip.text= list(cex = 1.5)),
             par.settings = list(layout.heights = list(strip = 1.5)),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,40,by=2), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()


########################################################################
#lon.dr <- ddply(.data=tmp,
#            .variables= "month",
#            .fun = summarise,
#             station.id=station.id,
#             fitted = loess(tmin ~ lon, degree=1, span=1/3)$fitted)
#
#lon.dr <- lon.dr[order(lon.dr$station.id, lon.dr$month),]
#tmp.lon <- cbind(tmp, loess=lon.dr$fitted)
#tmp1 <- tmp.lon[,c(1:7,10,8)]
#tmp2 <- tmp.lon[,c(1:7,10,11)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp.lon <- rbind(tmp2, tmp1)
#tmp.lon$group <- rep(c("loess","raw"), each=123600)
#########################################################################

##minimum temperature vs longitude conditional on month
trellis.device(postscript, file = paste("scatterplot_of_tmin_vs_lon_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmin ~ lon | month,
             data = tmp[order(tmp$lon),],
#             groups = group,
             xlab = list(label = "Longitude (degrees)", cex = 1.5),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.5),
             type = c("p"),
#             col = col[2:1],
#             distribute.type = TRUE,
             pch=16,
             cex=0.5,
             lwd =2,
             layout = c(1,1),
#             aspect= "xy",
             grib = TRUE,
             strip = strip.custom(par.strip.text= list(cex = 1.5)),
             par.settings = list(layout.heights = list(strip = 1.5)),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-40,40,by=2), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()


##############################################################################
#elev.dr <- ddply(.data=tmp,
#                .variables = "month",
#                .fun = summarise,
#                 station.id = station.id,
#                 fitted = loess(tmin~log2(elev), degree=1, span=1/3)$fitted)
#
#elev.dr <- elev.dr[order(elev.dr$station.id, elev.dr$month),]
#tmp.elev <- cbind(tmp, loess=elev.dr$fitted)
#tmp1 <- tmp.elev[, c(1:7, 10,8)]
#tmp2 <- tmp.elev[, c(1:7, 10,11)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp.elev <- rbind(tmp2, tmp1)
#tmp.elev$group <- rep(c("loess","raw"), each=123600)
################################################################################

##minimum temperature vs elevation conditional on month
trellis.device(postscript, file = paste("scatterplot_of_tmin_vs_elev_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmin ~ log2(elev) | month,
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
             strip = strip.custom(par.strip.text= list(cex = 1.5)),
             par.settings = list(layout.heights = list(strip = 1.5)),
             layout = c(1,1),
#             aspect= "xy",
             grib = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', cex=1, alternating=TRUE), x=list(tick.number=10, cex=1.2, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-50, 80,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  panel.loess(x, y, degree=1, span=2/3, col = col[2], ...)
             }
        )
        print(b)
dev.off()

