##############################################################################
##Precipitation and temperature for the 100 most observations stations
##Raw data AGAINST time (month index).
##Difference between this code and old code file is change the period to be 9.
##############################################################################

#Set up the directory and load the data.
library(lattice)
datadir <- "~/Projects/Spatial/NCAR/RData/"
outputdir <- "~/Projects/Spatial/NCAR/output/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))

#Count the observation number for each station.
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)
TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)

###################################################
##scatter plot and line plot for precipitation
####################################################
Pcount <- sort(Pcount, decreasing=TRUE)
stations <- head(names(Pcount),100)

tmp <- USppt[USppt[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$precip)),]


month <- tmp$month
levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$date),]
start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp[,1]),]
tmp$factor <- factor(rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times=100), levels=paste("Period", c(9:1)))
tmp$time <- c(rep(0:143,8), 0:83) 

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_precipitation_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){   
	b <- xyplot( precip ~ time | factor,
             data = tmp[tmp$station.id==i,],
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.2),
	     main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = "b",
             pch=16,
             cex=0.5,
	     layout = c(1,9),
#	     strip.left = TRUE,
	     strip = FALSE,
#             aspect= 0.1,
             grib = TRUE,
             xlim = c(0, 143),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
	     scales = list(y = list(relation = 'same', alternating=TRUE), x=list(at=seq(0,143,by=12), relation='same')),
             panel = function(...) {
                  panel.abline(h=seq(0,max(tmp$precip),by=10), v=seq(0,145,by=12), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
        )
        print(b)
     }
dev.off()

#trellis.device(postscript, file = paste("lineplot_of_precipitation_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
#    for(i in stations){
#        b <- xyplot( precip ~ date | factor,
#             data = tmp[tmp$station.id==i,],
#             xlab = list(label = "Time", cex = 1.5),
#             ylab = list(label = "Precipitation (millimeters)", cex = 1.5),
#	     main = list(label = paste("Precipitation for station ", i, sep=""), cex=1.5),
#             type = "l",
#             layout = c(1,9),
#	     strip.left= TRUE,
#	     strip = FALSE,
#             aspect= 0.06,
#             grib = TRUE,
#             scales = list(y = list(relation = 'same', cex=1, alternating=TRUE), x=list(draw=FALSE, relation='free')),
##            xlim = c(start, end),
##             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
#             panel = function(...) {
#                  panel.abline(h=seq(0,max(tmp$precip),by=10), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
#                  panel.xyplot(...)
#             }
#        )
#        print(b)
#    }
#dev.off()

#########################################################
##scatter plot and line plot for maximum temperature
#########################################################
TMaxcount <- sort(TMaxcount, decreasing=TRUE)
stations <- head(names(TMaxcount),100)

tmp <- UStemp[UStemp[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$tmax)),]


month <- tmp$month
levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$date),]
start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp[,1]),]
tmp$factor <- factor(rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times=100), levels=paste("Period", c(9:1)))
tmp$time <- c(rep(0:143,8), 0:83)


#trellis.device(postscript, file = paste("lineplot_of_tmax_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
#    for(i in stations){
#        b <- xyplot( tmax ~ date | factor,
#             data = tmp[tmp$station.id==i,],
#             xlab = list(label = "Time", cex = 1.5),
#             ylab = list(label = "Temperature (degrees centigrade)", cex = 1.5),
#	     main = list(label = paste("Max Temperature for station ", i, sep=""), cex=1.5),
#             type = "l",
#             layout = c(1,9),
#	     strip.left=TRUE,
#	     strip = FALSE,
#             aspect= 0.06,
#             grib = TRUE,
##             xlim = c(start, end),
##             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
#             scales = list(y = list(relation = 'same', cex=1, alternating=TRUE), x=list(draw=FALSE, relation='free')),
#             panel = function(...) {
#                  panel.abline(h=seq(0,max(tmp$tmax),by=5), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
#                  panel.xyplot(...)
#             }
#        )
#        print(b)
#    }
#dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmax_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        b <- xyplot( tmax ~ time | factor,
             data = tmp[tmp$station.id==i,],
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.2),
             main = list(label = paste("Station ", i, sep=""), cex=1.5),
	     type = "b",
             pch=16,
             cex=0.5,	     
             layout = c(1,9),
	     strip = FALSE,
#	     strip.left = TRUE,
             aspect= 0.06,
             grib = TRUE,
             xlim = c(0, 143),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(at = seq(0, 143, by=12), relation='same')),
             panel = function(...) {
                  panel.abline(h=seq(0,max(tmp$tmax),by=5), v=seq(0,145,by=12), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
        )
        print(b)
    }
dev.off()

########################################################
##scatter plot and line plot for minimum temperature
#######################################################
TMincount <- sort(TMincount, decreasing=TRUE)
stations <- head(names(TMincount),100)

tmp <- UStemp[UStemp[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$tmin)),]


month <- tmp$month
levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$date),]
start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp[,1]),]
tmp$factor <- factor(rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times=100), levels=paste("Period", c(9:1)))
tmp$time <- c(rep(0:143,8), 0:83)



trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmin_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
    for(i in stations){
        b <- xyplot( tmin ~ time | factor,
             data = tmp[tmp$station.id==i,],
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.2),
	     main = list(label = paste("Station ", i, sep=""), cex=1.5),
	     strip = FALSE,
#	     strip.left = TRUE,
             type = "b",
             pch=16,
             cex=0.5,
             layout = c(1,9),
             aspect= 0.06,
             grib = TRUE,
             xlim = c(0, 143),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(at = seq(0,143, by=12), relation='same')),
             panel = function(...) {
                  panel.abline(h=seq(min(tmp$tmin),max(tmp$tmin),by=5), v=seq(0,145,by=12), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
        )
        print(b)
    }
dev.off()

#trellis.device(postscript, file = paste("lineplot_of_tmin_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
#    for(i in stations){
#	b <- xyplot( tmin ~ date | factor,
#             data = tmp[tmp$station.id==i,],
#             xlab = list(label = "Time", cex = 1.5),
#             ylab = list(label = "Temperature (degrees centigrade)", cex = 1.5),
#             main = list(label = paste("Min Temperature for station ", i, sep=""), cex=1.5),
#             type = "l",
#             layout = c(1,9),
#	     strip = FALSE,
#	     strip.left = TRUE,
#             aspect= 0.06,
#             grib = TRUE,
#       #      xlim = c(start, end),
##             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
#             scales = list(y = list(relation = 'same', cex=1, alternating=TRUE), x=list(draw=FALSE, relation='free')),
#             panel = function(...) {
#                  panel.abline(h=seq(min(tmp$tmin),max(tmp$tmin),by=5), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
#                  panel.xyplot(...)
#             }
#        )
#        print(b)
#    }
#dev.off()


