##############################################################################
##Pot the precipitation and temperature for the 100 most observations stations
##############################################################################

#Set up the directory and load the data.
library(lattice)
datadir <- "~/Projects/Spatial/NCAR/RData/"
#load(paste(datadir,"USmonthlyMet.RData", sep=""))

#Count the observation number for each station.
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)
TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)

#############################################
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
tmp$factor <- factor(rep(paste("Period", 1:12), each=103,times=100), levels=paste("Period", c(3,2,1,6,5,4,9,8,7,12,11,10)))
 

trellis.device(postscript, file = paste("scatterplot_of_precipitation_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( precip ~ date | factor*as.factor(station.id),
             data = tmp,
             xlab = list(label = "Time", cex = 1.5),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.5),
             type = "p",
	     layout = c(1,3),
             aspect= "xy",
             grib = TRUE,
#             xlim = c(start, end),
             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             panel = function(...) {
                  panel.abline(h=seq(0,max(tmp$precip),by=10), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
        )
        print(b)
dev.off()

trellis.device(postscript, file = paste("lineplot_of_precipitation_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( precip ~ date | factor*as.factor(station.id),
             data = tmp,
             xlab = list(label = "Time", cex = 1.5),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.5),
             type = "l",
             layout = c(1,3),
             aspect= "xy",
             grib = TRUE,
#            xlim = c(start, end),
             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             panel = function(...) {
                  panel.abline(h=seq(0,max(tmp$precip),by=10), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
        )
        print(b)
dev.off()

################################################################################
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
tmp$factor <- factor(rep(paste("Period", 1:12), each=103,times=100), levels=paste("Period", c(3,2,1,6,5,4,9,8,7,12,11,10)))



trellis.device(postscript, file = paste("lineplot_of_tmax_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ date | factor*as.factor(station.id),
             data = tmp,
             xlab = list(label = "Time", cex = 1.5),
             ylab = list(label = "Temperature (degrees centigrade)", cex = 1.5),
             type = "l",
             layout = c(1,3),
             aspect= "xy",
             grib = TRUE,
#             xlim = c(start, end),
             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             panel = function(...) {
                  panel.abline(h=seq(0,max(tmp$tmax),by=5), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
        )
        print(b)
dev.off()
trellis.device(postscript, file = paste("scatterplot_of_tmax_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ date | factor*as.factor(station.id),
             data = tmp,
             xlab = list(label = "Time", cex = 1.5),
             ylab = list(label = "Temperature (degrees centigrade)", cex = 1.5),
             type = "p",
             layout = c(1,3),
             aspect= "xy",
             grib = TRUE,
#             xlim = c(start, end),
             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             panel = function(...) {
                  panel.abline(h=seq(0,max(tmp$tmax),by=5), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
        )
        print(b)
dev.off()

###################################################################
TMincount <- sort(TMincount, decreasing=TRUE)
stations <- head(names(TMincount),100)

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
tmp$factor <- factor(rep(paste("Period", 1:12), each=103,times=100)[1:123579], levels=paste("Period", c(3,2,1,6,5,4,9,8,7,12,11,10)))



trellis.device(postscript, file = paste("scatterplot_of_tmin_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ date | factor*as.factor(station.id),
             data = tmp,
             xlab = list(label = "Time", cex = 1.5),
             ylab = list(label = "Temperature (degrees centigrade)", cex = 1.5),
             type = "p",
             layout = c(1,3),
             aspect= "xy",
             grib = TRUE,
      #       xlim = c(start, end),
             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             panel = function(...) {
                  panel.abline(h=seq(0,max(tmp$tmax),by=5), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
        )
        print(b)
dev.off()

trellis.device(postscript, file = paste("lineplot_of_tmin_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( tmax ~ date | factor*as.factor(station.id),
             data = tmp,
             xlab = list(label = "Time", cex = 1.5),
             ylab = list(label = "Temperature (degrees centigrade)", cex = 1.5),
             type = "l",
             layout = c(1,3),
             aspect= "xy",
             grib = TRUE,
       #      xlim = c(start, end),
             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             panel = function(...) {
                  panel.abline(h=seq(0,max(tmp$tmax),by=5), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
        )
        print(b)
dev.off()


