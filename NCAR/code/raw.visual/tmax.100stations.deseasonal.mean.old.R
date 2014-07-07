##############################################################
##Max temperature, 100 stations, deseasonalizing
#############################################################

#load the data
library(lattice)
library(plyr)
datadir <- "~/Projects/Spatial/NCAR/RData/"
#load(paste(datadir,"USmonthlyMet.RData", sep=""))
TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)


#find the 100 stations with largest observation number for max temperature
TMaxcount <- sort(TMaxcount, decreasing=TRUE)
stations <- head(names(TMaxcount),100)

UStemp$month <- factor(UStemp$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))
tmp <- UStemp[UStemp[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$tmax)),]

month <- tmp$month
levels(month) <- c(1:12)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$date),]
tmp <- tmp[order(tmp[,1]),-8]
tmp$factor <- factor(rep(paste("Period", 1:12), each=103,times=100), levels=paste("Period", c(3,2,1,6,5,4,9,8,7,12,11,10)))

tmp <- tmp[order(tmp$station.id, tmp$month),]
a <- ddply(tmp, c("station.id","month"), function(r){c(average=mean(r$tmax))})
ave <- rep(a$average, each=103)

tmp$average <- ave
tmp$deseasonal <- tmp$tmax- tmp$average
tmp <- tmp[order(tmp$station.id,tmp$date),]

#start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
#end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
#
#trellis.device(postscript, file = paste("lineplot_of_deseasonal_tmax_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
#        b <- xyplot( deseasonal ~ date | factor*as.factor(station.id),
#             data = tmp,
#             xlab = list(label = "Time", cex = 1.5),
#             ylab = list(label = "Temperature (degrees centigrade)", cex = 1.5),
#             type = "l",
#             layout = c(1,3),
#             aspect= "xy",
#             grib = TRUE,
##            xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
#             panel = function(...) {
#                  panel.abline(h=seq(0,max(tmp$tmax),by=10), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
#                  panel.xyplot(...)
#             }
#        )
#        print(b)
#dev.off()
#
#trellis.device(postscript, file = paste("scatterplot_of_deseasonal_tmax_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
#        b <- xyplot( deseasonal ~ date | factor*as.factor(station.id),
#             data = tmp,
#             xlab = list(label = "Time", cex = 1.5),
#             ylab = list(label = "Temperature (degrees centigrade)", cex = 1.5),
#             type = "p",
#             layout = c(1,3),
#             aspect= "xy",
#             grib = TRUE,
##            xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
#             panel = function(...) {
#                  panel.abline(h=seq(0,max(tmp$tmax),by=10), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
#                  panel.xyplot(...)
#             }
#        )
#        print(b)
#dev.off()

##################################################################################
#QQ plot of Pooled deseasonal component of max temperature for 100 stations
##################################################################################

#Create the QQ plot of temperature conditional on month for one station
trellis.device(postscript, file = paste("QQ_plot_of_deseasonal_tmax_of_100_stations", ".ps", sep = ""), color=TRUE, paper="letter")
    for(i in stations){
        a <- qqmath(~ deseasonal,
                data = tmp[tmp[,1]==i,],
                distribution = qnorm,
                main = list(label= paste("QQ plot of max temperature for station ", i, sep=""), cex=2),
                xlab = list(label="Unit normal quantile", cex=1.5),
                ylab = list(label="Max temperature(degrees centigrade)", cex=1.5),
                scales = list(x = list(cex=1.5), y = list(cex=1.5)),
                prepanel = prepanel.qqmathline,
                panel = function(x, y,...) {
                        panel.grid()
                        panel.qqmathline(x, y=x)
                        panel.qqmath(x, y,...)
                }

        )
        print(a)
    }
dev.off()





