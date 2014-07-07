###################################################################################################
##Plot the precipitation and temperature for the 100 most observations stations conditional on month
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

#####################################################
##scatter plot for precipitation conditional on month
#####################################################
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

dd <- ddply(.data=tmp,
            .variables = c("station.id","month"),
            .fun = summarise,
             mean = mean(precip))
mm <- dd[rep(row.names(dd), each=103),]
tmp$precip <- tmp$precip - mm$mean


#dr <- ddply(.data=tmp,
#	    .variables= c("station.id","month"), 
#	    .fun = summarise,
#	     year = c(1895:1997),
#	     fitted = loess(precip~year, degree=1, span=2/3)$fitted)
#
#tmp <- cbind(tmp, loess=dr$fitted)
#tmp1 <- tmp[,c(1:7, 9,8)]
#tmp2 <- tmp[,c(1:7,9,10)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp <- rbind(tmp2, tmp1)
#tmp$group <- rep(c("loess","raw"), each=123600)


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


#########################################################
##scatter plot for maximum temperature conditional on month
#########################################################
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
start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$station.id, tmp$month, tmp$year),]

dd <- ddply(.data=tmp,
	    .variables = c("station.id","month"),
	    .fun = summarise,
	     mean = mean(tmax))
mm <- dd[rep(row.names(dd), each=103),]
tmp$tmax <- tmp$tmax - mm$mean

#dr <- ddply(.data=tmp,
#            .variables= c("station.id","month"),
#            .fun = summarise,
#             year = c(1895:1997),
#             fitted = loess(tmax~ year, degree = 1, span = 2/3)$fitted)
#
#tmp <- cbind(tmp, loess=dr$fitted)
#tmp1 <- tmp[,c(1:7, 10,9)]
#tmp2 <- tmp[,c(1:7,10,11)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp <- rbind(tmp2, tmp1)
#tmp$group <- rep(c("loess","raw"), each=123600)


trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmax_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        b <- xyplot( tmax ~ year | month,
             data = subset(tmp,station.id==i),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.2),
             main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = c("p"),
             pch=16,
             cex=0.5,	     
             layout = c(4,3),
	     strip = TRUE,
#             aspect= 1,
             grib = TRUE,
#             xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-10,10,by=2), v=seq(1900,2000,by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
		  panel.loess(x,y,span=2/3,degree=1,col=col[2],...)
             }
        )
        print(b)
    }
dev.off()

########################################################
##scatter plot for minimum temperature conditional on month
#######################################################
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
start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$station.id, tmp$month, tmp$year),]
#tmp$factor <- factor(rep(rep(paste("Period", 1:9), c(rep(137,8),140)), times=100), levels=paste("Period", c(9:1)))

dd <- ddply(.data=tmp,
            .variables = c("station.id","month"),
            .fun = summarise,
             mean = mean(tmin))
mm <- dd[rep(row.names(dd), each=103),]
tmp$tmin <- tmp$tmin - mm$mean


#dr <- ddply(.data=tmp,
#            .variables= c("station.id","month"),
#            .fun = summarise,
#             year = c(1895:1997),
#             fitted = loess(tmin~year, degree=1, span=2/3)$fitted)
#tmp <- cbind(tmp, loess=dr$fitted)
#tmp1 <- tmp[,c(1:7,10,8)]
#tmp2 <- tmp[,c(1:7,10,11)]
#names(tmp1)[9] <- names(tmp2)[9] <- "response"
#tmp <- rbind(tmp2, tmp1)
#tmp$group <- rep(c("loess","raw"), each=123600)


trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmin_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
    for(i in stations){
        b <- xyplot( tmin ~ year | month,
             data = subset(tmp,station.id==i),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = "Minimum Temperature (degrees centigrade)", cex = 1.2),
	     main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = c("p"),
	     strip = TRUE,
             pch=16,
             cex=0.5,
             layout = c(4,3),
#             aspect= "xy",
             grib = TRUE,
#             xlim = c(start, end),
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=seq(-10,10,by=2), v=seq(1900,2000, by=10), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
		  panel.loess(x,y,span=2/3,degree=1,col=col[2],...)
             }
        )
        print(b)
    }
dev.off()

