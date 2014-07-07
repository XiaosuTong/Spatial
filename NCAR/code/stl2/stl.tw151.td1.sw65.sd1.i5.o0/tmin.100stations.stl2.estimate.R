############################################################################
##Max temperature, 100 stations, after stl2, study about spatial property
##########################################################################

#load the data
library(lattice)
library(plyr)
library(stl2)
datadir <- "~/Projects/Spatial/NCAR/RData/"
#load(paste(datadir,"USmonthlyMet.RData", sep=""))
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)


#find the 100 stations with largest observation number for max temperature
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
tmp$factor <- factor(rep(rep(paste("Period", 1:9), c(rep(137,8),140)), times=100), levels=paste("Period", c(9:1)))


dr <- ddply(tmp, "station.id", function(r){stl2(r$tmin, r$date, n.p=12, s.window=15, inner = 5, outer = 1)$data})
tmp <- cbind(tmp, dr)

tmp <- tmp[order(tmp$station.id,tmp$month),]
tmp <- tmp[,c(1:10,13,14,15)]

estimate.min.QQ <- function(tmp, station, dist){

##Find one target station which will be predicted by using closed stations.
##unique(tmp[tmp$lon< -80&tmp$lon>-85&tmp$lat>38&tmp$lat<43,]$station.id)
##The station for this condition is 120200 with lon=-84.98 lat=41.63
        if(!(station %in% tmp$station.id))
        stop("The station is not in the 100 stations")
	
	longitude <- tmp[tmp$station.id == station,]$lon[1]
	latitude <- tmp[tmp$station.id == station,]$lat[1]
	tmp$distance <- sqrt((tmp$lon - longitude)^2 + (tmp$lat - latitude)^2)

	a <- cbind(station.id=as.character(tmp$station.id), station.name=as.character(tmp$station.name),distance=tmp$distance)
	rst <- as.data.frame(a[seq(1,123600,by=1236),], stringsAsFactors=FALSE)
	rst$distance <- as.numeric(rst$distance)
	rst <- rst[order(rst$distance),]
	thestations <- rst$station.id[rst$distance > 0 &rst$distance < dist]

        while(length(thestations)==0){
                dist <- dist + 1
                thestations <- rst$station.id[rst$distance > 0 &rst$distance < dist]
        }


	rst <- tmp[tmp$station.id %in% thestations,]
	rst <- rst[order(rst$date),]

	s <- ddply(rst, "date", function(r){c(est.seasonal = weighted.mean(r$seasonal, r$distance))})
	t <- ddply(rst, "date", function(r){c(est.trend = weighted.mean(r$trend, r$distance))})
	temp <- as.data.frame(cbind(s, est.trend=t$est.trend))
	data <- tmp[tmp$station.id == station,-14]
	temp <- cbind(data[order(data$date),], temp[,-1])
	remainder.est <- temp$tmin - temp$est.seasonal - temp$est.trend

       	a <- qqmath(~ remainder.est,
                distribution = qnorm,
                main = list(label= paste("QQ plot of remainder of min temperature for station ", station, sep=""), cex=2),
                xlab = list(label="Unit normal quantile", cex=1.5),
                ylab = list(label="Min temperature(degrees centigrade)", cex=1.5),
                prepanel = prepanel.qqmathline,
                panel = function(x, y,...) {
                        panel.grid()
                        panel.qqmathline(x, y=x)
                        panel.qqmath(x, y,...)
                }

        )
	print(a)
}


estimate.min.scatter <- function(tmp, station, dist){
	
	if(!(station %in% tmp$station.id))
	stop("The station is not in the 100 stations")

        longitude <- tmp[tmp$station.id == station,]$lon[1]
        latitude <- tmp[tmp$station.id == station,]$lat[1]
        tmp$distance <- sqrt((tmp$lon - longitude)^2 + (tmp$lat - latitude)^2)
	
        a <- cbind(station.id=as.character(tmp$station.id), station.name=as.character(tmp$station.name),distance=tmp$distance)
        rst <- as.data.frame(a[seq(1,123600,by=1236),], stringsAsFactors=FALSE)
        rst$distance <- as.numeric(rst$distance)
        rst <- rst[order(rst$distance),]
        thestations <- rst$station.id[rst$distance > 0 &rst$distance < dist]
	
	while(length(thestations)==0){
		dist <- dist + 1
		thestations <- rst$station.id[rst$distance > 0 &rst$distance < dist]
	}

        rst <- tmp[tmp$station.id %in% thestations,]	
        rst <- rst[order(rst$date),]
        
	s <- ddply(rst, "date", function(r){c(est.seasonal = weighted.mean(r$seasonal, r$distance))})
        t <- ddply(rst, "date", function(r){c(est.trend = weighted.mean(r$trend, r$distance))})
        temp <- as.data.frame(cbind(s, est.trend=t$est.trend))
        data <- tmp[tmp$station.id == station,-14]
        temp <- cbind(data[order(data$date),], temp[,-1])
	temp1 <- temp[c(1:7,9,10,8)]
	names(temp1)[10] <- "response"
	temp1$group <- rep("raw", 1236)
	temp2 <- temp[,c(1:7,9,10,14,15)]
	temp2$response <- temp2$est.seasonal + temp2$est.trend
	temp2 <- temp2[, -c(10,11)]
	temp2$group <- rep("estimate", 1236)
	temp <- rbind(temp1, temp2)
	
	start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
	end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")


	b <- xyplot( response ~ date | factor,
             data = temp,
             groups = group,
             xlab = list(label = "Time", cex = 1.5),
             ylab = list(label = "Temperature (degrees centigrade)", cex = 1.5),
             main = list(label = paste("Min Temperature for station ", station, sep=""), cex=1.5),
             type = c("p","l"),
             distribute.type= TRUE,
             layout = c(1,9),
             aspect= 0.06,
             strip = FALSE,
             strip.left = TRUE,
             grib = TRUE,
             scales = list(y = list(relation = 'same', cex=1, alternating=TRUE), x=list(draw=FALSE, relation='free')),
             panel = function(...) {
                  panel.abline(h=seq(min(temp$response), max(temp$response),by=5), v=seq(start, end, by="12 month"), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)

             }
        )
        print(b)
	return(dist)

}

trellis.device(postscript, file = paste("scatterplot_of_tmin_with_stl2_estimate_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
    for(i in unique(tmp$station.id)){
	estimate.min.scatter(tmp, i, 4)
    }
dev.off()

trellis.device(postscript, file = paste("QQ_plot_of_stl2_estimate_emainder_tmin_of_100_stations", ".ps", sep = ""), color=TRUE, paper="letter")
    for(i in unique(tmp$station.id)){
	estimate.min.QQ(tmp, i, 4)
    }
dev.off()
