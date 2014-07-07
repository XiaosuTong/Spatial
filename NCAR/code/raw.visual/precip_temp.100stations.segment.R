#############################################################
##Segment plot for each station to see the where are missing
#############################################################

#Set up the directory and load the data.
library(lattice)
datadir <- "~/Projects/Spatial/NCAR/RData/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))

#Count the observation number for each station.
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)
TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)

#Order the stations by the observation counts
#Pcount <- sort(Pcount, decreasing=TRUE)
#tmp.precip <- USppt
#tmp.precip$station.id <- factor(tmp.precip$station.id, levels = names(Pcount))
Pcount <- sort(Pcount, decreasing=TRUE)
stations <- tail(names(Pcount),10)

tmp.precip <- USppt[USppt[,1] %in% stations,]
#tmp.precip <- tmp.precip[!(is.na(tmp.precip$precip)),]

month <- tmp.precip$month
levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp.precip$year, month, "01", sep="-")
tmp.precip$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp.precip <- tmp.precip[order(tmp.precip$date),]
tmp.precip$time <- rep(1:1236,each=10)
tmp.precip$obs <- !is.na(tmp.precip$precip)

start <- as.POSIXct(strptime(paste(head(tmp.precip$year,1), "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
end <- as.POSIXct(strptime(paste("2000", "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")


trellis.device(postscript, file = paste("precipitation_measurement_status_against_month",".ps", sep = ""), color=TRUE, paper="legal")
        b <- dotplot( as.numeric(obs) ~ time | station.id,
             data = tmp.precip,
             xlab = list(label = "Month", cex = 1.5),
             ylab = list(label = "Measurement status", cex = 1.5),
             type = "p",
             col = c("black","red"),
             distribute.type= TRUE,
	     groups = obs,
             pch=16,
             cex=0.3,
             layout = c(1,10),
             strip.left = TRUE,
             strip = FALSE,
#             aspect= 0.1,
             grib = TRUE,
#             scales = list(y = list(relation = 'same', cex=1, labels=c("0","1"), alternating=TRUE), x=list(relation= 'same',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', cex=1, labels=c("0","1"), alternating=TRUE), x=list(tick.number = 10, relation='same', cex=1.2)),
             panel = function(x,y,...) {
                  panel.abline( v=seq(start, end, by="60 month"), color="lightgrey", lty=3, lwd=0.5)
		
		for(i in x){
		  panel.segments(x[i], y[i]-0.1, x[i], y[i]+0.1, lwd=0.5)
		}
                  panel.dotplot(x,y,...)
             }
        )
        print(b)
dev.off()

