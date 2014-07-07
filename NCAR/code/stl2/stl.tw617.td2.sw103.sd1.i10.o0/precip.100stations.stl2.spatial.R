#######################################################################################
##Precipitation, 100 stations, loess on each component from stl2 conditional on month
#######################################################################################

#load the data
library(lattice)
library(plyr)
library(stl2)
library(maps)
datadir <- "~/Projects/Spatial/NCAR/RData/"
outputdir <- "~/Projects/Spatial/NCAR/output/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)


#find the 100 stations with largest observation number for precipitation
Pcount <- sort(Pcount, decreasing=TRUE)
stations <- head(names(Pcount),100)
USppt$month <- factor(USppt$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))

tmp <- USppt[USppt[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$precip)),]

#Create date variable which is "POSIX1t" and "POSIXct" classes representing calender dates and time.
month <- tmp$month
levels(month) <- c(1:12)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$date),]

#Create start and end date for convenient
start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")

#stl2 for precipitation
dr <- ddply(tmp, "station.id", function(r){stl2(r$precip, r$date, n.p=12, s.window=65, t.window=151, inner = 5, outer = 0)$data})
tmp <- cbind(tmp, dr)
tmp <- tmp[order(tmp$station.id, tmp$month, tmp$year),]

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col


#################################################################
## seasonal component average cross years conditional on month
## then against spatial factors
#################################################################

tmp.seasonal <- tmp[, c(1:9,12)]

diff.seasonal <- ddply(.data = tmp.seasonal,
                        .variables = c("station.id","month"),
                        .fun = summarise,
                         lon = lon[1],
                         lat = lat[1],
                         elev = elev[1],
                         mean = mean(seasonal)
)

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_precip_seasonal_mean_vs_lat_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
       b <- xyplot( mean ~ lat | month,
             data = diff.seasonal,
             xlab = list(label = "Latitude (degree)", cex = 1.2),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.2),
             type = c("p"),
             pch=16,
             cex=0.5,
             layout = c(4,3),
             strip = TRUE,
#            aspect= "xy",
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
#                  panel.abline(h=seq(-20,20,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
             }
        )
        print(b)
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_precip_seasonal_mean_vs_lon_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
       b <- xyplot( mean ~ lon | month,
             data = diff.seasonal,
             xlab = list(label = "Longitude (degree)", cex = 1.2),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.2),
             type = c("p"),
             pch=16,
             cex=0.5,
             layout = c(4,3),
             strip = TRUE,
#            aspect= "xy",
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
#                  panel.abline(h=seq(-20,20,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
             }
        )
        print(b)
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_precip_seasonal_mean_vs_elev_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
       b <- xyplot( mean ~ log2(elev) | month,
             data = diff.seasonal,
             xlab = list(label = "Log of Elevation (log base 2 meter)", cex = 1.2),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.2),
             type = c("p"),
             pch=16,
             cex=0.5,
             layout = c(4,3),
             strip = TRUE,
#            aspect= "xy",
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
#                  panel.abline(h=seq(-20,20,by=5), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
             }
        )
        print(b)
dev.off()


