############################################################################
##Max temperature, 100 stations, seasonal+trend with raw data by stl2
##########################################################################

#load the data
library(lattice)
library(plyr)
library(stl2)
outputdir <- "~/Projects/Spatial/NCAR/output/"
datadir <- "~/Projects/Spatial/NCAR/RData/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))
load(paste(datadir,"stations.RData", sep=""))

dataset <- "tmax"
parameter1 <- list(sw=51, sd=1, tw=1141, td=1, fcw=1141, fcd=1, inner=10, outer=0)

if(dataset == "tmax"){
ylab <- "Maximum Temperature (degrees centigrade)"
}else if(dataset == "tmin"){
        ylab <- "Minimum Temperature (degrees centigrade)"
}else {
        ylab <- "Precipitation (millimeters)"
}
if(dataset %in% c("tmax", "tmin")){ 
	data <- UStemp
	datainfo <- UStinfo
}else{
	data <- USppt
	ddatainfo <- USpinfo
}
rm(list=grep("US", ls(), value=TRUE))
#find the 100 stations with largest observation number for max temperature
stations <- get(paste("stations", dataset, sep="."))
tmp <- data[with(data, station.id %in% stations),]
if(any(with(tmp, is.na(get(dataset))) == TRUE)) stop("The first 100 stations have NA")
month <- tmp$month
levels(month) <- c(1:12)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
column <- c("station.id", "year","month","date", dataset)
tmp <- tmp[with(tmp, order(station.id, date)), column]
names(tmp)[dim(tmp)[2]] <- "response"

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

dr <- ddply(tmp, "station.id", function(r){
	do.call("cbind", stl2(r$response, r$date, n.p=12, 
			s.window=parameter1$sw, s.degree = parameter1$sd, 
			t.window=parameter1$tw, t.degree = parameter1$td,
			fc.window=parameter1$fcw, fc.degree = parameter1$fcd,
			inner = parameter1$inner, outer = parameter1$outer)[c("data", "fc")])
})
tmp <- cbind(tmp, dr)
names(tmp)[grep("fc.fc", names(tmp))] <- "fc.trend"
tmp$ds <- tmp$response - tmp$data.seasonal
tmp$time <- rep(0:1235, 100)

info.100 <- datainfo[with(datainfo, station.id %in% stations), ]
info.100$state <- substring(info.100$station.id, 1,2)
info.100$id <- substring(info.100$station.id, 3)
stations <- as.character(info.100[with(info.100, order(state, elev)),]$station.id)
tmp$station.id <- factor(tmp$station.id, levels=stations)
##################################################
##Deseasonal component 
##################################################
trellis.device(postscript, file = paste(outputdir, dataset, "_with_stl2_deseasonal_conditional_month_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        a <- xyplot( ds ~ time | month,
             data = subset(tmp, station.id==i),
             xlab = list(label = "Month", cex = 1),
             ylab = list(label = ylab, cex = 1),
             main = list(label = paste("Station ", i, sep=""), cex=1),
             type = "b",
             pch = 16,
	     	 layout = c(2,6),
             cex = 0.5,
#             aspect = "xy",
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', alternating=TRUE), x=list(at=seq(0,1235, by=120), relation='same')),
             panel = function(...) {
#                  panel.abline(h=seq(10,40,by=1), v=seq(0,1236, by=120), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)

             }
        )
        print(a)
   }
dev.off()


trellis.device(postscript, file = paste(outputdir, dataset, "_deseasonal_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
#   for(i in stations){
        a <- xyplot( fc.remainder ~ time | station.id,
             data = tmp,
             xlab = list(label = "Month", cex = 1),
             ylab = list(label = ylab, cex = 1),
             type = "p",
             pch = 16,
			 layout = c(5,2),
             cex = 0.3,,
             aspect = "xy",
			 ylim = c(-1,1),
			 prepanel = function(x,y,...){
			 	prepanel.loess(x,y,degree=2, span=3/4, evaluation = 100, ...)
			 },
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list( relation='same', at=seq(0,1235,by=600))),
             panel = function(x,y,...) {
#                  panel.abline( v=seq(0,1235,by=24), color="lightgrey", lty=3, lwd=0.5)
#                  panel.xyplot(x,y,...)
		  		  panel.loess(x,y,span=3/4,degree=2, evaluation = 100, col="red",...)

             }
        )
        print(a)
#   }
dev.off()


trellis.device(postscript, file = paste(outputdir, "QQ_plot_", dataset, "_deseasonal_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        a <- qqmath( ~ ds | station.id,
             data = tmp,
			 distribution = qnorm,
             xlab = list(label="Unit normal quantile", cex=1.2),
             ylab = list(label = ylab, cex = 1),
             type = "p",
             pch = 16,
             layout = c(5,2),
             cex = 0.3,,
             aspect = 1,
             panel = function(x, y,...) {
             	panel.grid(lty=3, lwd=0.5, col="black",...)
                panel.qqmathline(x, y=x)
                panel.qqmath(x, y,...)
             }
        )
        print(a)
dev.off()



#fit the second seasonal component
dr.2 <- ddply(tmp, "station.id", function(r){
    do.call("cbind", loess(r$fc.remainder ~ r$time, span=2/3, degree=2)[c("fitted","residuals")] 
)})
dr.2$time <- rep(0:1235, 100)
#substract the second seasonal component from the deseasonal
#tmp$second.trend <- tmp$ds - dr.2$fitted
#for(z in 1){
#	if(z != 1){
#		tmp$second.trend <- tmp$ds - dr.4$fitted
#	}
#	#fit the second trend component
#	dr.3 <- ddply(tmp, "station.id", function(r){
#	    do.call("cbind", loess(r$second.trend ~ r$time, span=1141/1236, degree=1)[c("fitted","residuals")]
#	)})
#	dr.3$time <- rep(0:1235, 100)
#	tmp$second.seasonal <- tmp$ds - dr.3$fitted
#	dr.4 <- ddply(tmp, "station.id", function(r){
#	    do.call("cbind", loess(r$second.seasonal ~ r$time, span=2/3, degree=2)[c("fitted","residuals")]
#	)})
#	dr.4$time <- rep(0:1235, 100)
##dr.5 <- ddply(tmp, "station.id", function(r){
##    do.call("cbind", loess(r$fitted ~ r$time, span=3/4, degree=1)[c("fitted","residuals")]
##)})
#}
#tmp$second.trend <- dr.3$fitted
tmp$second.seasonal <- dr.2$fitted
tmp$second.remainder <- dr.2$residuals

trellis.device(postscript, file = paste(outputdir, dataset, "_secondremainder_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        a <- xyplot( second.remainder ~ time | station.id,
             data = tmp, 
             xlab = list(label = "Month", cex = 1),
             ylab = list(label = ylab, cex = 1),
			 layout = c(1,1),
#            aspect = "xy",
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', alternating=TRUE), x=list(at=seq(0,1235, by=240), relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=0, color="black", lty=1, lwd=0.5)
                  panel.xyplot(x,y,type="p", cex =0.3,...)
				  panel.loess(x,y,degree=1, span=0.35, col="red", type="l", lwd=0.5,...)
				  panel.loess(x,y,degree=2, span=0.15, col="green", type="l", lwd=0.5,...)
             }
        )
        print(a)
dev.off()

trellis.device(postscript, file = paste(outputdir, dataset, "_secondtrend_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        a <- xyplot( fc.trend ~ time | station.id,
             data = tmp,
             xlab = list(label = "Month", cex = 1),
             ylab = list(label = ylab, cex = 1),
             type = "p",
             pch = 16,
			 layout = c(4,3),
             cex = 0.3,
#			 ylim = c(-1,1),
#             aspect = "xy",
             scales = list(y = list(relation = 'free', alternating=TRUE), x=list(at=seq(0,1235, by=240), relation='same')),
#			 prepanel = function(x,y,...){
#			 	prepanel.loess(x,y,degree=1, span=3/4, evaluation = 100, ...)
#			 },
             panel = function(x,y,...) {
#                  panel.abline(h=seq(10,40,by=1), v=seq(0,1236, by=120), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
#				  panel.loess(x,y,span=3/4, degree=1, col="red",...)
             }
        )
        print(a)
dev.off()

trellis.device(postscript, file = paste(outputdir, dataset, "_secondseasonal_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
        a <- xyplot( second.seasonal ~ time | station.id,
             data = tmp,
             xlab = list(label = "Month", cex = 1),
             ylab = list(label = ylab, cex = 1),
             type = "p",
             pch = 16,
             layout = c(4,3),
             cex = 0.3,
#            ylim = c(-1,1),
#             aspect = "xy",
             scales = list(y = list(relation = 'free', alternating=TRUE), x=list(at=seq(0,1235, by=240), relation='same')),
#            prepanel = function(x,y,...){
#               prepanel.loess(x,y,degree=1, span=3/4, evaluation = 100, ...)
#            },
             panel = function(x,y,...) {
#                  panel.abline(h=seq(10,40,by=1), v=seq(0,1236, by=120), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
#                 panel.loess(x,y,span=3/4, degree=1, col="red",...)
             }
        )
        print(a)
dev.off()


ACF <- ddply(.data=tmp,
        .variables="station.id",
        .fun= summarise,
           correlation = c(acf(second.remainder, plot=FALSE)$acf),
           lag = c(acf(ds, plot=FALSE)$lag)
)

trellis.device(postscript, file = paste(outputdir, "acf_of_tmax_deseasonal_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        b <- xyplot( correlation ~ lag,
             data = subset(ACF, station.id==i & lag!=0),
             xlab = list(label = "Lag", cex = 1.2),
             ylab = list(label = paste("Station", i, "ACF"), cex = 1.2),
#             main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = "h",
             panel = function(x,y,...) {
                  panel.abline(h=0)
                  panel.xyplot(x,y,...)
             }
        )
        print(b)
   }
dev.off()

