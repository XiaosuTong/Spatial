#######################################################################################
##100 stations, loess on each component from stl2 conditional on month
#######################################################################################

#load the data
library(lattice)
library(plyr)
library(stl2)
library(maps)
datadir <- "~/Projects/Spatial/NCAR/RData/"
outputdir <- "~/Projects/Spatial/NCAR/output/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))
load(paste(datadir,"stations.RData", sep=""))

dataset <- "tmax"
#parameter <- list(sw="periodic", sd=1, tw=1141, td=1, inner=10, outer=0, flag=FALSE)
#parameter <- list(sw="periodic", sd=1, tw=1855, td=1, fcw=c(1855,121), fcd=c(1,2), inner=10, outer=0, flag=TRUE)
parameter <- list(sw="periodic", sd=1, tw=1855, td=1, fcw=c(1855,241), fcd=c(1,1), inner=10, outer=0, flag=TRUE)

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
tmp <- data[data[,1] %in% stations,]
if(any(with(tmp, is.na(get(dataset))) == TRUE)) stop("The first 100 stations have NA")

#Create date variabel which is "POSIX1t" and "POSIXct" classes representing calender dates and time.
month <- tmp$month
levels(month) <- c(1:12)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
column <- c("station.id", "year","month","date", dataset)
tmp <- tmp[with(tmp, order(station.id, date)), column]
names(tmp)[dim(tmp)[2]] <- "response"

#stl2
if(parameter$flag){
dr <- ddply(tmp, "station.id", function(r){
    do.call("cbind", stl2(r$response, r$date,
    n.p=12,
    s.window = parameter$sw,
    s.degree = parameter$sd,
    t.window = parameter$tw,
    t.degree = parameter$td,
    fc.window = parameter$fcw,
    fc.degree = parameter$fcd,
    inner = parameter$inner,
    outer = parameter$outer)[c("data","fc")])
})
tmp <- cbind(tmp, dr[, c(!(names(dr) %in% c("station.id", "data.raw", "data.sub.labels")))])
names(tmp)[grep("fc.fc", names(tmp))] <- c("fc.trend", "fc.second")
}else{
dr <- ddply(tmp, "station.id", function(r){
    stl2(r$response, r$date,
    n.p = 12,
    s.window = parameter$sw,
    s.degree = parameter$sd,
    t.window = parameter$tw,
    t.degree = parameter$td,
    inner = parameter$inner,
    outer = parameter$outer)$data
})
tmp <- cbind(tmp, dr[, c(!(names(dr) %in% c("station.id", "raw", "sub.labels")))])
}
tmp$year <- tmp$year - 1894

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

###########################################
##Seasonal component conditional on month
###########################################
tmp.seasonal <- tmp[, c(!(names(tmp) %in% c("data.weights","data.trend","data.remainder")))]

order <- ddply(.data = tmp.seasonal,
                .variables = "station.id",
                .fun = summarise,
                 mean = mean(response)
)
order.st <- as.character(order[order(order$mean, decreasing=TRUE), ]$station.id)
tmp.seasonal$station.id <- factor(tmp.seasonal$station.id, levels=order.st)

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_seasonal_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in order.st){
        b <- xyplot( data.seasonal ~ year | month,
             data = subset(tmp.seasonal,station.id==i),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = paste("Station", i, ylab, sep=" "), cex = 1.2),
             type = c("p"),
             pch=16,
             cex=0.5,
             layout = c(12,1),
             strip = TRUE,
#            scales = list(y = list(relation='free', cex=1.5), x=list(relation='free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=0, color="black", lty=1)
                  panel.xyplot(x,y,...)
             }
        )
	print(b)
     }
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_seasonal2_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
        b <- xyplot( data.seasonal ~ month | station.id,
             data = subset(tmp.seasonal,year==1),
             xlab = list(label = "Month", cex = 1.2),
             ylab = list(label = ylab, cex = 1.2),
             type = "b",
             pch=16,
			 aspect = "xy",
             cex=0.5,
             layout = c(5,4),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(at=c(1, 3, 5, 7, 9, 11), relation='same')),
             panel = function(x,y,...) {
                  panel.xyplot(x,y,...)
             }
        )
        print(b)
dev.off()

rm(tmp.seasonal)


############################################
##Remainder component conditional on month
############################################
tmp.remainder <- tmp[, c(!(names(tmp) %in% c("data.weights","data.trend","data.seasonal")))]

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_remainder_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
        b <- xyplot( fc.remainder ~ year | month,
             data = subset(tmp.remainder,station.id==i),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = paste("Station", i, ylab, sep=" "), cex = 1.2),
             pch=16,
             cex=0.5,
             layout = c(12,1),
             strip = TRUE,
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...){
             	panel.abline(h=0, color="black", lty=1)
             	panel.xyplot(x,y,...)
		     	panel.loess(x,y,span=3/4, degree=1, col=col[2],...)
             }
        )
        print(b)
     }
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_remainder2_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
        b <- xyplot( fc.remainder ~ year | month,
             data = subset(tmp.remainder,station.id==i),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = paste("Station",i,ylab, sep=" "), cex = 1.2),
	       	 type= "b",	     
             pch=16,
             cex=0.5,
             layout = c(2,6),
             strip = TRUE,
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=0, color="black", lty=1)
                  panel.xyplot(x,y,...)
             }
        )
        print(b)
     }
dev.off()


trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_", dataset, "_remainder_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
    for(i in stations){
       a <- qqmath(~ fc.remainder | month,
                data = subset(tmp.remainder,station.id==i),
                distribution = qnorm,
                aspect = 1,
                pch = 16,
                cex = 0.3,
                layout = c(4,3),
                xlab = list(label="Unit normal quantile", cex=1.2),
                ylab = list(label = paste("Station",i,ylab, sep=" "), cex = 1.2),
                prepanel = prepanel.qqmathline,
                panel = function(x, y,...) {
                        panel.grid(lty=3, lwd=0.5, col="black",...)
                        panel.qqmathline(x, y=x)
                        panel.qqmath(x, y,...)
                }
        )
        print(a)
    }
dev.off()


rm(tmp.remainder)

###########################################
##Seasonal+remainder conditional on month
##########################################
tmp.sr <- tmp[, c(!(names(tmp) %in% c("weights","trend")))]
tmp.sr$sr <- tmp.sr$seasonal + tmp.sr$remainder

#just the original seasonal and seasonal+remainder
#tmp1 <- tmp.sr[, c(!(names(tmp.sr) %in% c("sr","remainder")))]
#tmp2 <- tmp.sr[, c(!(names(tmp.sr) %in% c("seasonal","remainder")))]
#names(tmp1)[dim(tmp1)[2]] <- names(tmp2)[dim(tmp2)[2]] <- "resp"
#tmp.sr.raw <- rbind(tmp1, tmp2)
#tmp.sr.raw$group <- rep(c("seasonal","seasonal+remainder"), each=123600)

#substruct the mean of the monthly seasonal from seasonal and seasonal+remainder
#dd <- ddply(.data=tmp.sr,
#            .variables = c("station.id", "month"),
#            .fun = summarise,
#             mean = mean(seasonal))
#mm <- dd[rep(row.names(dd), each=103),]
#mm$year <- rep(1:103, 120)
#mm <- mm[with(mm, order(station.id, year)), ]
#tmp.sr.mean <- tmp.sr
#tmp.sr.mean$seasonal <- tmp.sr.mean$seasonal - mm$mean
#tmp.sr.mean$sr <- tmp.sr.mean$sr - mm$mean

#tmp1 <- tmp.sr.mean[, c(!(names(tmp.sr) %in% c("sr","remainder")))]
#tmp2 <- tmp.sr.mean[, c(!(names(tmp.sr) %in% c("seasonal","remainder")))]
#names(tmp1)[dim(tmp1)[2]] <- names(tmp2)[dim(tmp2)[2]] <- "resp"
#tmp.sr.mean <- rbind(tmp1, tmp2)
#tmp.sr.mean$group <- rep(c("seasonal","seasonal+remainder"), each=123600)


#with mean with same scale
#trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmax_seasonal+remainder_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
#     for(i in stations){
#        b <- xyplot( response ~ year | month,
#             data = subset(tmp.sr.raw, station.id==i),
#             groups = group,
#             xlab = list(label = "Year", cex = 1.2),
#             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.2),
#             main = list(label = paste("Station ", i, sep=""), cex=1.5),
#             type = c("l","p"),
#             col = col[2:1],
#             distribute.type= TRUE,
#             pch=16,
#             cex=0.5,
#	     key=list(type="p", text=list(label=c("seasonal","seasonal+remainder")), points=list(cex=1, col=col[2:1], pch=16), columns=2),
#             layout = c(6,1),
#             strip = TRUE,
##             aspect= "xy",
#             grib = TRUE,
##             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
#             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
#             panel = function(...) {
#                  panel.abline(h=seq(-20,20,by=5), v=seq(1900,2000,by=10), color="lightgrey", lty=3, lwd=0.5)
#                  panel.xyplot(...)
#             }
#        )
#        print(b)
#     }
#dev.off()

#no mean with same scale
trellis.device(postscript, file = paste(outputdir, "scatterplot_of_", dataset, "_seasonal+remainder_nomean_same_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
        b <- xyplot( sr ~ year | month,
             data = subset(tmp.sr, station.id==i),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = paste("Station",i, ylab, sep=" "), cex = 1.2),
             key=list(text=list(label=c("seasonal","seasonal+remainder")), lines=list(pch=c("","."), cex=4, lwd=1.5, type=c("l","p"), col=col[2:1]), columns=2),
             layout = c(12,1),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
			 prepanel = function(x,y,subscripts,...){
			    v <- subset(tmp.sr, station.id==i)[subscripts,]
                mean <- v$seasonal
				prepanel.default.xyplot(x,y-mean,...)
			 },
             panel = function(x,y,subscripts,...) {
             	panel.abline(h=0, color="black", lty=1)
			    v <- subset(tmp.sr, station.id==i)[subscripts,]
			 	mean <- v$seasonal
			 	panel.xyplot(x, y-mean, type="p", col=col[1], pch=16, cex=0.5, ...)
			 	panel.xyplot(x, v$seasonal-mean, type="l", col=col[2], ...)
             }
        )
        print(b)
     }
dev.off()

#no mean with free scale
#trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmax_seasonal+remainder_nomean_free_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
#     for(i in stations){
#        b <- xyplot( response ~ year | month,
#             data = subset(tmp.sr.mean, station.id==i),
#             groups = group,
#             xlab = list(label = "Year", cex = 1.2),
#             ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.2),
#             main = list(label = paste("Station ", i, sep=""), cex=1.5),
#             type = c("l","p"),
#             col = col[2:1],
#             distribute.type= TRUE,
#             pch=16,
#             cex=0.5,
#             key=list(type=c("p"), text=list(label=c("seasonal","seasonal+remainder")), points=list(cex=1, col=col[2:1], pch=16), columns=2),
#             layout = c(6,1),
#             strip = TRUE,
##             aspect= "xy",
#             grib = TRUE, 
##             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
#             scales = list(y = list(relation = 'free', alternating=TRUE), x=list(tick.number=10, relation='same')),
#             panel = function(...) {
#                  panel.abline(h=seq(-20,20,by=5), v=seq(1900,2000,by=10), color="lightgrey", lty=3, lwd=0.5)
#                  panel.xyplot(...)
#             }
#        )
#        print(b)
#     }
#dev.off()


