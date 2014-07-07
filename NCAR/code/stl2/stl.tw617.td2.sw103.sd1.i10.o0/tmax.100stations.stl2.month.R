#######################################################################################
##Max temperature, 100 stations, loess on each component from stl2 conditional on month
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

#find the 100 stations with largest observation number for max temperature
stations <- stations.tmax
tmp <- UStemp[UStemp[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$tmax)),]

#Create date variabel which is "POSIX1t" and "POSIXct" classes representing calender dates and time.
month <- tmp$month
levels(month) <- c(1:12)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp <- tmp[order(tmp$date),]
tmp <- tmp[order(tmp[,1]),-8]

#Create start and end date for convenient.
start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")

#stl2 for max temperature.
dr <- ddply(tmp, "station.id", function(r){stl2(r$tmax, r$date, n.p=12, t.window = 617, t.degree=2, s.window=103, s.degree=1, inner = 10, outer = 0)$data})
tmp <- cbind(tmp, dr)
tmp <- tmp[order(tmp$station.id, tmp$month, tmp$year),]
tmp$year <- tmp$year - 1894
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

###########################################
##Seasonal component conditional on month
###########################################

tmp.seasonal <- tmp[, c(1:9,12)]

#order the 100 stations of Max temperature by the elevation.
#order.elev <- tmp.seasonal[order(tmp.seasonal$year, tmp.seasonal$month),][1:100,1:2]
#order.elev <- order.elev[order(order.elev$elev),1]
#elev.od <- as.character(order.elev)

#Instead of calculating loess estimate manually, we just call panel.loess function in xyplot to get the loess curve.
#dd <- ddply(.data=tmp.seasonal,
#            .variables = c("station.id","month"),
#            .fun = summarise,
#             mean = mean(seasonal))
#mm <- dd[rep(row.names(dd), each=103),]
#tmp.seasonal$seasonal <- tmp.seasonal$seasonal - mm$mean

#seasonal.dr <- ddply(.data=tmp.seasonal,
#            	.variables= c("station.id","month"),
#            	.fun = summarise,
#            	 year = c(1895:1997),
#            	 fitted = loess(seasonal~year, span=2/3, degree=1)$fitted)
#
#tmp.seasonal <- cbind(tmp.seasonal, loess=seasonal.dr$fitted)

#tmp1 <- tmp.seasonal[,c(1:10)]
#tmp2 <- tmp.seasonal[,c(1:9,11)]
#names(tmp1)[10] <- names(tmp2)[10] <- "response"
#tmp.seasonal <- rbind(tmp1, tmp2)
#tmp.seasonal$group <- rep(c("seasonal","loess"), each=123600)

#Attend to add a small map in the corner of the plot.
#us.map <- map('state', plot = FALSE, fill = TRUE)

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmax_seasonal_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
        b <- xyplot( seasonal ~ year | month,
             data = subset(tmp.seasonal,station.id==i),
#             main = list(label= paste("Station ", i, sep=""), cex=1.5),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = paste("Station", i, "Maximum Temperature (degrees centigrade)"), cex = 1.2),
             type = c("p"),
             pch=16,
             cex=0.5,
             layout = c(12,1),
             strip = TRUE,
#             aspect= 5,
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=0, color="black", lty=1)
                  panel.xyplot(x,y,...)
             }
        )
	print(b)
#        a <- xyplot(lat ~ lon,
#                data = subset(tmp.seasonal, station.id==i),
#                xlab = NULL,
#                ylab = NULL,
#                pch = 16,
#                cex = 0.4,
#                xlim = c(-127,-65),
#                ylim = c(23,50),
#                col = "red",
##               scales = list(x = list(draw=FALSE), y = list(draw=FALSE)),
##               strip = strip.custom(par.strip.text= list(cex = 1.5)),
##               par.settings = list(layout.heights = list(strip = 1.5)),
#                panel = function(x,y,...) {
#                        panel.polygon(us.map$x,us.map$y,lwd=0.2)
#                        panel.xyplot(x,y,...)
#                }
#        )
#        plot.new()
#        title(paste("Station ", i, sep=""), cex=1.5)
#        print(b, pos=c(0,0,1,0.65), newpage=FALSE, more=TRUE)
#        print(a, pos=c(0,0.6,0.4,1), more=FALSE)
     }
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmax_seasonal2_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
        b <- xyplot( seasonal ~ year | month,
             data = subset(tmp.seasonal,station.id==i),
#             main = list(label= paste("Station ", i, sep=""), cex=1.5),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = paste("Station", i, "Maximum Temperature (degrees centigrade)"), cex = 1.2),
             type = c("p"),
             pch=16,
             cex=0.5,
             layout = c(2,6),
             strip = TRUE,
#             aspect= 5,
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'free', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
#                  panel.abline(h=0, color="black", lty=1)
                  panel.xyplot(x,y,...)
             }
        )
        print(b)
     }
dev.off()

rm(tmp.seasonal)


############################################
##Remainder component conditional on month
############################################
tmp.remainder <- tmp[, c(1:9,14)]

#remainder.dr <- ddply(.data=tmp.remainder,
#                .variables= c("station.id","month"),
#                .fun = summarise,
#                 year = c(1895:1997),
#                 fitted = loess(remainder~year, span=2/3, degree=1)$fitted)
#
#tmp.remainder <- cbind(tmp.remainder, loess=remainder.dr$fitted)
#
#tmp1 <- tmp.remainder[,c(1:10)]
#tmp2 <- tmp.remainder[,c(1:9,11)]
#names(tmp1)[10] <- names(tmp2)[10] <- "response"
#tmp.remainder <- rbind(tmp1, tmp2)
#tmp.remainder$group <- rep(c("remaind","loess"), each=123600)

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmax_remainder_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
        b <- xyplot( remainder ~ year | month,
             data = subset(tmp.remainder,station.id==i),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = paste("Station", i, "Maximum Temperature (degrees centigrade)"), cex = 1.2),
#             main = list(label = paste("Station ", i, sep=""), cex=1.5),
             pch=16,
             cex=0.5,
             layout = c(12,1),
             strip = TRUE,
#             aspect= 5,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=0, color="black", lty=1)
                  panel.xyplot(x,y,...)
		  panel.loess(x,y,span=3/4, degree=1, col=col[2],...)
             }
        )
        print(b)
     }
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmax_remainder2_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
        b <- xyplot( remainder ~ year | month,
             data = subset(tmp.remainder,station.id==i),
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = paste("Station",i,"Maximum Temperature (degrees centigrade)"), cex = 1.2),
#             main = list(label = paste("Station ", i, sep=""), cex=1.5),
	     type= "b",	     
             pch=16,
             cex=0.5,
             layout = c(2,6),
             strip = TRUE,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(x,y,...) {
                  panel.abline(h=0, color="black", lty=1)
                  panel.xyplot(x,y,...)
             }
        )
        print(b)
     }
dev.off()


trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_tmax_remainder_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
    for(i in stations){
       a <- qqmath(~ remainder | month,
                data = subset(tmp.remainder,station.id==i),
                distribution = qnorm,
                aspect = 1,
                pch = 16,
                cex = 0.5,
                layout = c(4,3),
                main = list(label= paste("Station ", i, sep=""), cex=1.5),
                xlab = list(label="Unit normal quantile", cex=1.2),
                ylab = list(label="Maximum Temperature(degrees centigrade)", cex=1.2),
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

tmp.sr <- tmp[, c(1:9,12,14)]
tmp.sr$sr <- tmp.sr$seasonal + tmp.sr$remainder

#just the original seasonal and seasonal+remainder
tmp1 <- tmp.sr[,c(1:10)]
tmp2 <- tmp.sr[,c(1:9,12)]
names(tmp1)[10] <- names(tmp2)[10] <- "response"
tmp.sr.raw <- rbind(tmp1, tmp2)
tmp.sr.raw$group <- rep(c("seasonal","seasonal+remainder"), each=123600)

#substruct the mean of the monthly seasonal from seasonal and seasonal+remainder
dd <- ddply(.data=tmp.sr,
            .variables = c("station.id", "month"),
            .fun = summarise,
             mean = mean(seasonal))
mm <- dd[rep(row.names(dd), each=103),]
mm <- mm[order(mm$station.id, mm$month),]
tmp.sr.mean <- tmp.sr
tmp.sr.mean$seasonal <- tmp.sr.mean$seasonal - mm$mean
tmp.sr.mean$sr <- tmp.sr.mean$sr - mm$mean

tmp1 <- tmp.sr.mean[,c(1:10)]
tmp2 <- tmp.sr.mean[,c(1:9,12)]
names(tmp1)[10] <- names(tmp2)[10] <- "response"
tmp.sr.mean <- rbind(tmp1, tmp2)
tmp.sr.mean$group <- rep(c("seasonal","seasonal+remainder"), each=123600)


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
trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmax_seasonal+remainder_nomean_same_for_100_stations_conditional_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
        b <- xyplot( response ~ year | month,
             data = subset(tmp.sr.mean, station.id==i),
             groups = group,
             xlab = list(label = "Year", cex = 1.2),
             ylab = list(label = paste("Station",i,"Maximum Temperature (degrees centigrade)"), cex = 1.2),
#             main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = c("l","p"),
             col = col[2:1],
             distribute.type= TRUE,
             pch=16,
             cex=0.5,
             key=list(type="p", text=list(label=c("seasonal","seasonal+remainder")), points=list(cex=1, col=col[2:1], pch=16), columns=2),
             layout = c(12,1),
             strip = TRUE,
#             aspect= 5,
#             scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same')),
             panel = function(...) {
                  panel.abline(h=0, color="black", lty=1)
                  panel.xyplot(...)
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


