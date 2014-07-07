####################################################
##Check the time series plot for each station
####################################################

#load the data
library(lattice)
datadir <- "~/Project/Spatial/NCAR/RData/"
#load(paste(datadir,"USmonthlyMet.RData", sep=""))

#find the tempreture data for one station

timeseries <- function(USppt, UStemp, station){

	tmp <- UStemp[UStemp[,1] == station,]
	tmp <- tmp[!(is.na(tmp$tmax)&is.na(tmp$tmin)),]


	month <- tmp$month
	levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
	month <- as.numeric(factor(month, levels=c(1:12)))
	date <- paste(tmp$year, month, "01", sep="-")
	tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
	tmp <- tmp[order(tmp$date),]
	start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
	end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")

	tmp1 <- tmp[,-8]
	tmp1$group <- rep("tmax", nrow(tmp))
	tmp2 <- tmp[,-9]
	tmp2$group <- rep("tmin", nrow(tmp))
	names(tmp1)[8] <- "temp"
	names(tmp2)[8] <- "temp"
	tmp <- rbind(tmp1, tmp2)

	lattice.theme <- trellis.par.get()
	col <- lattice.theme$superpose.symbol$col
	trellis.device(postscript, file = paste("lineplot_of_tempreture_for_", station, ".ps", sep = ""), color=TRUE, paper="legal")
	b <- xyplot( temp ~ date,
             data = tmp,
	     lwd = 2,
	     group = group,
	     xlab = list(label = "Time", cex = 1.5),
	     ylab = list(label = "Tempreture (degrees centigrade)", cex = 1.5),
             type = "l",
             grib = TRUE,
             xlim = c(start, end),
             scales = list(y = list(relation = 'free', cex=1.5), x=list(format = "%b %Y", tick.number=10), cex=1.2),
             key = list(columns = 3, text = list(c("Max", "Min"), cex = 1.2), lines=list(lwd = 3, col = col[1:2])),	     
	     panel = function(...) {
                  panel.abline(h=seq(0,40,by=5), v=seq(start, end, by="48 month"), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
	)
	print(b)
	dev.off()


	#find the precipitation data for one station
	tmp <- USppt[USppt[,1] == station,]
	tmp <- tmp[!(is.na(tmp$precip)),]

	month <- tmp$month
	levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
	month <- as.numeric(factor(month, levels=c(1:12)))
	date <- paste(tmp$year, month, "01", sep="-")
	tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
	tmp <- tmp[order(tmp$date),]
	start <- as.POSIXct(strptime(paste(head(tmp$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
	end <- as.POSIXct(strptime(paste(tail(tmp$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")

        trellis.device(postscript, file = paste("lineplot_of_precipitation_for_", station, ".ps", sep = ""), color=TRUE, paper="legal")
	b <- xyplot( precip ~ date,
             data = tmp,
	     lwd = 2,
             xlab = list(label = "Time", cex = 1.5),
             ylab = list(label = "Precipitation (millimeters)", cex = 1.5),
             type = "l",
             grib = TRUE,
             xlim = c(start, end),
             scales = list(y = list(relation = 'free', cex=1.5), x=list(format = "%b %Y", tick.number=10), cex=1.2),
             panel = function(...) {
                  panel.abline(h=seq(0,max(tmp$precip),by=5), v=seq(start, end, by="48 month"), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(...)
             }
	)
	print(b)
	dev.off()
}



