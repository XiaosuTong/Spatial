###############################################################
##Count the observation number in precipitation and tempreture
##############################################################

#Set up the directory and load the data.
library(lattice)
library(plyr)
library(maps)
datadir <- "~/Projects/Spatial/NCAR/RData/"
outputdir <- "~/Projects/Spatial/NCAR/output/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))

sum(unique(USppt$station.id) %in% unique(UStemp$station.id))
##There are 4267 stations only for precipitation, 
##7651 stations for both, 
##and 474 stations only for tempreture,
##totally, there are 12,392 stations.

dim(USppt[!is.na(USppt[,8]),])
dim(UStemp[!is.na(UStemp[,8]),])
##There are 6,204,442 observations for precipitation
##and 4,285,841 observations for tempreture.
##Totally, there should be 11918*12*103=14,730,648 observations for precipitation
##and should be 8125*12*103=10,042,500 observations for tempreture.


#Count the observation number for each station.
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)
TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)



########################################################################
##QQ plot of observation count of each station for three types of response.
######################################################################
data <- c(Pcount, TMaxcount, TMincount)
names(data) <- NULL

data <- data.frame(response = data, type= rep(c("Precip","tmax","tmin"), c(length(Pcount), length(TMaxcount), length(TMincount))))

trellis.device(postscript, file = paste(outputdir, "Distn_of_Observation_count", ".ps", sep = ""), color=TRUE, paper="legal")
	a <- qqmath(~ log2(response) | type,	
		data = data,
		distribution = qunif,
		grid = TRUE,
		pch = 16,
		cex = 0.5,
        	xlab = list(label="f-value", cex=1.5),
	        ylab = list(label="Log of Number of Observation (log base 2 number)", cex=1.5),
	        scales = list(x = list(cex=1.5), y = list(cex=1.5)),
	)
	print(a)
dev.off()

trellis.device(postscript, file = paste(outputdir, "Rate_of_Observation_count", ".ps", sep = ""), color=TRUE, paper="legal")
    a <- qqmath(~ response/1236 | type,
        data = data,
        distribution = qunif,
        grid = TRUE,
        pch = 16,
        cex = 0.5,
            xlab = list(label="f-value", cex=1.5),
            ylab = list(label="Rate of Valid Observation", cex=1.5),
            scales = list(x = list(cex=1.5), y = list(cex=1.5)),
    )
    print(a)
dev.off()



################################################
##Line plot of observation count vs. time
###############################################
pcount <- ddply(.data = USppt,
		.variable = c("year","month"),
		.fun = summarize,
			count = sum(!is.na(precip))
)
pcount$group <- "precip"
tmincount <- ddply(.data = UStemp,
			.variable = c("year", "month"),
			.fun = summarize,
			   count = sum(!is.na(tmin))
)
tmincount$group <- "tmin"
tmaxcount <- ddply(.data = UStemp,
            .variable = c("year", "month"),
            .fun = summarize,
               count = sum(!is.na(tmax))
)
tmaxcount$group <- "tmax"
data <- do.call(rbind, list(pcount, tmaxcount, tmincount))
data$month <- as.numeric(data$month)
date <- paste(data$year, data$month, "01", sep="-")
data$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
data <- data[order(data$date),]
data$total <- c(11918, 8125, 8125)

start <- as.POSIXct(strptime(paste(head(data$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
end <- as.POSIXct(strptime(paste(tail(data$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col
trellis.device(postscript, file = paste(outputdir, "lineplot_of_observation_count", ".ps", sep = ""), color=TRUE, paper="legal")
b <- xyplot( log2(count) ~ date,
	data = data,
	type = c("l"),
	grib = TRUE,
	groups = group,
	xlab = list(label = "Time", cex = 1.5),
	xlim = c(start, end),
	ylab = list(label="Log of Number of Observation (log base 2 number)",cex=1.5),
	scales = list(y = list(relation = 'free', cex = 1.5), x = list(format = "%b %Y", tick.number = 14), cex = 1.2),                
    key = list(columns = 3, text = list(c("Precip", "Tmax", "Tmin"), cex = 1.2), lines=list(lwd = 3, col = col[1:3])),
	panel = function(x,y,...) {
		panel.abline(h=seq(0,13,by=1), v=seq(start, end, by="48 month"), color="lightgrey", lty=3, lwd=0.5)
		panel.xyplot(x,y,...)
	}
)
        print(b)
dev.off()

trellis.device(postscript, file = paste(outputdir, "na_rate_of_observation_count", ".ps", sep = ""), color=TRUE, paper="legal")
b <- xyplot( (1 - count/total) ~ date,
    data = data,
    type = c("l"),
    grib = TRUE,
    groups = group,
    xlab = list(label = "Time", cex = 1.5),
    xlim = c(start, end),
    ylab = list(label="Rate of Missing Value",cex=1.5),
    scales = list(y = list(relation = 'free', cex = 1.5), x = list(format = "%b %Y", tick.number = 14), cex = 1.2),
    key = list(columns = 3, text = list(c("Precip", "Tmax", "Tmin"), cex = 1.2), lines=list(lwd = 3, col = col[1:3])),
    panel = function(x,y,...) {
        panel.abline(h=seq(0,13,by=1), v=seq(start, end, by="48 month"), color="lightgrey", lty=3, lwd=0.5)
        panel.xyplot(x,y,...)
        panel.text(x[which.max(y)], max(y), round(max(y),2), adj=c(0,0))
		panel.text(x[which.min(y)], min(y), round(min(y),2), adj=c(1,1))
    }
)
        print(b)
dev.off()

#####################################
##measurement status on the US map
#####################################
us.map <- map('state', plot = FALSE, fill = TRUE)

data.tmax <- subset(UStemp, !is.na(tmax))[, !(names(UStemp) %in% c("station.id", "elev", "tmin", "station.name"))]
data.tmax <- subset(data.tmax, year %in% c("1895", "1969"))
data.tmin <- subset(UStemp, !is.na(tmin))[, !(names(UStemp) %in% c("station.id", "elev", "tmax", "station.name"))]
data.tmin <- subset(data.tmin, year %in% c("1895", "1969"))
data.precip <- subset(USppt, !is.na(precip))[, !(names(USppt) %in% c("station.id", "elev", "station.name"))]
data.precip <- subset(data.precip, year %in% c("1895", "1969"))

trellis.device(postscript, file = paste(outputdir, "tmax.status.condition.time", ".ps", sep = ""), color=TRUE, paper="legal")
	a <- xyplot(lat ~ lon | month*year,
		data = data.tmax,
        main = "Maximum Temperature Valid Observations",
		ylab = "Latitude",
		xlab = "Longitude",
		pch = 16,
		col = "red",
		cex = 0.2,
		layout = c(4,3),
		panel = function(x, y, ...){
			panel.polygon(us.map$x, us.map$y, lwd = 0.2)
			panel.xyplot(x, y, ...)
		}
	)
	print(a)
dev.off()
trellis.device(postscript, file = paste(outputdir, "tmin.status.condition.time", ".ps", sep = ""), color=TRUE, paper="legal")
    a <- xyplot(lat ~ lon | month*year,
        data = data.tmin,
        main = "Minimum Temperature Valid Observations",
        ylab = "Latitude",
        xlab = "Longitude",
        pch = 16,
        col = "red",
        cex = 0.2,
        layout = c(4,3),
        panel = function(x, y, ...){
            panel.polygon(us.map$x, us.map$y, lwd = 0.2)
            panel.xyplot(x, y, ...)
        }
    )
    print(a)
dev.off()
trellis.device(postscript, file = paste(outputdir, "precip.status.condition.time", ".ps", sep = ""), color=TRUE, paper="legal")
    a <- xyplot(lat ~ lon | month*year,
        data = data.precip,
		main = "Precipitation Valid Observations",
        ylab = "Latitude",
        xlab = "Longitude",
        pch = 16,
        col = "red",
        cex = 0.2,
        layout = c(4,3),
        panel = function(x, y, ...){
            panel.polygon(us.map$x, us.map$y, lwd = 0.2)
            panel.xyplot(x, y, ...)
        }
    )
    print(a)
dev.off()

