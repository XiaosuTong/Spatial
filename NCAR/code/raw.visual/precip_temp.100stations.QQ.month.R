###################################################
##Precipitation and temperature, 100 stations
##QQ plot of raw data conditional on month
###################################################

#load the data
library(lattice)
library(maps)
datadir <- "~/Projects/Spatial/NCAR/RData/"
outputdir <- "~/Projects/Spatial/NCAR/output/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))

#Count the observation number for each station.
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)
TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)

#Attend to add a small map in the corner of the QQ plot
#us.map <- map('state', plot = FALSE, fill = TRUE)


#find the temperature data for the 100 stations
#remove all NA observations
TMaxcount <- sort(TMaxcount, decreasing=TRUE)
stations <- head(names(TMaxcount),100)
UStemp$month <- factor(UStemp$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))

tmp <- UStemp[UStemp[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$tmax)),]

#Create the QQ plot of temperature conditional on month for one station
trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_tmax_of_month", ".ps", sep = ""), color=TRUE, paper="legal")	
    for(i in stations){	
	a <- qqmath(~ tmax | month,
		data = tmp[tmp[,1]==i,],
        distribution = qnorm,
		aspect = "xy",
		layout = c(12,1),
		pch = 16,
		cex = 0.5, 
		main = list(label= paste("Station ", i, sep=""), cex=1.5),
        xlab = list(label="Unit normal quantile", cex=1.2),
        ylab = list(label="Max Temperature(degrees centigrade)", cex=1.2),
#       scales = list(x = list(cex=1.5), y = list(cex=1.5)),
		prepanel = prepanel.qqmathline,
  		panel = function(x, y,...) {
    			panel.grid()
    			panel.qqmathline(x, y=x)
    			panel.qqmath(x, y,...)
  		}

	)
	print(a)
#	b <- xyplot(lat ~ lon,
#                data = tmp[tmp[,1]==i,],
#                xlab = NULL,
#                ylab = NULL,
##                layout = c(1,1),
#                pch = 16,
#                cex = 0.3,
#		xlim = c(-125,-65),
#		ylim = c(23,50),
#                col = "red",
##		scales = list(x = list(draw=FALSE), y = list(draw=FALSE)),
##                strip = strip.custom(par.strip.text= list(cex = 1.5)),
##                par.settings = list(layout.heights = list(strip = 1.5)),
#                panel = function(x,y,...) {
#                        panel.polygon(us.map$x,us.map$y,lwd=0.2)
#                        panel.xyplot(x,y,...)
#                }
#        )
#	plot.new()
#	title(paste("Station ", i, sep=""), cex=1.5)
#	print(a, pos=c(0,0,1,0.65), newpage=FALSE, more=TRUE)
#	print(b, pos=c(0,0.6,0.4,1), more=FALSE)
    }
dev.off()

#find the temperature data for the 100 stations
#remove all NA observations
TMincount <- sort(TMincount, decreasing=TRUE)
stations <- head(names(TMincount),100)
	
tmp <- UStemp[UStemp[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$tmin)),]

#Create the QQ plot of temperature conditional on month for one station
trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_tmin_of_month", ".ps", sep = ""), color=TRUE, paper="legal")                 
    for(i in stations){
        a <- qqmath(~ tmin | month,
                data = tmp[tmp[,1]==i,],
                distribution = qnorm,
		aspect = "xy",
		pch = 16,
		cex = 0.5,
		layout = c(12,1),
                main = list(label= paste("Station ", i, sep=""), cex=1.5),
                xlab = list(label="Unit normal quantile", cex=1.2),
                ylab = list(label="Min Temperature(degrees centigrade)", cex=1.2),
#                scales = list(x = list(cex=1.5), y = list(cex=1.5)),
                prepanel = prepanel.qqmathline,
                panel = function(x, y,...) {
                        panel.grid()
                        panel.qqmathline(x, y=x)
                        panel.qqmath(x, y,...)
                }

        )
	print(a)
#        b <- xyplot(lat ~ lon,
#                data = tmp[tmp[,1]==i,],
#                xlab = NULL,
#                ylab = NULL,
##                layout = c(1,1),
#                pch = 16,
#                cex = 0.3,
#                xlim = c(-125,-65),
#                ylim = c(23,50),
#                col = "red",
##               scales = list(x = list(draw=FALSE), y = list(draw=FALSE)),
##                strip = strip.custom(par.strip.text= list(cex = 1.5)),
##                par.settings = list(layout.heights = list(strip = 1.5)),
#                panel = function(x,y,...) {
#                        panel.polygon(us.map$x,us.map$y,lwd=0.2)
#                        panel.xyplot(x,y,...)
#                }
#        )
#        plot.new()
#        title(paste("Station ", i, sep=""), cex=1.5)
#        print(a, pos=c(0,0,1,0.65), newpage=FALSE, more=TRUE)
#        print(b, pos=c(0,0.6,0.4,1), more=FALSE)
    }
dev.off()


###################
##Precipitation
###################
Pcount <- sort(Pcount, decreasing=TRUE)
stations <- head(names(Pcount),100)
USppt$month <- factor(USppt$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))

tmp <- USppt[USppt[,1] %in% stations,]
tmp <- tmp[!(is.na(tmp$precip)),]

#Create the QQ plot of temperature conditional on month for one station
trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_precipitation_of_month", ".ps", sep = ""), color=TRUE, paper="legal")
    for(i in stations){
        a <- qqmath(~ precip | month,
                data = tmp[tmp[,1]==i,],
                distribution = qnorm,
                aspect = "xy",
		pch = 16,
		cex = 0.5,
                layout = c(12,1),
                main = list(label= paste("Station ", i, sep=""), cex=1.5),
                xlab = list(label="Unit normal quantile", cex=1.2),
            	ylab = list(label = "Precipitation (millimeters)", cex = 1.2),
#                scales = list(x = list(cex=1.5), y = list(cex=1.5)),
                prepanel = prepanel.qqmathline,
                panel = function(x, y,...) {
                        panel.grid()
                        panel.qqmathline(x, y=x)
                        panel.qqmath(x, y,...)
                }

        )
	print(a)
#        b <- xyplot(lat ~ lon,
#                data = tmp[tmp[,1]==i,],
#                xlab = NULL,
#                ylab = NULL,
##                layout = c(1,1),
#                pch = 16,
#                cex = 0.3,
#                xlim = c(-125,-65),
#                ylim = c(23,50),
#                col = "red",
##               scales = list(x = list(draw=FALSE), y = list(draw=FALSE)),
##                strip = strip.custom(par.strip.text= list(cex = 1.5)),
##                par.settings = list(layout.heights = list(strip = 1.5)),
#                panel = function(x,y,...) {
#                        panel.polygon(us.map$x,us.map$y,lwd=0.2)
#                        panel.xyplot(x,y,...)
#                }
#        )
#        plot.new()
#        title(paste("Station ", i, sep=""), cex=1.5)
#        print(a, pos=c(0,0,1,0.65), newpage=FALSE, more=TRUE)
#        print(b, pos=c(0,0.6,0.4,1), more=FALSE)
    }
dev.off()

