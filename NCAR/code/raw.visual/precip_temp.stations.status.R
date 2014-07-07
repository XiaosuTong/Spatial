#############################################################
##Measurement status plot for each station to see the where are missing
#############################################################

#Set up the directory and load the data.
library(lattice)
library(plyr)
datadir <- "~/Projects/Spatial/NCAR/RData/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))
outputdir <- "~/Projects/Spatial/NCAR/output/"

###################
##Precipitation
###################
USppt$month <- factor(USppt$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))

#Count the observation number for each station.
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)

#Order the stations by the observation counts
Pcount <- sort(Pcount, decreasing=TRUE)
tmp.precip <- USppt
tmp.precip$station.id <- factor(tmp.precip$station.id, levels = names(Pcount))

#month <- tmp.precip$month
#levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
#month <- as.numeric(factor(month, levels=c(1:12)))
#date <- paste(tmp.precip$year, month, "01", sep="-")
#tmp.precip$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
tmp.precip <- tmp.precip[order(tmp.precip$station.id, tmp.precip$year, tmp.precip$month),]
tmp.precip$time <- rep(0:1235,11918)
tmp.precip$obs <- as.numeric(!is.na(tmp.precip$precip))

#each page has 50 stations, totally 239 pages, 50*238+18=11918
tmp.precip$page <- c(rep(1:238,each=1236*50), rep(239, 1236*18))

#indx is the variable that assigned the "y" value to each station. 
#Two station is different by 2 in "y" scale which is the blank between two stations.
tmp.precip$indx <- c(rep(rep(seq(1, 100, by=2), each=1236),238), rep(seq(1,36,by=2), each=1236))
#rst is the final "y" value for each observation. It is consist of the "y" value of the station and 0.6*obs.
tmp.precip$rst <- tmp.precip$obs*0.6 + tmp.precip$indx


#start <- as.POSIXct(strptime(paste(head(tmp.precip$year,1), "01", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")
#end <- as.POSIXct(strptime(paste("2000", "12", "01", sep="-"), format = "%Y-%m-%d"), format='%Y%m%d', tz="")

trellis.device(postscript, file = paste(outputdir, "precipitation_measurement_status_against_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in 1:239){
	data <- subset(tmp.precip, page==i)
	if(sum(data$obs)==1236*50){
        b <- xyplot( rst ~ time ,
             data = data,
             xlab = list(label = "Month", cex = 1.5),
             ylab = list(label = "Measurement status", cex = 1.5),
             type = "p",
             col = c("red","black"),
             distribute.type= TRUE,
	     groups = as.factor(obs),
             pch=16,
             cex=0.2,
#             scales = list(y = list(relation = 'same', cex=1, labels=c("0","1"), alternating=TRUE), x=list(relation= 'same',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list( y =list(draw = FALSE), x=list(tick.number = 10, relation='same', cex=1.2)),
             panel = function(x,y,...) {
                  panel.abline( v= seq(0,1235, by=60), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
		  }
        )
        print(b)
        }
	else{
        b <- xyplot( rst ~ time ,
             data = data,
             xlab = list(label = "Month", cex = 1.5),
             ylab = list(label = "Measurement status", cex = 1.5),
             type = "p",
             col = c("black","red"),
             distribute.type= TRUE,
             groups = as.factor(obs),
             pch=16,
             cex=0.2,
#             scales = list(y = list(relation = 'same', cex=1, labels=c("0","1"), alternating=TRUE), x=list(relation= 'same',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list( y =list(draw = FALSE), x=list(tick.number = 10, relation='same', cex=1.2)),
             panel = function(x,y,...) {
                  panel.abline( v= seq(0,1235, by=60), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  }
        )
        print(b)
        }
     }
dev.off()

##########################
##Maximum Temperature
##########################
UStemp$month <- factor(UStemp$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))

TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)
TMaxcount <- sort(TMaxcount, decreasing=TRUE)
tmp.tmax <- UStemp
tmp.tmax$station.id <- factor(tmp.tmax$station.id, levels = names(TMaxcount))

tmp.tmax <- tmp.tmax[order(tmp.tmax$station.id, tmp.tmax$year, tmp.tmax$month),]
tmp.tmax$time <- rep(0:1235,8125)
tmp.tmax$obs <- as.numeric(!is.na(tmp.tmax$tmax))

#each page has 50 stations, totally 163 pages, 50*162+25=8125
tmp.tmax$page <- c(rep(1:162,each=1236*50), rep(163, 1236*25))

#indx is the variable that assigned the "y" value to each station. 
#Two station is different by 2 in "y" scale which is the blank between two stations.
tmp.tmax$indx <- c(rep(rep(seq(1,100,by=2), each=1236),162), rep(seq(1,50,by=2), each=1236))
#rst is the final "y" value for each observation. It is consist of the "y" value of the station and 0.6*obs.
tmp.tmax$rst <- tmp.tmax$obs*0.6 + tmp.tmax$indx

trellis.device(postscript, file = paste(outputdir, "tmax_measurement_status_against_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in 1:163){
        data <- subset(tmp.tmax, page==i)
        if(sum(data$obs)==1236*50){
        b <- xyplot( rst ~ time ,
             data = data,
             xlab = list(label = "Month", cex = 1.5),
             ylab = list(label = "Measurement status", cex = 1.5),
             type = "p",
             col = c("red","black"),
             distribute.type= TRUE,
             groups = as.factor(obs),
             pch=16,
             cex=0.2,
#             scales = list(y = list(relation = 'same', cex=1, labels=c("0","1"), alternating=TRUE), x=list(relation= 'same',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list( y =list(draw = FALSE), x=list(tick.number = 10, relation='same', cex=1.2)),
             panel = function(x,y,...) {
                  panel.abline( v= seq(0,1235, by=60), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  }
        )
        print(b)
        }
        else{
        b <- xyplot( rst ~ time ,
             data = data,
             xlab = list(label = "Month", cex = 1.5),
             ylab = list(label = "Measurement status", cex = 1.5),
             type = "p",
             col = c("black","red"),
             distribute.type= TRUE,
             groups = as.factor(obs),
             pch=16,
             cex=0.2,
#             scales = list(y = list(relation = 'same', cex=1, labels=c("0","1"), alternating=TRUE), x=list(relation= 'same',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list( y =list(draw = FALSE), x=list(tick.number = 10, relation='same', cex=1.2)),
             panel = function(x,y,...) {
                  panel.abline( v= seq(0,1235, by=60), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  }
        )
        print(b)
        }
     }
dev.off()


####################################
##Minimum Temperature
##################################
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)
TMincount <- sort(TMincount, decreasing=TRUE)
tmp.tmin <- UStemp
tmp.tmin$station.id <- factor(tmp.tmin$station.id, levels = names(TMincount))

tmp.tmin <- tmp.tmin[order(tmp.tmin$station.id, tmp.tmin$year, tmp.tmin$month),]
tmp.tmin$time <- rep(0:1235,8125)
tmp.tmin$obs <- as.numeric(!is.na(tmp.tmin$tmin))

#each page has 50 stations, totally 163 pages, 50*162+25=8125
tmp.tmin$page <- c(rep(1:162,each=1236*50), rep(163, 1236*25))

#indx is the variable that assigned the "y" value to each station. 
#Two station is different by 2 in "y" scale which is the blank between two stations.
tmp.tmin$indx <- c(rep(rep(seq(1,100,by=2), each=1236),162), rep(seq(1,50,by=2), each=1236))
#rst is the final "y" value for each observation. It is consist of the "y" value of the station and 0.6*obs.
tmp.tmin$rst <- tmp.tmin$obs*0.6 + tmp.tmin$indx

trellis.device(postscript, file = paste(outputdir, "tmin_measurement_status_against_month",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in 1:163){
        data <- subset(tmp.tmin, page==i)
        if(sum(data$obs)==1236*50){
        b <- xyplot( rst ~ time ,
             data = data,
             xlab = list(label = "Month", cex = 1.5),
             ylab = list(label = "Measurement status", cex = 1.5),
             type = "p",
             col = c("red","black"),
             distribute.type= TRUE,
             groups = as.factor(obs),
             pch=".",
             cex=0.2,
#             scales = list(y = list(relation = 'same', cex=1, labels=c("0","1"), alternating=TRUE), x=list(relation= 'same',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list( y =list(draw = FALSE), x=list(tick.number = 10, relation='same', cex=1.2)),
             panel = function(x,y,...) {
                  panel.abline( v= seq(0,1235, by=60), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  }
        )
        print(b)
        }
        else{
        b <- xyplot( rst ~ time ,
             data = data,
             xlab = list(label = "Month", cex = 1.5),
             ylab = list(label = "Measurement status", cex = 1.5),
             type = "p",
             col = c("black","red"),
             distribute.type= TRUE,
             groups = as.factor(obs),
             pch=".",
             cex=0.2,
#             scales = list(y = list(relation = 'same', cex=1, labels=c("0","1"), alternating=TRUE), x=list(relation= 'same',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list( y =list(draw = FALSE), x=list(tick.number = 10, relation='same', cex=1.2)),
             panel = function(x,y,...) {
                  panel.abline( v= seq(0,1235, by=60), color="lightgrey", lty=3, lwd=0.5)
                  panel.xyplot(x,y,...)
                  }
        )
        print(b)
        }
     }
dev.off()

