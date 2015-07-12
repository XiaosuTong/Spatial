##############################################################################
##Precipitation and temperature for the 100 most observations stations
##Raw data AGAINST time (month index).
##Difference between this code and old code file is change the period to be 9.
##############################################################################

#Set up the directory and load the data.
library(lattice)
library(plyr)

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

datadir <- "~/Projects/Spatial/NCAR/RData"
outputdir <- "~/Projects/Spatial/NCAR/output"
load(file.path(datadir, "USmonthlyMet.RData"))
load(file.path(datadir, "stations.100.RData"))

#############################################
## scatter plot of raw obs for given stations
## 
## scatterRaw function
##   - input:
##       stations: station.id vector
##       
##       data: UStemp or USppt
##
##       target: which variable to be plotted, tmax, tmin, or precip
##
##       size: the size of ps file, letter or legal
##
##       test: if it is TRUE, only this first station in stations vector will be plotted
##
###############################################
scatterRaw <- function(stations=stations.tmax, data=UStemp, target = "tmax", size = "letter", test = TRUE){
  
  data <- subset(data, station.id %in% stations)
  
  data  <- arrange(data, station.id, year, month)

  data$factor <- factor(
    x = rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times=100), 
    levels = paste("Period", c(9:1))
  )
  data$time <- c(rep(0:143,8), 0:83) 
  
  if (target == "tmax") {
    ylab <- "Maximum Temperature (degrees centigrade)"
  } else if (target == "tmin") {
    ylab <- "Minimum Temperature (degrees centigrade)"
  } else {
    ylab <- "Precipitation (millimeters)"
  }

  if(test) {
    stations <- stations[1]
  }

  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("scatterplot", target,"100stations.ps", sep=".")), 
    color=TRUE, 
    paper=size
  )
    for(i in stations){   
      b <- xyplot( get(target) ~ time | factor
        , data = data
        , subset = station.id == i
        , xlab = list(label = "Month", cex = 1.2)
        , ylab = list(label = ylab, cex = 1.2)
        , main = list(label=paste("Station ", i, sep=""), cex=1)
        , type = "b"
        , pch = 16
        , cex = 0.5
        , aspect = "xy"
        , layout = c(1,9)
        , strip = FALSE
        , xlim = c(0, 143)
        , scales = list(
            y = list(relation = 'same', alternating=TRUE), 
            x = list(at=seq(0,143,by=12), relation='same')
          )
        , panel = function(x, y, ...) {
            panel.abline( v=seq(0,143,by=12), color="lightgrey", lty=3, lwd=0.5)
            panel.xyplot(x, y, ...)
          }
      )
      print(b)
    }
  dev.off()

}


######################################################
## scatter plot of raw obs for given stations conditional on month
## 
## monthRaw function
##   - input:
##       stations: station.id vector
##       
##       data: UStemp or USppt
##
##       target: which variable to be plotted, tmax, tmin, or precip
##
##       size: the size of ps file, letter or legal
##
##       test: if it is TRUE, only this first station in stations vector will be plotted
##
#####################################################
monthRaw <- function(stations=stations.tmax, data=UStemp, target = "tmax", size = "letter", test = TRUE) {

  library(plyr)

  data <- subset(data, station.id %in% stations)

  data <- arrange(data, station.id, month, year)
  
  dd <- ddply(
    .data = data,
    .vari = c("station.id","month"),
    .fun = summarise,
    mean = mean(get(target))
  )
  data$mean <- dd[rep(row.names(dd), each=103),]$mean

  if (target == "tmax") {
    ylab <- "Maximum Temperature (degrees centigrade)"
  } else if (target == "tmin") {
    ylab <- "Minimum Temperature (degrees centigrade)"
  } else {
    ylab <- "Precipitation (millimeters)"
  }

  if(test) {
    stations <- stations[1]
  }

trellis.device(
  device = postscript, 
  file = file.path(outputdir, paste("scatter", target,"100stations.cond.month.ps", sep=".")), 
  color = TRUE, 
  paper = size
)
  for(i in stations){
    b <- xyplot( get(target)-mean ~ year | month
      , data = data
      , subset = station.id == i
      , xlab = list(label = "Year", cex = 1.2)
      , ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.2)
      , main = list(label = paste("Station ", i, sep=""))
      , type = "p"
      , pch = 16
      , cex = 0.5      
      , layout = c(4,3)
      , strip = TRUE
      , scales = list(
          y = list(relation = 'same', alternating=TRUE), 
          x = list(tick.number=10, relation='same')
        )
      , key=list(
        text = list(label=c("observation","loess smoothing")), 
        lines = list(pch=16, cex=0.7, lwd=1.5, type=c("p","l"), col=col[1:2]),
        columns=2
      )
      , panel = function(x,y,...) {
          panel.abline(h=seq(-10,10,by=5), v=seq(1900,2000,by=20), color="lightgrey", lty=3, lwd=0.5)
          panel.xyplot(x,y,...)
          panel.loess(x,y,span=2/3,degree=1,family="symmetric",evaluation = 100,col=col[2],...)
        }
    )
    print(b)
  }
dev.off()  

}