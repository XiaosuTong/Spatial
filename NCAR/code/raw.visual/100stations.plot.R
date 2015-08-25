#############################################
## scatter plot of raw obs for given stations
## 
## scatterRaw function
##   - input:
##       data: data.frame read from HDFS, which has 100 stations obs. 
##
##       outputdir: output path for plot
##
##       target: which variable to be plotted, tmax, tmin, or precip
##
##       size: the size of ps file, letter or legal
##
##       test: if it is TRUE, only this first station in stations vector will be plotted
##
###############################################
scatterRaw <- function(data=rst, outputdir=local.output, target="tmax", size = "letter", test = TRUE){

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

  stations <- unique(data$station.id)

  if(test) {
    stations <- stations[1]
  }

  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("scatter", "100stations", target, "ps", sep=".")), 
    color=TRUE, 
    paper=size
  )
    for(i in stations){   
      b <- xyplot( resp ~ time | factor
        , data = data
        , subset = station.id == i
        , xlab = list(label = "Month", cex = 1.5)
        , ylab = list(label = ylab, cex = 1.5)
        , main = list(label=paste("Station ", i, sep=""), cex=1)
        , type = "b"
        , pch = 16
        , cex = 0.5
        , aspect = "xy"
        , layout = c(1,9)
        , strip = FALSE
        , xlim = c(0, 143)
        , scales = list(
            y = list(alternating=TRUE, tick.number = 4), 
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
##       data: data.frame read from HDFS, which has 100 stations obs
## 
##       outputdir: output path for the plot
##
##       target: which variable to be plotted, tmax, tmin, or precip
##
##       size: the size of ps file, letter or legal
##
##       test: if it is TRUE, only this first station in stations vector will be plotted
##
#####################################################
monthRaw <- function(data=rst, outputdir=local.output, target = "tmax", size = "letter", test = TRUE) {

  library(plyr)
  
  dd <- ddply(
    .data = data,
    .vari = c("station.id","month"),
    .fun = summarise,
    mean = mean(resp)
  )

  data <- arrange(data, station.id, month)
  data$mean <- dd[rep(row.names(dd), each=103),]$mean

  stations <- unique(data$station.id)

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
    file = file.path(outputdir, paste("scatter", "100stations", target,"cond.month.ps", sep=".")), 
    color = TRUE, 
    paper = size
  )
    for(i in stations){
      b <- xyplot( resp-mean ~ as.numeric(year) | month
        , data = data
        , subset = station.id == i
        , xlab = list(label = "Year", cex = 1.5)
        , ylab = list(label = "Maximum Temperature (degrees centigrade)", cex = 1.5)
        , sub = list(label = paste("Station ", i, sep=""))
        , type = "p"
        , pch = 16
        , cex = 0.5      
        , layout = c(4,3)
        , strip = TRUE
        , scales = list(
            y = list(alternating=TRUE, cex = 1.2), 
            x = list(tick.number=5, cex = 1.2)
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


monthQQ <- function(data=rst, outputdir=local.output, target = "tmax", size = "letter", test = TRUE) {

  stations <- unique(data$station.id)

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
    file = file.path(outputdir, paste("QQ", "100stations", target,"cond.month.ps", sep=".")),
    color = TRUE, 
    paper=size
  ) 
    for(i in stations){ 
      a <- qqmath(~ resp | factor(month, levels = month.abb)
        , data = data
        , subset = station.id == i
        , distribution = qnorm
        , aspect = "xy"
        , layout = c(12,1)
        , pch = 16
        , cex = 0.5
        , scales = list(
            y = list(alternating=TRUE, cex = 1.2), 
            x = list(tick.number=3, cex = 1.2)
          )
        , sub = list(label= paste("Station ", i, sep=""))
        , xlab = list(label="Unit normal quantile", cex=1.5)
        , ylab = list(label=ylab, cex=1.5)
        , prepanel = prepanel.qqmathline
        , panel = function(x, y,...) {
            panel.grid()
            panel.qqmathline(x, y=x)
            panel.qqmath(x, y,...)
          }
      )
      print(a)
    } 
  dev.off()

}


monthSpatial <- function(data=rst, outputdir=local.output, target = "tmax", vari="lon", size = "letter") {

  if (target == "tmax") {
    ylab <- "Maximum Temperature (degrees centigrade)"
  } else if (target == "tmin") {
    ylab <- "Minimum Temperature (degrees centigrade)"
  } else {
    ylab <- "Precipitation (millimeters)"
  }

  if (vari == "lon") {
    xlab <- "Longitude"
  } else if (vari == "lat") {
    xlab <- "Latitude"
  } else {
    xlab <- "Log Base 2 Elevation"
  }

  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("scatter", "100stations", target, "vs", vari, "cond.month.ps", sep=".")),    
    color = TRUE, 
    paper = size
  )
    b <- xyplot( resp ~ get(vari) | factor(month, levels = month.abb)
      , data = rst
      , xlab = list(label = xlab, cex = 1.5)
      , ylab = list(label = ylab, cex = 1.5)
      #, strip = strip.custom(par.strip.text= list(cex = 1.5))
      #, par.settings = list(layout.heights = list(strip = 1.5))
      , scales = list(
          y = list(relation = 'free', cex=1.2, tick.number=4), 
          x = list(tick.number=5, cex=1.2, relation='same')
        )
      , prepanel = function(x,y,...) {
          if(vari == "elev") {
            prepanel.default.xyplot(log2(x), y, ...)
          } else {
            prepanel.default.xyplot(x, y, ...)
          }
        }
      , panel = function(x,y,...) {
          if(vari == "elev") {
            panel.xyplot(log2(x),y,pch=16, cex=0.2, col=col[1], ...)
            panel.loess(log2(x), y, degree=1, span=2/3, col = col[2], lwd = 2, ...)
          } else {
            panel.xyplot(x,y,pch=16, cex=0.2, col=col[1], ...)
            panel.loess(x, y, degree=1, span=2/3, col = col[2], lwd = 2, ...)
          }
          
        }
    )
    print(b)
  dev.off()

}