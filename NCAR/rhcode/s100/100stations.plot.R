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


#########################################################
##  Diagnostic plots for components from stl2 fitting  ## 
#########################################################
#########################################
##  fitted value and obs against time  ##
#########################################
fitRaw <- function(data=rst, outputdir, target="tmax", size = "letter", test = TRUE){

  data$factor <- factor(
    x = rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times=100),
    levels = paste("Period", c(9:1))
  )
  data$time <- c(rep(0:143,8), 0:83) 

  stations <- unique(data$station.id)

  if(test) {
    stations <- stations[1]
  }

  if (target == "tmax") {
    ylab <- "Maximum Temperature (degrees centigrade)"
  } else if (target == "tmin") {
    ylab <- "Minimum Temperature (degrees centigrade)"
  } else {
    ylab <- "Precipitation (millimeters)"
  }

  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("fitted.time", "100stations", target, "ps", sep=".")),
    color = TRUE, 
    paper = size
  )
    for(i in stations){
      sub <- subset(data, station.id==i)
      b <- xyplot( resp ~ time | factor,
        , data = sub
        , xlab = list(label = "Month", cex = 1.5)
        , ylab = list(label = ylab, cex = 1.5)
        , main = list(label=paste("Station ", i, sep=""), cex=1)
        , layout = c(1,9)
        , aspect= "xy"
        , strip = FALSE,
        , xlim = c(0, 143)
        , key=list(
          text = list(label=c("raw","fitted")), 
          lines = list(pch=16, cex=0.7, lwd=1.5, type=c("p","l"), col=col[1:2]),
          columns=2
        )
        , scales = list(
            y = list(tick.number=4), 
            x = list(at=seq(0, 143, by=12), relation='same')
          )
        , panel = function(x,y,subscripts,...) {
            panel.abline(v=seq(0,145, by=12), color="lightgrey", lty=3, lwd=0.5)
            panel.xyplot(x, y, type="p", col=col[1], pch=16, cex=0.5, ...)
            if (!any(grepl("fc", names(rst)))) {
              panel.xyplot(sub[subscripts,]$time, (sub[subscripts,]$trend+sub[subscripts,]$seasonal), type="l", col=col[2], lwd=1, ...)            
            } else {
              panel.xyplot(sub[subscripts,]$time, (sub[subscripts,]$data.seasonal+sub[subscripts,]$fc.first+sub[subscripts,]$fc.second), type="l", col=col[2], lwd=1, ...)
            }
          }
      )
      print(b)
    }
  dev.off()

}

################################################################
##  diagnostic plots for remainders, QQ, scatter, line plots  ## 
################################################################
remainderDiag <- function(data=rst, outputdir, target="tmax", size = "letter", test=TRUE) {

  stations <- unique(data$station.id)

  if (target == "tmax") {
    ylab <- "Maximum Temperature (degrees centigrade)"
  } else if (target == "tmin") {
    ylab <- "Minimum Temperature (degrees centigrade)"
  } else {
    ylab <- "Precipitation (millimeters)"
  }
  
  idx <- grep("remainder", names(data))
  
  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("QQ.remainder", "100stations", target, "ps", sep=".")), 
    color = TRUE, 
    paper = size
  )
    a <- qqmath( ~ data[, idx] | unique(station.id),
      , data = data
      , distribution = qnorm
      , aspect = 1
      , layout = c(5,3)
      , xlab = list(label="Unit normal quantile", cex=1.5)
      , ylab = list(label=ylab, cex=1.5)
      , scales = list(x = list(cex=1.2), y = list(cex=1.2))
      , prepanel = prepanel.qqmathline
      , panel = function(x, y,...) {
          panel.abline(h= seq(-15,15,by=5), v=c(-2,0,2), lty=1, lwd=1.5, col="lightgrey")
          panel.qqmathline(x, y=x)
          panel.qqmath(x, y, pch=16, cex=0.3, ...)
        }
    )
    print(a)
  dev.off()

  if(test) {
    stations <- stations[1]
  }
  
  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("remainder.time", "100stations", target, "ps", sep=".")),  
    color = TRUE, 
    paper=size
  )
  for(i in stations){
    b <- xyplot( data[, idx] ~ date,
      , data = data
      , subset = station.id == i
      , main = list(label=paste("Station ", i, sep=""), cex=1)
      , xlab = list(label = "Month", cex = 1.5)
      , ylab = list(label = ylab, cex = 1.5)
      , xlim = c(0, 1235)
      , key = list(
          text=list(label=c("remainder","degree=2,span=0.15","degree=1,span=0.35")), 
          lines=list(pch=16, cex=1, lwd=2, type=c("p","l", "l"), col=col[1:3]), 
          columns=3
        )
      , scales = list(
          y = list(relation = 'same', alternating=TRUE, cex = 1.2), 
          x = list(at=seq(0, 1235, by=120), cex = 1.2, relation='same')
        )
      , panel = function(x,y,...) {
          panel.abline(h=0)
          panel.xyplot(x, y, pch=16, cex = 0.6, ...)
          panel.loess(x,y,degree=2,span=0.15, col=col[2], lwd = 2, evaluation=200,...)
          panel.loess(x,y,degree=1,span=0.35, col=col[3], lwd = 2, evaluation=200,...)
        }
    )
    print(b)
  }
  dev.off()

  data$factor <- factor(
    x = rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times=100),
    levels = paste("Period", c(3,2,1,6,5,4,9,8,7))
  )
  data$time <- c(rep(0:143,8), 0:83)

  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("remainder.time2", "100stations", target, "ps", sep=".")), 
    color = TRUE, 
    paper = size
  )
    for(i in stations){
      b <- xyplot( data[, idx] ~ time | factor,
        , data = data
        , subset = station.id == i
        , xlab = list(label = "Month", cex = 1.5),
        , ylab = list(label = ylab, cex = 1.5),
        , main = list(label=paste("Station ", i, sep=""), cex=1)
        , layout = c(1,3),
        , xlim = c(0, 143),
        , strip = FALSE,
        , scales = list(
            y = list(at = seq(-10,10,5), cex=1.2), 
            x = list(at=seq(0, 143, by=12), relation='same', cex=1.2)
          )
        , panel = function(x,y,...) {
            panel.abline(h=0, v=seq(0,143, by=12), color="lightgrey", lty=3, lwd=0.5)
            panel.xyplot(x,y, pch=16, cex=0.6, ...)
            panel.loess(x,y, degree=2,span=1/4, col=col[2], lwd=2, evaluation=200, ...)
          }
      )
      print(b)
   }
  dev.off()

  ACF <- ddply(
    .data = data,
    .variables = "station.id",
    .fun = function(r) {
      corr <- acf(r[, idx], plot=FALSE)
      data.frame(
        correlation = corr$acf,
        lag = corr$lag 
      )
    }
  )

  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("remainder.acf", "100stations", target, "ps", sep=".")), 
    color = TRUE, 
    paper = size
  )
    b <- xyplot( correlation ~ lag | station.id
      , data = ACF
      , subset = lag != 0
      , layout = c(2,2)
      , xlab = list(label = "Lag", cex = 1.5)
      , ylab = list(label = "ACF", cex = 1.5)
      , scales = list(
          y = list(cex=1.2), 
          x = list(relation='same', cex=1.2)
        )
      , panel = function(x,y,...) {
          panel.abline(h=0)
          panel.xyplot(x,y, type="h", lwd=1.5,...)
        }
    )
    print(b)
  dev.off()

 trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("remainder.month", "100stations", target, "ps", sep=".")),
    color = TRUE, 
    paper=size
  )
    for(i in stations){
      b <- xyplot( data[,idx] ~ as.numeric(year) | factor(month, levels = month.abb)
        , data = data
        , subset = station.id == i
        , xlab = list(label = "Year", cex = 1.5)
        , ylab = list(label = ylab, cex = 1.5)
        , main = list(label=paste("Station ", i, sep=""), cex=1)
        , layout = c(12,1)
        , scales = list(
           y = list(cex=1.2), 
           x = list(tick.number=2, cex=1.2)
          )
        , key = list(
            text=list(label=c("remainder","loess smoothing")), 
            lines=list(pch=16, cex=1, lwd=2, type=c("p","l"), col=col[1:2]), 
            columns=2
          )
        , panel = function(x,y,...){
            panel.abline(h=0, color="black", lty=1)
            panel.xyplot(x,y, cex = 0.6, pch=16, ...)
            panel.loess(x,y,span=3/4, degree=1, lwd = 2, col=col[2],...)
          }
      )
      print(b)
    }
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("remainder.month2", "100stations", target, "ps", sep=".")),
    color = TRUE, 
    paper = size
  )
    for(i in stations){
      b <- xyplot( data[,idx] ~ as.numeric(year) | factor(month, levels = month.abb)
        , data = data
        , subset = station.id==i
        , xlab = list(label = "Year", cex = 1.5)
        , ylab = list(label = ylab, cex = 1.5)
        , layout = c(2,6)
        , main = list(label=paste("Station ", i, sep=""), cex=1)
        , scales = list(
            y = list(cex=1.2, at=seq(-10,10,5)), 
            x = list(tick.number=10, relation='same', cex=1.2)
          )
        , panel = function(x,y,...) {
            panel.abline(h=0, color="black", lty=1)
            panel.xyplot(x, y, type="b", pch=16, cex=0.5, ...)
          }
      )
      print(b)
    }
  dev.off()


  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("QQ.remainder.month", "100stations", target, "ps", sep=".")),
    color = TRUE, 
    paper = size
  )
    for(i in stations){
      a <- qqmath( ~ data[, idx] | factor(month, levels = month.abb),
        , data = data
        , subset = station.id==i
        , distribution = qnorm
        , aspect = 1
        , layout = c(4,3),
        , scales = list(
            y = list(cex=1.2, at=seq(-10,10,5)), 
            x = list(cex=1.2)
          )
        , main = list(label=paste("Station ", i, sep=""), cex=1)
        , xlab = list(label="Unit normal quantile", cex=1.5)
        , ylab = list(label = ylab, cex = 1.5)
        , prepanel = prepanel.qqmathline,
        , panel = function(x, y,...) {
            panel.abline(h= seq(-15,15,by=5), v=c(-2,0,2), lty=1, lwd=1.5, col="lightgrey")
            panel.qqmathline(x, y=x)
            panel.qqmath(x, y, pch=16, cex=0.4,...)
          }
      )
      print(a)
    }
  dev.off()

}


############################################################
##  diagnostic plots for remainders, scatter, line plots  ##
############################################################
seasonalDiag <- function(data=rst, outputdir, target="tmax", size = "letter", test=TRUE) {

  rorder <- ddply(
    .data = rst,
    .variables = "station.id",
    .fun = summarise,
    mean = mean(resp)
  )
  order.st <- arrange(rorder, mean)$station.id

  idx <- grep("seasonal", names(data))

  if(test) {
    order.st <- order.st[1]
  }
  
  if (target == "tmax") {
    ylab <- "Maximum Temperature (degrees centigrade)"
  } else if (target == "tmin") {
    ylab <- "Minimum Temperature (degrees centigrade)"
  } else {
    ylab <- "Precipitation (millimeters)"
  }

  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("seasonal.month", "100stations", target, "ps", sep=".")),
    color = TRUE, 
    paper = size
  )
    for(i in order.st){
      b <- xyplot( data[, idx] ~ as.numeric(year) | factor(month, levels = month.abb)
        , data = data
        , subset = station.id == i
        , main = list(label=paste("Station ", i, sep=""), cex=1)
        , xlab = list(label = "Year", cex = 1.5)
        , ylab = list(label = ylab, cex = 1.5)
        , layout = c(12,1)
        , scales = list(
            y = list(cex=1.2), 
            x = list(tick.number=2, relation='same', cex=1.2)
          )
        , panel = function(x,y,...) {
            panel.abline(h=0, color="black", lty=1)
            panel.xyplot(x,y,cex=0.4, pch=16, ...)
          }
      )
      print(b)
    }
  dev.off()

  mmean <- ddply(
    .data = data,
    .variable = c("station.id","month"),
    .fun = summarise,
    mean = mean(resp)
  )
  mmean <- arrange(mmean, match(month, month.abb))

  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("seasonalperiod.month", "100stations", target, "ps", sep=".")),
    color = TRUE, 
    paper = size
  )
    b <- xyplot( mean ~ match(month, month.abb) | station.id,
      , data = mmean
      , xlab = list(label = "Month", cex = 1.5)
      , ylab = list(label = ylab, cex = 1.5)
      , layout = c(5,4)
      , scales = list(
          y = list(cex=1.2), 
          x = list(at=c(1, 3, 5, 7, 9, 11), relation='same', cex=1.2)
        )
      , panel = function(x,y,...) {
          panel.xyplot(x,y,, type="b", cex=0.5, pch=16, ...)
        }
    )
    print(b)
  dev.off()
  
  smean <- ddply(
    .data = rst,
    .variables = c("station.id","month"),
    .fun = function(r) {
      data.frame(mean = mean(r[, idx]))
    }
  )
  searmd <- merge(x=data, y=smean, by=c("station.id","month"), all.x=TRUE)
  ridx <- grep("remainder", names(data))
  
  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("sea+remaind.month", "100stations", target, "ps", sep=".")), 
    color = TRUE, 
    paper = size
  )
    for(i in order.st){
      sub <- arrange(subset(searmd, station.id == i), date)
      b <- xyplot( (sub[, idx] - sub$mean + sub[, ridx]) ~ as.numeric(year) | factor(month, levels = month.abb),
        , data = sub
        , xlab = list(label = "Year", cex = 1.5)
        , ylab = list(label = ylab, cex = 1.5)
        , main = list(label=paste("Station ", i, sep=""), cex=1)
        , key=list(
            text = list(label=c("seasonal","seasonal+remainder")), 
            lines = list(pch=16, cex=1, lwd=2, type=c("l","p"), col=col[2:1]), 
            columns = 2
          )
        , layout = c(12,1)
        , scales = list(
            y = list(relation = 'same', alternating=TRUE, cex=1.2), 
            x = list(tick.number=2, relation='same', cex=1.2)
          )
        , panel = function(x,y,subscripts,...) {
            panel.abline(h=0, color="black", lty=1)
            panel.xyplot(x, y, col=col[1], pch=16, cex=0.6, ...)
            panel.xyplot(as.numeric(sub$year)[subscripts], (sub[subscripts,idx]-sub$mean[subscripts]), type="l", lwd=2, col=col[2], ...)
          }
      )
      print(b)
    }
  dev.off()

}


#######################################################
##  diagnostic plots for trend, scatter, line plots  ##
#######################################################
trendDiag <- function(data=rst, outputdir, target="tmax", size = "letter", test=TRUE, fc = FALSE) {
  
  if (target == "tmax") {
    ylab <- "Maximum Temperature (degrees centigrade)"
  } else if (target == "tmin") {
    ylab <- "Minimum Temperature (degrees centigrade)"
  } else {
    ylab <- "Precipitation (millimeters)"
  }
  
  if (fc) {
    vari <- "fc.first"
    klab <- "low frequency component"
  } else {
    klab <- "trend component"
    vari <- "trend"
  }
  
  tmean <- ddply(
    .data = data,
    .variables = c("station.id","year"),
    .fun = function(r) {
      data.frame(mean = mean(r$resp))
    }
  )
  tmv <- ddply(
    .data = arrange(tmean, year),
    .variables = "station.id",
    .fun = summarise,
    mv = as.numeric(filter(mean, rep(1/10,10), sides=2)),
    year = year
  )
  trendmean <- merge(x=data, y=tmv, by=c("station.id","year"), all.x=TRUE)

  trellis.device(
    device = postscript, 
    file = file.path(outputdir, paste("trend", "100stations", target, "ps", sep=".")),
    color = TRUE, 
    paper = size
  )
    b <- xyplot( get(vari) ~ date | station.id, 
      , data = trendmean
      , xlab = list(label = "Month", cex = 1.5)
      , ylab = list(label = ylab, cex = 1.5)
      , layout = c(4,3)
      , key=list(
          text = list(label=c(klab,"moving average of yearly mean")), 
          lines = list(pch=16, cex=1, lwd=2, type=c("l","p"), col=col[1:2]),
          columns=2
        )
      , prepanel = function(x,y,subscripts,...){ 
          v <- trendmean[subscripts,] 
          ylim <- range(v$mv, na.rm=TRUE) 
          ans <- prepanel.default.xyplot(v$date, v[, vari], ...) 
          ans$ylim <- range(c(ans$ylim, ylim), na.rm=TRUE) 
          ans 
        } 
      , scales = list(
          y = list(relation = 'free', cex=1.2), 
          x = list(at=c(0,500,1000), relation = 'same', cex=1.2)
        )
      , panel = function(x, y, subscripts, ...) {
          v <- trendmean[subscripts,]
          panel.xyplot(v$date[seq(1,1236,12)], v$mv[seq(1,1236,12)], type="p", pch=16, cex=0.5, col=col[2], ...)
          panel.xyplot(x, y, type="l", lwd=2, col=col[1], ...)
        }
    )
    print(b)
  dev.off()

  if (fc) {

    trendmean <- merge(x=data, y=tmean, by=c("station.id","year"), all.x=TRUE)

    tmv <- ddply(
      .data = arrange(trendmean, year),
      .variables = "station.id",
      .fun = summarise,
      mv = as.numeric(filter(mean-fc.first, rep(1/10,10), sides=2)),
      year = year
    )
    trendmean <- cbind(arrange(data, station.id, year), mv=tmv[, "mv"])

    trellis.device(
      device = postscript, 
      file = file.path(outputdir, paste("fc", "100stations", target, "ps", sep=".")),
      color = TRUE, 
      paper = size
    )
      b <- xyplot( fc.second ~ date | station.id,
        , data = trendmean,
        , xlab = list(label = "Month", cex = 1.5),
        , ylab = list(label = ylab, cex = 1.5),
        , xlim = c(0, 1236),
        , layout = c(2,2),
        , between = list(x=0.5)
        , key=list(
            text = list(label=c(
                "middle frequency component",
                "moving average of (yearly mean-low frequency component)"
            )),
            lines = list(pch=16, cex=1, lwd=2, type=c("l","p"), col=col[1:2]), 
            columns = 2
          )
        , scales = list(
            y = list(relation = 'same', cex=1.2), 
            x = list(at=seq(0,1236,480), relation = 'same', cex=1.2)
          )
        , prepanel = function(x,y,subscripts,...){
            v <- trendmean[subscripts,]
            ylim <- range(v$mv, na.rm=TRUE) 
            ans <- prepanel.default.xyplot(v$date, v$fc.second, ...) 
            ans$ylim <- range(c(ans$ylim, ylim), na.rm=TRUE) 
            ans 
          }
        , panel = function(x, y, subscripts, ...) {
            panel.abline(v=seq(0,1236, by=480), h=0, color="lightgrey", lty=1, lwd=0.5) 
            v <- trendmean[subscripts, ]
            panel.xyplot(v$date[seq(6,1236,12)], v$mv[seq(6,1236,12)], type="p", pch=16, cex=0.5, col=col[2], ...)
            panel.xyplot(x, y, type="l", lwd=2, col=col[1], ...)        
          }
      )
      print(b)
    dev.off()

  }

}


