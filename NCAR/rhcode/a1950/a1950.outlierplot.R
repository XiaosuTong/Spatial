outlierVisAll <- function(dataset, byname=NULL, byvari=NULL) {

  tmp <- unique(dataset[, c("lon","lat")])
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersLoc.ps"),
    color = TRUE, 
    paper = "letter"
  )
  b <- xyplot(lat ~ lon
    , data = tmp
    , xlab = list(label="Logitude", cex=1.5)
    , ylab = list(label="Latitude", cex=1.5)
    , scale = list(cex=1.2)
    , panel = function(x,y,...) {
        panel.polygon(us.map$x,us.map$y)   
        panel.xyplot(x,y,...)
    }
  )
  print(b)
  dev.off()

  if (is.null(byname)) {
    foml <- abs(remainder-spafit) ~ elev2
  } else if (byname == "Winter") {
    foml <- abs(remainder-spafit) ~ elev2 | get(byvari)
  }
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("a1950.outliersElev", byname, ".ps", sep="")),
    color = TRUE, 
    paper = "letter"
  )
  b <- xyplot(foml
    , data = dataset
    , strip = strip.custom(factor.levels=c("Other months", "Winter months"))
    , xlab = list(label="Log Base 2 (Elevation + 128)", cex=1.5)
    , ylab = list(label="Absolute Value of Residual", cex=1.5)
    , scale = list(cex=1.2, y=list(at=c(seq(0,15,5), 2.5)))
    , panel = function(x,y,...) {
        panel.abline(h=c(seq(0,15,5),2.5), v=seq(6,12,1), col="lightgray", lwd=0.5)
        panel.xyplot(x,y,...)
        panel.loess(x,y,col="red", lwd=1, span=0.25, degree=1, evaluation = 200)
    }
  )
  print(b)
  dev.off()

  if (is.null(byname)) {
    keys <- NULL
  } else if (byname == "Winter") {
    keys <- list(
      text = list(label=c("Other months", "Winter months"), cex=1.2), 
      points = list(pch=1, cex=1, col=col[1:2]),
      columns=2
    )  
  }
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("a1950.outliersQuant", byname, ".ps", sep="")),
    color = TRUE, 
    paper = "letter"
  )
  if (is.null(byname)) {
    b <- qqmath(~(remainder-spafit)
      , data = dataset
      , key = keys
      , distribution = qunif
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Residual", cex=1.5)
      , scale = list(cex=1.2, y=list(at=seq(-15,15,5)))
      , panel = function(x,...) {
          panel.abline(h=seq(-15,15,5), v=seq(0,1,0.2), col="lightgray", lwd=1)
          panel.qqmath(x, pch = 16, cex = 0.5, ...)
      }
    )
  } else if (byname == "Winter") {
    b <- qqmath(~(remainder-spafit)
      , data = dataset
      , group = get(byvari)
      , key = keys
      , distribution = qunif
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Residual", cex=1.5)
      , scale = list(cex=1.2, y=list(at=seq(-15,15,5)))
      , panel = function(x,...) {
          panel.abline(h=seq(-15,15,5), v=seq(0,1,0.2), col="lightgray", lwd=1)
          panel.qqmath(x, pch = 16, cex = 0.5, ...)
      }
    )
  }
  print(b)
  dev.off()


  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersResidvsFit.ps"),
    color = TRUE, 
    paper = "letter"
  )
  b <- xyplot((remainder-spafit) ~ spafit
    , data = dataset
    , ylab = list(label="Residual", cex=1.5)
    , xlab = list(label="Spatial Fitted Value", cex=1.5)
    , xlim = c(-16,11)
    , scale = list(cex=1.2, y=list(at=seq(-15,15,5)), x=list(at=seq(-15,10,5)))
    , panel = function(x,y,...) {
        panel.abline(h=seq(-15,15,5), v=seq(-15,10,5), col="lightgray", lwd=1)
        panel.xyplot(x,y,...)
    }
  )
  print(b)
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersResidvsRemain.ps"),
    color = TRUE, 
    paper = "letter"
  )
  b <- xyplot((remainder-spafit) ~ remainder
    , data = outliers
    , ylab = list(label="Residual", cex=1.5)
    , xlab = list(label="Remainder", cex=1.5)
    , scale = list(cex=1.2)
    , panel = function(x,y,...) {
        panel.abline(h=seq(-15,15,5), v=seq(-15,10,5), col="lightgray", lwd=1)
        panel.xyplot(x,y,...)
    }
  )
  print(b)
  dev.off()

}

outlierVisStat <- function(dataset, byname=NULL, byvari=NULL) {

  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersCountbystation.ps"),
    color = TRUE, 
    paper = "letter"
  )
  b <- qqmath(~log2(out + 1)
    , data=dataset
    , distribution = qunif
    , xlab = list(label = "f-value", cex=1.5)
    , ylab = list(label = "Log Base 2 (Number of Outliers + 1)", cex=1.5)
    , ylim = c(-0.5, 8.5)
    , scale = list(cex=1.2, y = list(at=seq(0,8,2)))
    , panel = function(x,...) {
        panel.abline(h=seq(0,8,2), v=seq(0,1,0.2), col="lightgray", lwd=1)
        panel.qqmath(x, pch = 16, cex=0.5, ...)
    }
  )
  print(b)
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersCountbystation.ps"),
    color = TRUE, 
    paper = "letter"
  )
  b <- 

  )
  print(b)
  dev.off()  

}


outlierVisMonth <- function(dataset, byname=NULL, byvari=NULL) {

  if (is.null(byname)) {
    keys <- NULL
  } else if (byname == "Winter") {
    keys <- list(
      text = list(label=c("Other months", "Winter months"), cex=1.2), 
      points = list(pch=1, cex=1, col=col[1:2]),
      columns=2
    )  
  }
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("a1950.outliersCountmonth", byname, ".ps", sep="")),
    color = TRUE, 
    paper = "letter"
  )
  if (is.null(byname)) {
    b <- qqmath(~log2(out + 1)
      , data=dataset
      , key = keys
      , distribution = qunif
      , cex = 0.5
      , xlab = list(label = "f-value", cex=1.5)
      , ylab = list(label = "Log Base 2 (Number of Outliers + 1)", cex=1.5)
      , ylim = c(3.5, 8)
      , scale = list(cex=1.2, y = list(at=seq(3,7,1)))
      , panel = function(x,...) {
          panel.abline(h=seq(3,7,1), v=seq(0,1,0.2), col="lightgray", lwd=1)
          panel.qqmath(x,...)
      }
    )
  } else if (byname == "Winter") {
    b <- qqmath(~log2(out + 1)
      , data=dataset
      , key = keys
      , group = get(byvari)
      , distribution = qunif
      , cex = 0.5
      , xlab = list(label = "f-value", cex=1.5)
      , ylab = list(label = "Log Base 2 (Number of Outliers + 1)", cex=1.5)
      , ylim = c(3.5, 8)
      , scale = list(cex=1.2, y = list(at=seq(3,7,1)))
      , panel = function(x,...) {
          panel.abline(h=seq(3,7,1), v=seq(0,1,0.2), col="lightgray", lwd=1)
          panel.qqmath(x,...)
      }
    )
  }
  print(b)
  dev.off()
  
  tmp <- ddply(.data=dataset, .vari="month", .fun=summarise, out=sum(out))
  tmp$month <- factor(tmp$month, levels=arrange(tmp, out)$month)
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersCountCondmonth.ps"),
    color = TRUE, 
    paper = "letter"
  )
  b <- dotplot( month ~ out
    , data = tmp
    , cex = 1
    , xlab = list(label = "Number of Outliers", cex=1.5)
    , ylab = list(label = "Month", cex=1.5)
    , scale = list(cex=1.2)
    , panel = function(x,y,...) {
        panel.dotplot(x,y,...)
    }
  )
  print(b)
  dev.off()

  date <- paste(dataset$year, match(dataset$month, month.abb), "01", sep="-")
  dataset$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
  dataset <- arrange(dataset, year, match(month, month.abb))
  start <- head(dataset$date, 1)
  end <- tail(dataset$date, 1)
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersCountvsTime.ps"),
    color = TRUE, 
    paper = "letter"
  )
  b <- xyplot( out ~ date
    , data = dataset
    , type = "l"
    , ylim = c(-5, 205)
    , ylab = list(label = "Number of Outliers", cex=1.5)
    , xlab = list(label = "Month", cex=1.5)
    , scale = list(cex=1.2,
        y = list(at = seq(0, 200, 50)),
        x = list(format = "%b %Y", at = c(seq(start, end, by="96 month"), end))
      )
    , panel = function(x,y,...) {
        panel.abline(v=c(seq(start, end, by="96 month"), end), col = "lightgray", lwd=0.5)
        panel.abline(h=seq(0, 200, 50), col = "lightgray", lwd=0.5)
        panel.xyplot(x,y,...)
        panel.loess(x,y,span=0.2, degree=1, evaluation=200, col = "red", ... )
    }
  )
  print(b)
  dev.off()

}

plotEng.rawOutlier <- function(data, station, lim=2) {

  data <- arrange(data, date)
  data$flag <- with(data, remainder - spafit <= -lim | remainder - spafit >= lim)
  data$factor <- factor(
    x = rep(paste("Period", 1:6), c(rep(108,5), 36)),
    levels = paste("Period", c(6:1))
  )
  data$time <- c(rep(0:107, times = 5), 0:35) 

  b <- xyplot( resp ~ time | factor
    , data = data
    , xlab = list(label = "Month", cex = 1.5)
    , ylab = list(label = ylab, cex = 1.5)
    , sub = list(label=paste("Station", station), cex=1.2)
    , layout = c(1,6)
    , strip = FALSE,
    , xlim = c(0, 107)
    , ylim = c(min(c(data$resp, data$fitted), na.rm=TRUE)-1, max(c(data$resp, data$fitted), na.rm=TRUE)+1)
    , key=list(
        cex = 1.2,
        text = list(label=c("raw", "outlier", "fitted")), 
        lines = list(pch=16, cex=0.7, lwd=1.5, type=c("p","p","l"), col=c(col[1], "red", col[2])),
        columns=3
      )
    , scales = list(
        y = list(tick.number=4, cex=1.2), 
        x = list(at=seq(0, 107, by=12), relation='same', cex=1.2)
      )
    , panel = function(x,y, subscripts,...) {
        panel.abline(v=seq(0,108, by=12), color="lightgrey", lty=3, lwd=0.5)
        panel.superpose(x, y, subscripts, groups = data$flag, panel.group = "panel.xyplot", type="p", pch=16, cex=c(0.5,0.6), col=c(col[1],"red"),...)
        if (!any(grepl("fc", names(data)))) {
          panel.xyplot(data[subscripts,]$time, (data[subscripts,]$trend+data[subscripts,]$seasonal), type="l", col=col[2], lwd=1, ...)            
        } else {
          panel.xyplot(data[subscripts,]$time, (data[subscripts,]$data.seasonal+data[subscripts,]$fc.first+data[subscripts,]$fc.second), type="l", col=col[2], lwd=1, ...)
        }
      }
  )
  return(b)

}

plotEng.outlierMonth <- function(dataset, ym, lim=2) {

#  hexbinplot(rnorm(10000)~rnorm(10000), colramp= function(n){heat.ob (n, beg=1, end=256)})

  colcold <- BTC(4, beg=50, end=140)
  colhot <- magent(4, beg=140, end=50)
  col <- c(colcold, colhot)
  dataset <- subset(dataset, remainder-spafit <= -lim | remainder-spafit >= lim)
  b <- xyplot( lat ~ lon
    , data = dataset
    , xlab = list(label="Longitude", cex = 1.5)
    , ylab = list(label="Latitude", cex = 1.5)
    , sub = paste(ym[1], ym[2])
    , xlim = c(-125,-66)
    , ylim = c(24, 50)
    , key = list(
        type = "p", 
        text = list(label=c("<-15","[-15, -10)","[-10, -5)", "[-5, -2)", "(2, 5]", "(5, 10]", "(10, 15]",">15")),  
        points = list(cex=1, pch=16, col=col[1:8]), 
        columns = 4
      )
    , pch = 16
    , cex = 1
    , scales = list(
        y = list(cex = 1.2),x = list(cex = 1.2)
      )
    , panel = function(x,y,...) {
        panel.polygon(us.map$x,us.map$y) 
        panel.text(x=-123, y=26, labels=paste("Outliers:", nrow(dataset)), adj=c(0,0))
        for(i in 1:4) {
          sub <- subset(dataset, remainder-spafit > (0 + (i-1)*5) & remainder-spafit <= ( 5 + (i-1)*5))
          panel.xyplot(sub$lon, sub$lat, col=col[i+4], ...)
          sub <- subset(dataset, remainder-spafit < (0 - (i-1)*5) & remainder-spafit >= (-5 - (i-1)*5))
          panel.xyplot(sub$lon, sub$lat, col=col[5-i], ...)
        }
      }
  )
  return(b)

}

outlierMonth <- function(input, output, name, plotEng=plotEng.outlierMonth) {
 
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      mindex <- (as.numeric(map.keys[[r]][1]) - 1950)*12 + match(map.keys[[r]][2], month.abb)
      if (mindex %in% outliers.a1950.month) {
        rhcounter("map", "month", 1)
        value <- subset(map.values[[r]], flag == 1, select = c(spafit, remainder, lon, lat, elev2, station.id))
        pp <- plotEng(dataset=value, ym=map.keys[[r]])
        rhcollect(mindex, serialize(pp, NULL))
      }
    })
  }) 
  job$setup <- expression(
    map = {
      load("outliersMonthTop.a1950.RData")
      library(maps)
      library(lattice)
      library(hexbin, lib.loc=lib.loc)
      us.map <- map('state', plot = FALSE, fill = TRUE)
    }
  )
  job$shared <- c(
    file.path(file.path(rh.root, par$dataset, "a1950", "Rdata", "outliersMonthTop.a1950.RData"))
  )
  job$parameters <- list(
    plotEng = plotEng
  )
  job$mapred <- list(
    mapred.reduce.tasks = 1,
    mapred.tasktimeout = 0
  )
  job$input <- rhfmt(input, type="sequence")
  job$output <- rhfmt(output, type="sequence")
  job$mon.sec <- 10
  job$jobname <- output
  job$readback <- TRUE
  job.mr <- do.call("rhwatch", job) 

  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("a1950", name, "ps", sep=".")),
    color = TRUE, 
    paper = "letter"
  )
  for(ii in order(unlist(lapply(job.mr, "[[", 1)))) {
    print(unserialize(job.mr[[ii]][[2]]))
  }
  dev.off()

}

