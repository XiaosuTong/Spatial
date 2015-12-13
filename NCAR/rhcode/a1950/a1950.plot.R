intpolat.visual <- function(size = "letter", surf, SPsize, check=NULL) {
  
  rst1 <- rhread(file.path(rh.root, par$dataset, "a1950", "bymonth.fit", "symmetric", surf, "1", "MSE"))[[1]][[2]]
  rst2 <- rhread(file.path(rh.root, par$dataset, "a1950", "bymonth.fit", "symmetric", surf, "2", "MSE"))[[1]][[2]]
  rst <- rbind(rst1, rst2)
  rst$degree <- rep(c(1,2), each = nrow(rst1))

  if (is.null(check)) {

    if(SPsize == "median") {
      sub <- subset(rst, span %in% c(0.008, 0.01, 0.015, 0.02, 0.04))
    } else if (SPsize == "large") {
      sub <- subset(rst, span %in% seq(0.06, 0.1, 0.01))
    } else if (SPsize == "small"){
      sub <- subset(rst, span %in% seq(0.005, 0.009, 0.001))
    }

    trellis.device(
      device = postscript, 
      file = file.path(local.root, "output", paste("QuanMSE", "a1950", par$dataset, "span", SPsize, "ps", sep=".")), 
      color=TRUE, 
      paper=size
    )
    b <- qqmath(~mse|as.factor(span)
      , data = sub
      , subset = mse <= 4
      , group = degree
      , dist = qunif
      , cex = 0.5
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Mean Square Error", cex=1.5)
      , key=list(
          text = list(label=c("degree=1","degree=2")),
          lines = list(pch=1, cex=1, type="p", col=col[1:2]), 
          columns = 2
        )
      , layout = c(length(unique(sub$span)), 1)
      , scale = list(cex=1.2)
      , panel = function(x,...) {
          panel.abline(h=seq(0,4,0.5), v=seq(0,1,0.25), col="lightgray")
          panel.qqmath(x,...)
        }
    )
    print(b)
    dev.off()
    
    trellis.device(
      device = postscript, 
      file = file.path(local.root, "output", paste("QuanMSE", "a1950", par$dataset, "degree", SPsize ,"ps", sep=".")), 
      color=TRUE, 
      paper=size
    )
    b <- qqmath(~mse|as.factor(degree)
      , data = sub
      , group = span
      , subset = mse <= 4
      , dist = qunif
      , cex = 0.5
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Mean Square Error", cex=1.5)
      , key=list(
          text = list(label=paste("span=", sort(unique(sub$span)), sep="")),
          lines = list(pch=1, cex=1, type="p", col=col[1:length(unique(sub$span))]), 
          columns = length(unique(sub$span))
        )
      , layout = c(2, 1)
      , scale = list(cex=1.2)
      , panel = function(x,...) {
          panel.abline(h=seq(0,4,0.5), v=seq(0,1,0.2), col="lightgray")
          panel.qqmath(x,...)
        }  

    )
    print(b)
    dev.off()

  } else {

    sub <- subset(rst, (span == check[[1]][1] & degree == check[[1]][2]) | (span == check[[2]][1] & degree == check[[2]][2]))
    trellis.device(
      device = postscript, 
      file = file.path(local.root, "output", paste("QuanMSE", "a1950", par$dataset, "span", "check", "ps", sep=".")), 
      color=TRUE, 
      paper=size
    )
    b <- qqmath(~mse
      , data = sub
      , group = degree
      , dist = qunif
      , cex = 0.5
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Mean Square Error", cex=1.5)
      , key=list(
          text = list(label=c("degree=1, span=0.008","degree=2, span=0.015")),
          lines = list(pch=1, cex=1, type="p", col=col[1:2]), 
          columns = length(unique(sub$span))
        )
      , scale = list(cex=1.2)
      , panel = function(x,...) {
          panel.abline(h=seq(0,4,0.5), v=seq(0,1,0.2), col="lightgray")
          panel.qqmath(x,...)
        }
    )
    print(b)
    dev.off()

  } 
  
}


intpolat.visualNew <- function(size = "letter", surf="direct") {

  rst1 <- rhread(file.path(rh.root, par$dataset, "a1950", "bymonth.fit.new", "symmetric", surf, "1", "MABSE"))[[1]][[2]]
  rst2 <- rhread(file.path(rh.root, par$dataset, "a1950", "bymonth.fit.new", "symmetric", surf, "2", "MABSE"))[[1]][[2]]
  rst <- rbind(rst1, rst2)
  rst$degree <- rep(c(1,2), each = nrow(rst1))

    trellis.device(
      device = postscript, 
      file = file.path(local.root, "output", paste("QuanMABSE", "a1950", par$dataset, "span", "ps", sep=".")), 
      color=TRUE, 
      paper=size
    )
    b <- qqmath(~mse|as.factor(span)
      , data = rst
      , group = degree
      , dist = qunif
      , cex = 0.5
      , layout = c(4,1)
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Mean Square Error", cex=1.5)
      , key=list(
          text = list(label=c("degree=1","degree=2")),
          lines = list(pch=1, cex=1, type="p", col=col[1:2]), 
          columns = 2
        )
      , scale = list(cex=1.2)
      , panel = function(x,...) {
          panel.abline(h=seq(0,2,0.2), v=seq(0,1,0.25), col="lightgray")
          panel.qqmath(x,...)
        }
    )
    print(b)
    dev.off()

  for(i in c("small","median","large")) {
    if(i == "small") {
      sub <- subset(rst, as.numeric(span) <= 0.035)
    }else if (i == "median") {
      sub <- subset(rst, as.numeric(span)<=0.085 & as.numeric(span) > 0.035)
    }else {
      sub <- subset(rst, as.numeric(span)>0.085)
    }

    trellis.device(
      device = postscript, 
      file = file.path(local.root, "output", paste("QuanMABSE", "a1950", par$dataset, i, "degree","ps", sep=".")), 
      color=TRUE, 
      paper=size
    )
    b <- qqmath(~mse|as.factor(degree)
      , data = sub
      , group = span
      , dist = qunif
      , cex = 0.5
      , layout = c(2,1)
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Mean Square Error", cex=1.5)
      , key=list(
          text = list(label=paste("span=", sort(unique(sub$span)), sep="")),
          lines = list(pch=1, cex=1, type="p", col=col[1:length(unique(sub$span))]), 
          columns = length(unique(sub$span))
        )
      , scale = list(cex=1.2)
      , panel = function(x,...) {
          panel.abline(h=seq(0,2,0.2), v=seq(0,1,0.25), col="lightgray")
          panel.qqmath(x,...)
        }
    )
    print(b)
    dev.off()
  }
} 












#########################################################
##  Diagnostic plots for components from stl2 fitting  ## 
#########################################################
#########################################
##  fitted value and obs against time  ##
#########################################
a1950.fitRaw <- function(data=rst, outputdir, target="tmax", size = "letter", St.num = 128, test = TRUE){

  data$factor <- factor(
    x = rep( rep(paste("Period", 1:6), c(rep(108,5), 36)), times=St.num),
    levels = paste("Period", c(6:1))
  )
  data$time <- c(rep(0:107, times = 5), 0:35) 

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
    file = file.path(outputdir, paste("fitted.time", "a1950", target, "ps", sep=".")),
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
        , layout = c(1,6)
        , strip = FALSE,
        , xlim = c(0, 107)
		, ylim = c(min(c(sub$resp, sub$fitted), na.rm=TRUE), max(c(sub$resp, sub$fitted), na.rm=TRUE))
        , key=list(
          text = list(label=c("raw", "interpolate", "fitted")), 
          lines = list(pch=16, cex=0.7, lwd=1.5, type=c("p","p","l"), col=col[c(1,3,2)]),
          columns=3
        )
        , scales = list(
            y = list(tick.number=4), 
            x = list(at=seq(0, 107, by=12), relation='same')
          )
        , panel = function(x,y,subscripts,...) {
            panel.abline(v=seq(0,108, by=12), color="lightgrey", lty=3, lwd=0.5)
			fit <- subset(sub[subscripts,], flag == 0)
			obs <- subset(sub[subscripts,], flag == 1)
			panel.xyplot(obs$time, obs$resp, type="p", col=col[1], pch=16, cex=0.5, ...)
            panel.xyplot(fit$time, fit$fitted, type="p", col=col[3], pch=16, cex=0.5, ...)
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
