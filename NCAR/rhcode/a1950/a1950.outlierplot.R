outlierVisAll <- function(dataset, by=NULL) {

  tmp <- unique(outliers[, c("lon","lat")])
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersLoc.ps"),
    color = TRUE, 
    paper = "letter"
  )
  xyplot(lat ~ lon
    , data = tmp
    , xlab = list(label="Logitude", cex=1.5)
    , ylab = list(label="Latitude", cex=1.5)
    , scale = list(cex=1.2)
    , panel = function(x,y,...) {
        panel.polygon(us.map$x,us.map$y)   
        panel.xyplot(x,y,...)
    }
  )
  dev.off()
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersElevWinter.ps"),
    color = TRUE, 
    paper = "letter"
  )
  xyplot(abs(remainder-spafit) ~ elev2 | mFlag
    , data = outliers
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
  dev.off()
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersQuant.ps"),
    color = TRUE, 
    paper = "letter"
  )
  qqmath(~(remainder-spafit)
    , data = outliers
    , distribution = qunif
    , xlab = list(label="f-value", cex=1.5)
    , ylab = list(label="Residual", cex=1.5)
    , scale = list(cex=1.2, y=list(at=seq(-15,15,5)))
    , panel = function(x,...) {
        panel.abline(h=seq(-15,15,5), v=seq(0,1,0.2), col="lightgray", lwd=1)
        panel.qqmath(x,...)
    }
  )
  dev.off()
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersQuantWinter.ps"),
    color = TRUE, 
    paper = "letter"
  )
  qqmath(~(remainder-spafit)
    , data = outliers
    , group = mFlag
    , key=list(
        text = list(label=c("Other months", "Winter months"), cex=1.2), 
        points = list(pch=1, cex=1, col=col[1:2]),
        columns=2
      )  
    , distribution = qunif
    , xlab = list(label="f-value", cex=1.5)
    , ylab = list(label="Residual", cex=1.5)
    , scale = list(cex=1.2, y=list(at=seq(-15,15,5)))
    , panel = function(x,...) {
        panel.abline(h=seq(-15,15,5), v=seq(0,1,0.2), col="lightgray", lwd=1)
        panel.qqmath(x,...)
    }
  )
  dev.off()
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersResidvsFit.ps"),
    color = TRUE, 
    paper = "letter"
  )
  xyplot((remainder-spafit) ~ spafit
    , data = outliers
    , ylab = list(label="Residual", cex=1.5)
    , xlab = list(label="Spatial Fitted Value", cex=1.5)
    , xlim = c(-16,11)
    , scale = list(cex=1.2, y=list(at=seq(-15,15,5)), x=list(at=seq(-15,10,5)))
    , panel = function(x,y,...) {
        panel.abline(h=seq(-15,15,5), v=seq(-15,10,5), col="lightgray", lwd=1)
        panel.xyplot(x,y,...)
    }
  )
  dev.off()

}

outlierVisStat <- function(dataset, by=NULL) {

  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersCountbystation.ps"),
    color = TRUE, 
    paper = "letter"
  )
  qqmath(~log2(out + 1)
    , data=outlierBystation
    , distribution = qunif
    , cex = 0.5
    , xlab = list(label = "f-value", cex=1.5)
    , ylab = list(label = "Log Base 2 (Number of Outliers + 1)", cex=1.5)
    , ylim = c(-0.5, 8.5)
    , scale = list(cex=1.2, y = list(at=seq(0,8,2)))
    , panel = function(x,...) {
        panel.abline(h=seq(0,8,2), v=seq(0,1,0.2), col="lightgray", lwd=1)
        panel.qqmath(x,...)
    }
  )
  dev.off()

}


outlierVisMonth <- function(dataset, by=NULL) {

  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersCountmonthWinter.ps"),
    color = TRUE, 
    paper = "letter"
  )
  qqmath(~log2(out + 1)
    , data=outlierBymonth
    , group = mFlag
    , key=list(
        text = list(label=c("Other months", "Winter months"), cex=1.2), 
        points = list(pch=1, cex=1, col=col[1:2]),
        columns=2
      )
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
  dev.off()
  tmp <- ddply(.data=outlierBymonth, .vari="month", .fun=summarise, out=sum(out))
  tmp$month <- factor(tmp$month, levels=arrange(tmp, out)$month)
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.outliersCountCondmonth.ps"),
    color = TRUE, 
    paper = "letter"
  )
  dotplot( month ~ out
    , data = tmp
    , cex = 1
    , xlab = list(label = "Number of Outliers", cex=1.5)
    , ylab = list(label = "Month", cex=1.5)
    , scale = list(cex=1.2)
    , panel = function(x,y,...) {
        panel.dotplot(x,y,...)
    }
  )
  dev.off()

}