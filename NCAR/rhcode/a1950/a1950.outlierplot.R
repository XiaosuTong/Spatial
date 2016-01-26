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

}


xyplot( remainder-spafit ~ lon | equal.count(lat, 20, overlap=0)
  , data = outliers
)

xyplot( lat ~ lon
  , data = outliers
  , subset = station.id %in% subset(outlierBystation, out>=2^4)$station.id
  , panel = function(x,y,...) {
      panel.polygon(us.map$x,us.map$y)   
      panel.xyplot(x,y,...)
  }
)


