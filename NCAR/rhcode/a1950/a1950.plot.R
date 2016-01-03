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


imputeCrossValid <- function(surf="direct", Edeg = TRUE) {

  if (Edeg) {
    layout <- c(2,1)
    fo <- ~ mse | factor(degree)
    rst1 <- rhread(file.path(rh.root, par$dataset, "a1950", "bymonth.fit.cv", "symmetric", surf, "1", "MABSE"))[[1]][[2]]
    rst2 <- rhread(file.path(rh.root, par$dataset, "a1950", "bymonth.fit.cv", "symmetric", surf, "2", "MABSE"))[[1]][[2]]
    rst <- rbind(rst1, rst2)
    rst$degree <- rep(c(1,2), each = nrow(rst1))
    sub <- subset(rst, span %in% c(0.003, 0.005, 0.015, 0.035,0.095) & mse <=10)

    trellis.device(
      device = postscript, 
      file = file.path(local.root, "output", paste("QuanMABSE", "a1950", par$dataset, "span", "ps", sep=".")), 
      color=TRUE, 
      paper="letter"
    )
    b <- qqmath(~mse|as.factor(span)
      , data = sub
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
      , scale = list(cex=1.2, y=list(relation="free"))
      , panel = function(x,...) {
          if(max(x)<5) {
            panel.abline(h=seq(0,3,0.5),v=seq(0,1,0.2), col="lightgray")
          } else {
            panel.abline(h=seq(0,10,2),v=seq(0,1,0.2), col="lightgray")
          }
          panel.qqmath(x,...)
        }
    )
    print(b)
    dev.off()

  } else {
    rst <- rhread(file.path(rh.root, par$dataset, "a1950", "bymonth.fit.cv", "symmetric", surf, "0", "MABSE"))[[1]][[2]]
    layout <- c(1,1)
    fo <- ~ mse
  }

    sub <- subset(rst, span %in% c(0.003, 0.005, 0.015, 0.035,0.095) & mse <=10)

    trellis.device(
      device = postscript, 
      file = file.path(local.root, "output", paste("QuanMABSE", "a1950", par$dataset, "degree","ps", sep=".")), 
      color=TRUE, 
      paper="letter"
    )
    b <- qqmath( fo
      , data = sub
      , group = span
      , dist = qunif
      , cex = 0.5
      , layout = layout
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Mean Square Error", cex=1.5)
      , key=list(
          text = list(label=paste("span=", sort(unique(sub$span)), sep="")),
          lines = list(pch=1, cex=1, type="p", col=col[1:length(unique(sub$span))]), 
          columns = length(unique(sub$span)),
          cex = 1.2
        )
      , scale = list(cex=1.2, y=list(relation="free"))
      , panel = function(x,...) {
          panel.abline(h=seq(0,10,2), v=seq(0,1,0.2), col="lightgray")
          panel.qqmath(x,...)
        }
    )
    print(b)
    dev.off()

} 


###################################################################
##  Visualization plots for the spatial loess fit for missing    ##
##  value imputing.                                              ##
##  input files are /a1950/bymonth.fit or /a1950/bymonth.fit.new ##
###################################################################
Qsample <- function(x, num) {

  xx <- x[!is.na(x)]
  xx <- sort(xx)
  print(xx)
  len <- length(xx)
  idx <- round(seq(1, len, length.out = num))
  f.value <- (idx - 0.5) / len
  value <- data.frame(
    resid = xx[idx],
    fv = f.value
  )
  value

}

a1950.spaImputeVisual <- function(family = "symmetric", surf = "direct", Edeg = 2, span = 0.05) {
  
  FileInput <- file.path(
    rh.root, par$dataset, "a1950", "bymonth.fit", 
    family, surf, Edeg, paste("sp", span, sep="")
  )

  atlevels <- c(seq(-4,-1,0.2),seq(-0.8,0.8,0.1), seq(1,4,0.2))

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      v <- map.values[[r]]
      value <- Qsample(v$resp - v$fitted, 1000)
      value$month <- map.keys[[r]][2]
      value$year <- as.numeric(map.keys[[r]][1])
      rhcollect(1, value)
    })
  })
  job$reduce <- expression(
    pre = {
      combine <- data.frame()
    },
    reduce = {
      combine <- rbind(combine, do.call("rbind", reduce.values))
    },
    post = {
      rhcollect(reduce.key, combine)
    }
  )  
  job$parameters <- list(
    Qsample = Qsample
  )
  job$input <- rhfmt(FileInput, type = "sequence")
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bymonth.fit.plot", family, surf, Edeg, paste("sp",span, sep="")), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 1,  #cdh3,4
    mapreduce.job.reduces = 1  #cdh5
  )
  job$readback <- TRUE
  job$combiner <- TRUE
  job$jobname <- file.path(rh.root, par$dataset, "a1950", "bymonth.fit.plot", family, surf, Edeg, paste("sp",span, sep=""))
  job.mr <- do.call("rhwatch", job)

  ## normal quantile plots which includes all quantiles inside and outside of [0.015, 0.985]
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("a1950", "spaloessResid", "bytime", "ps", sep=".")), 
    color=TRUE, 
    paper="letter"
  )
  for(i in c(1950:1997)) {
    sub <- subset(job.mr[[1]][[2]], year == i)
    b <- xyplot(resid ~ qnorm(fv) | factor(month, levels=month.abb)
      , data = sub
      , layout = c(6, 2)
      , sub = paste("Year", i)
      , ylab = list(label="Spatial Residuals", cex=1.5)
      , xlab = list(label="Normal Quantiles", cex=1.5)
      , scale = list(cex=1.2)
      , key=list(
          text = list(label=c("inside [0.015, 0.985]", "outside [0.015, 0.985]"), cex=1.2), 
          lines = list(pch=16, cex=1, type=c("p"), col=col[c(1,2)]),
          columns = 2
        )
      , panel = function(x, y, ...) {
          idx <- which(pnorm(y) <= 0.985 & pnorm(y) >= 0.015)
          idx2 <- which(pnorm(y) > 0.985 | pnorm(y) < 0.015)
          panel.abline(h=0, col="black", lwd=0.5, lty=1)
          panel.xyplot(x[idx], y[idx], col = col[1], pch = 16, cex =0.5, ...)
          panel.xyplot(x[idx2], y[idx2], col = col[2], pch = 16, cex =0.5, ...)
      }
    )
    print(b)
  }
  dev.off()

  ## Normal quantiles plot which only includes quantiles between 0.015 and 0.985
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("a1950", "spaResidcenter", "bytime", "ps", sep=".")), 
    color=TRUE, 
    paper="letter"
  )
  for(i in c(1950:1997)) {
    sub <- subset(job.mr[[1]][[2]], year == i & fv <= 0.985 & fv >= 0.015)
    b <- xyplot(resid ~ qnorm(fv) | factor(month, levels=month.abb)
      , data = sub
      , layout = c(4, 3)
      , sub = paste("Year", i)
      , ylab = list(label="Spatial Residuals", cex=1.5)
      , xlab = list(label="Normal Quantiles", cex=1.5)
      , scale = list(cex=1.2)
      , aspect = 1
      , panel = function(x, y, ...) {
          panel.qqmathline(y,y=y,...)
          panel.xyplot(x, y, col = col[1], pch = 16, cex =0.3, ...)
      }
    )
    print(b)
  }
  dev.off()

  ## The third plot is the contourplot of the smoothed residuals over the US map
  new.grid <- expand.grid(
    lon = seq(-124, -67, by = 0.25),
    lat = seq(25, 49, by = 0.25)
  )
  instate <- !is.na(map.where("state", new.grid$lon, new.grid$lat))
  new.grid <- new.grid[instate, ]
  if(Edeg != 0) {
    load(file.path(local.root, "RData", "stations.a1950.RData"))
    load(file.path(local.root, "RData", "info.RData"))
    loc <- subset(UStinfo, station.id %in% stations.a1950.tmax)[,-5]
    elev.fit <- spaloess( elev ~ lon + lat,
      data = loc,
      degree = 2, 
      span = 0.05,
      distance = "Latlong",
      normalize = FALSE,
      napred = FALSE
    )
    grid.fit <- predloess(
      object = elev.fit,
      newdata = data.frame(
        lon = new.grid$lon,
        lat = new.grid$lat
      )
    )
    new.grid$elev2 <- log2(grid.fit + 128)
  }

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      v <- map.values[[r]]
      if(Edeg == 2) {
        resid.fit <- spaloess( resp-fitted ~ lon + lat + elev2, 
          data    = v, 
          degree  = 2, 
          span    = 0.05,
          para    = "elev2",
          family  = "symmetric",
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = "direct"),
          napred = FALSE
        )
      } else if(Edeg == 1) {
        resid.fit <- spaloess( resp-fitted ~ lon + lat + elev2, 
          data    = v, 
          degree  = 2, 
          span    = 0.05,
          drop    = "elev2",
          para    = "elev2",
          family  = "symmetric",
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = "direct"),
          napred = FALSE
        )
      } else if (Edeg == 0) {
        resid.fit <- spaloess( resp-fitted ~ lon + lat, 
          data    = v, 
          degree  = 2, 
          span    = 0.05,
          family  = "symmetric",
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = "direct"),
          napred = FALSE
        )
      }
      if (Edeg != 0) {
        grid.fit <- predloess(
          object = resid.fit,
          newdata = data.frame(lon = new.grid$lon, lat = new.grid$lat, elev2 = new.grid$elev2)
        )
      } else {
        grid.fit <- predloess(
          object = resid.fit,
          newdata = data.frame(lon = new.grid$lon, lat = new.grid$lat)
        )        
      }
      new.grid$smooth <- grid.fit
      rhcollect(map.keys[[r]], new.grid)
    })
  })
  job$parameters <- list(
    new.grid = new.grid,
    Edeg = Edeg
  )
  job$input <- rhfmt(FileInput, type = "sequence")
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bymonth.fit.plot", family, surf, Edeg, paste("sp", span, ".contour", sep="")), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 50,  #cdh3,4
    mapreduce.job.reduces = 50  #cdh5
  )
  job$setup <- expression(
    map = {library(Spaloess, lib.loc=lib.loc)}
  )
  job$readback <- TRUE
  job$jobname <- file.path(
    rh.root, par$dataset, "a1950", "bymonth.fit.plot", family, surf, Edeg, paste("sp",span, ".contour", sep="")
  )
  job.mr <- do.call("rhwatch", job)

  year <- data.frame(do.call(rbind, lapply(job.mr, `[[`, 1)))
  time <- arrange(mutate(year, idx = row.names(year)), X1, match(X2, month.abb))
  
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.spaResidcontour.bytime.ps"),
    color = TRUE,
    paper = "letter"
  )
  for(i in 1:576){
    b <- levelplot( smooth ~ lon * lat
      , data = job.mr[[as.numeric(time$idx[i])]][[2]]
      , region = TRUE
      , at = atlevels
      , col.regions = colorRampPalette(c("blue", "yellow","red"))
      , xlab = list(label="Longitude", cex=1.5)
      , ylab = list(label="Latitude", cex=1.5)
      , xlim = c(-125, -66.5)
      , scale = list(cex=1.2)
      , sub = paste(time[i, 1], time[i, 2])
      , panel = function(x, y, z, ...) {
          panel.levelplot(x,y,z,...) 
          panel.polygon(us.map$x,us.map$y, border = "black")
        }
    )
    print(b)
  }
  dev.off()

  ## This last 6 plots are scatter plot of residuals against one of 
  ## spatial factor conditional on another
  rst <- rhread(FileInput)
  tmp <- do.call(rbind, lapply(rst,"[[",1)) %>% data.frame(stringsAsFactors=FALSE) 
  tmp$X2 <- tmp$X2 %>% substr(1,3) %>% match(month.abb)
  mod <- tmp %>% with(tmp[order(X1, X2),]) %>% row.names() %>% as.numeric()
  
  if (Edeg == 0) {
    spaPara <- data.frame(permutations(2, 2, c("lon","lat")))
  } else {
    spaPara <- data.frame(permutations(3, 2, c("lon","lat","elev")))
  }
  
  for (ii in 1:nrow(spaPara)) {
    against <- spaPara[ii, 1]
    condit <- spaPara[ii, 2] 
    varname <- grep(condit, c("Latitude","Longitude","Elevation"), ignore.case=TRUE, value=TRUE)
    xlab <- grep(against, c("Latitude","Longitude","Elevation"), ignore.case=TRUE, value=TRUE)
    if (against == "elev") {
      against <- "elev2"
      xlab <- "Log Base 2 (Elevation + 128)"
    }
    trellis.device(
      device = postscript, 
      file = file.path(local.root, "output", paste("a1950.spaResid.vs", against, condit, "ps", sep=".")),
      color = TRUE, 
      paper = "letter"
    )
    for (r in mod) {
      b <- xyplot( resp - fitted ~ get(against) | equal.count(get(condit), 20, overlap=0)
        , data = rst[[r]][[2]]
        , strip=strip.custom(var.name = varname, strip.levels=rep(FALSE, 2))
        , pch = 16
        , cex = 0.3
        , scale = list(
            y = list(relation = "free", tick.number = 5, cex=1.2),
            x = list(relation = "free", tick.number = 3, cex=1.2)
          )
        , layout = c(5,4)
        , xlab = list(label=xlab, cex=1.5)
        , ylab = list(label="Residual of Spatail Imputing", cex=1.5)
        , sub = paste(rst[[r]][[1]][1], rst[[r]][[1]][2])
        , panel = function(x,y,...) {
            panel.abline(h=0, lwd=0.5, col="black")
            panel.xyplot(x,y,...)
            panel.loess(x,y, span=0.25, degree=1, col=col[2], evaluation=100,...)
        }
      )
      print(b)
    }
    dev.off()
  }

}


#########################################################
##  Diagnostic plots for components from stl2 fitting  ## 
#########################################################
#########################################
##  fitted value and obs against time  ##
#########################################
a1950.stlFitRaw <- function(data=rst, St.num = 128, test = TRUE){

  data$factor <- factor(
    x = rep( rep(paste("Period", 1:6), c(rep(108,5), 36)), times=St.num),
    levels = paste("Period", c(6:1))
  )
  data$time <- c(rep(0:107, times = 5), 0:35) 

  stations <- unique(data$station.id)

  if(test) {
    stations <- stations[1]
  }

  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("fitted.time", "a1950", "tmax", "ps", sep=".")),
    color = TRUE, 
    paper = "letter"
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
