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


imputeCrossValid <- function(input, Edeg = TRUE) {

  if (Edeg) {
    layout <- c(2,1)
    fo <- ~ mse | factor(degree)
    rst1 <- rhread(file.path(input, "1", "MABSE"))[[1]][[2]]
    rst2 <- rhread(file.path(input, "2", "MABSE"))[[1]][[2]]


    rst <- rbind(rst1, rst2)
    rst$degree <- rep(c(1,2), each = nrow(rst1))
    sub <- subset(rst, span %in% c(0.005, 0.015, 0.035, 0.055, 0.085) & mse <=10)

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
      , scale = list(cex=1.2, y=list(relation="same"))
      , panel = function(x,...) {
          panel.abline(h=seq(0,1,0.2),v=seq(0,1,0.2), col="lightgray")
          if(max(x)<5) {
            #panel.abline(h=seq(0,3,0.5),v=seq(0,1,0.2), col="lightgray")
          } else {
            panel.abline(h=seq(0,10,2),v=seq(0,1,0.2), col="lightgray")
          }
          panel.qqmath(x,...)
        }
    )
    print(b)
    dev.off()

  } else {
    rst <- rhread(file.path(input, "0", "MABSE"))[[1]][[2]]
    layout <- c(1,1)
    fo <- ~ mse
  }

    sub <- subset(rst, span %in% c(0.003, 0.005, 0.015, 0.035, 0.055, 0.085) & mse <=10)

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
      , scale = list(cex=1.2, y=list(relation="same"))
      , panel = function(x,...) {
          panel.abline(h=seq(0,1,0.2), v=seq(0,1,0.2), col="lightgray")
          panel.qqmath(x,...)
        }
    )
    print(b)
    dev.off()

} 


###################################################################
##  Visualization plots for the spatial loess fit for missing    ##
##  value imputing.                                              ##
##  input files are /a1950/bymonth.fit                           ##
###################################################################
Qsample <- function(x, num) {

  xx <- x[!is.na(x)]
  xx <- sort(xx)
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

a1950.spaImputeVisual <- function(input, plotEng, name, sample=FALSE, multiple=FALSE, byvari=FALSE) {

  output <- paste(input, "plot", sep=".")

  job <- list()
  if (multiple) {
    job$map <- expression({
      lapply(seq_along(map.keys), function(r) {
        if (!sample) {
          map.values[[r]]$month <- map.keys[[r]][2]
          rhcollect(as.numeric(map.keys[[r]][1]), map.values[[r]])    
        } else {
          if(byvari) {
            value <- Qsample(with(subset(map.values[[r]], elev<=58|elev>=1647), resp - fitted), 1000)
            value$flag <- 0
            value2 <- Qsample(with(subset(map.values[[r]], elev>58&elev<1647), resp - fitted), 1000)
            value2$flag <- 1
            value <- rbind(value, value2)
          } else {
            value <- Qsample(with(map.values[[r]], resp - fitted), 1000)
          }
          value$month <- map.keys[[r]][2]
          rhcollect(as.numeric(map.keys[[r]][1]), value)
        }
      })
    })
  } else {
    job$map <- expression({
      lapply(seq_along(map.keys), function(r){
        index <- (as.numeric(map.keys[[r]][1]) - 1950)*12 + match(map.keys[[r]][2], month.abb)
        pp <- plotEng(map.values[[r]], map.keys[[r]][1], map.keys[[r]][2])
        rhcollect(index, serialize(pp, NULL))
      })
    })
  }
  if (!is.null(multiple)) {
    job$reduce <- expression(
      pre = {
        combine <- data.frame()
      },
      reduce = {
        combine <- rbind(combine, do.call("rbind", reduce.values))
      },
      post = {
        pp <- plotEng(combine, reduce.key, byvari)
        rhcollect(reduce.key, serialize(pp, NULL))
      }
    )
  }
  job$setup <- expression(
    map = {
      library(lattice)
      library(plyr)
    },
    reduce = {
      library(lattice)
      library(plyr)
    }
  )
  job$parameters <- list(
    plotEng = plotEng,
    multiple = multiple,
    byvari = byvari,
    col = col,
    ylab = ylab,
    Qsample = Qsample,
    sample = sample
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(output, type = "sequence")
  job$mapred <- list(
    mapred.reduce.tasks = 20,  #cdh3,4
    mapreduce.job.reduces = 20  #cdh5
  )
  job$readback <- TRUE
  job$jobname <- output
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

plotEng.impResvsElev <- function(data, year, month) {

  b <- xyplot( resp-fitted ~ elev2
    , data = data
    , xlab = list(label = "Log Base 2 (Elevation + 128)", cex = 1.5)
    , ylab = list(label = "Residual of Spatail Imputing", cex = 1.5)
    , sub = paste(year, month)
    , scales = list(cex = 1.2)
    , panel = function(x,y,...) {
        panel.abline(h=0, lwd=0.5, col="black")
        panel.xyplot(x,y,, type="p", cex=0.5, pch=16, ...)
      }
  )
  return(b)

}

plotEng.QQimpRes <- function(data, year, by) {
  
  if(by) {
    fo <- resid ~ qnorm(fv) | factor(flag)*factor(month, levels=month.abb)
  } else {
    fo <- resid ~ qnorm(fv) | factor(month, levels=month.abb)
  }
  b <- xyplot( fo
    , data = data
    , layout = c(8, 3)
    , sub = paste("Year", year)
    , ylab = list(label="Spatial Residuals", cex=1.5)
    , xlab = list(label="Normal quantiles", cex=1.5)
    , scale = list(cex=1.2)
    , panel = function(x, y, ...) {
        #panel.qqmathline(y,y=y,...)
        panel.xyplot(x, y, col = col[1], pch = 16, cex =0.5, ...)
    }
  )
  return(b)

}

############################################################
##  dotplot for the range of each levels of 20 equal cut  ##
##  intervel for each of three spatial factors.           ##
############################################################
spacutRange <- function(input) {

  output <- paste(input, "plot", sep=".")

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      value <- data.frame()
      for (vari in c("lon","lat","elev")) {
        inval <- unlist(attributes(equal.count(map.values[[r]][, vari], 20, overlap=0))$level)
        month <- (as.numeric(map.keys[[r]][1]) -1950)*12 + match(map.keys[[r]][2], month.abb)
        tmp <- data.frame(
          cut=inval, type=rep(c("low","up"), times=20), 
          idx=rep(1:20, each=2), month=month, vari=vari
        )
        value <- rbind(value, tmp)
      }
      rhcollect(1, value)
    })
  })
  job$reduce <- expression(
    pre = {
      combine <- data.frame()
      value <- data.frame()
    },
    reduce = {
      combine <- rbind(combine, do.call("rbind", reduce.values))
    },
    post = {
      value <- ddply(
        .data = combine,
        .vari = c("vari","type", "idx"),
        .fun = summarize,
        cut = mean(cut)
      )
      rhcollect(reduce.key, value)
    }
  )
  job$setup <- expression(
    map = {
      library(lattice)
      library(plyr)
    },
    reduce = {
      library(plyr)
    }
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(output, type = "sequence")
  job$mapred <- list(
    mapred.reduce.tasks = 20,  #cdh3,4
    mapreduce.job.reduces = 20  #cdh5
  )
  job$readback <- TRUE
  job$jobname <- output
  job.mr <- do.call("rhwatch", job)
  
  tmp <- ddply(
    .data = job.mr[[1]][[2]],
    .vari = c("vari","idx"),
    .fun = function(r) {
      r$label <- rep(paste(round(min(r$cut),1), round(max(r$cut),1), sep="~"),2)
      r
    }
  )  

  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "equalcut.ps"),
    color = TRUE, 
    paper = "letter"
  )
  dotplot( factor(label, levels = unique(label)) ~ cut | vari
    , data = tmp
    , group = type
    , auto.key = TRUE
    , layout=c(1,1)
    , cex=1
    , pch = 16
    , key=list(
        cex = 1.2,
        text = list(label=c("lower bound", "upper bound")), 
        lines = list(pch=16, cex=1, type="p", col=col[1:2]),
        columns=2
      )
    , scale = list(relation="free", cex=1.2)
  )
  dev.off()

}


############################################################
##  Diagnostic plots for components from stlplus fitting  ## 
############################################################
a1950.STLvisual <- function(paras, input, plotEng, name, sample = TRUE, Alloutlier=TRUE, multiple=NULL){

  output <- file.path(
    rh.root, par$dataset, "a1950", "STL.plot", paste("t",paras$tw, "td", paras$td, "_s", paras$sw, 
    "sd", paras$sd, "_f", paras$fcw, "fd", paras$fcd, sep="")
  )
  job <- list()
  if (!is.null(multiple)) {
    job$map <- expression({
      lapply(seq_along(map.keys), function(r) {
        if (sample) {
          if (map.keys[[r]] %in% sample.a1950$station.id) {
            index <- as.numeric(sample.a1950$leaf[which(sample.a1950$station.id==map.keys[[r]])])
            num <- multiple[1] * multiple[2]
            map.values[[r]]$leaf <- index
            rhcollect(ceiling((index-0.5)/num), map.values[[r]])
          }
        } else {
          if (map.keys[[r]] %in% outliers.a1950.stations) {
            num <- multiple[1] * multiple[2]
            index <- as.numeric(which(outliers.a1950.stations==map.keys[[r]]))
            map.values[[r]]$station.id <- map.keys[[r]]
            rhcollect(ceiling((index-0.5)/num), map.values[[r]])
          }
        }        
      })
    })
  } else {
    job$map <- expression({
      lapply(seq_along(map.keys), function(r){
        if (sample) {
          if (map.keys[[r]] %in% sample.a1950$station.id) {
            index <- sample.a1950$leaf[which(sample.a1950$station.id==map.keys[[r]])]
            pp <- plotEng(map.values[[r]], map.keys[[r]], index)
            rhcollect(as.numeric(index), serialize(pp, NULL))
          }
        } else {
          if (map.keys[[r]] %in% outliers.a1950.stations) {
            pp <- plotEng(map.values[[r]], map.keys[[r]])
            rhcollect(map.keys[[r]], serialize(pp, NULL))
          }
        }
      })
    })
  }
  if (!is.null(multiple)) {
    job$reduce <- expression(
      pre = {
        combine <- data.frame()
      },
      reduce = {
        combine <- rbind(combine, do.call("rbind", reduce.values))
      },
      post = {
        pp <- plotEng(combine, multiple)
        rhcollect(reduce.key, serialize(pp, NULL))
      }
    )
  }
  if (Alloutlier) {
    job$shared <- c(
      file.path(rh.root, par$dataset, "a1950", "Rdata", "sample.a1950.RData"), 
      file.path(rh.root, par$dataset, "a1950", "Rdata", "outliers.a1950.RData")
    )
  } else {
    job$shared <- c(
      file.path(rh.root, par$dataset, "a1950", "Rdata", "sample.a1950.RData"), 
      file.path(rh.root, par$dataset, "a1950", "Rdata", "outliersTop.a1950.RData")
    )
  }
  job$setup <- expression(
    map = {
      load("sample.a1950.RData")
      if(Alloutlier) {
        load("outliers.a1950.RData")
        rhcounter("setup", "mapall", 1)
      } else {
        load("outliersTop.a1950.RData")
        rhcounter("setup", "mapTop", 1)
      }
      library(lattice)
      library(plyr)
    },
    reduce = {
      library(plyr)
      library(lattice)
    }
  )
  job$parameters <- list(
    sample = sample,
    plotEng = plotEng,
    col = col,
    ylab = ylab,
    Alloutlier = Alloutlier,
    multiple = multiple
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(output, type = "sequence")
  job$mapred <- list(
    mapred.reduce.tasks = 20,  #cdh3,4
    mapreduce.job.reduces = 20  #cdh5
  )
  job$readback <- TRUE
  job$jobname <- output
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

plotEng.raw <- function(data, station, leaf) {

  data$factor <- factor(
    x = rep(paste("Period", 1:6), c(rep(108,5), 36)),
    levels = paste("Period", c(6:1))
  )
  data$time <- c(rep(0:107, times = 5), 0:35) 

  b <- xyplot( resp ~ time | factor
    , data = data
    , xlab = list(label = "Month", cex = 1.5)
    , ylab = list(label = ylab, cex = 1.5)
    , sub = list(label=paste("Station ", station, "from cell", leaf), cex=1.2)
    , layout = c(1,6)
    , strip = FALSE,
    , xlim = c(0, 107)
    , ylim = c(min(c(data$resp, data$fitted), na.rm=TRUE), max(c(data$resp, data$fitted), na.rm=TRUE))
    , key=list(
        cex = 1.2,
        text = list(label=c("spatial smoothed value", "temporal fitted value")), 
        lines = list(pch=16, cex=0.7, lwd=1.5, type=c("p","l"), col=col[c(1:2)]),
        columns=2
      )
    , scales = list(
        y = list(tick.number=4, cex=1.2), 
        x = list(at=seq(0, 107, by=12), relation='same', cex=1.2)
      )
    , panel = function(x,y,subscripts,...) {
        panel.abline(v=seq(0,108, by=12), color="lightgrey", lty=3, lwd=0.5)
        fit <- subset(data[subscripts,], flag == 0)
        obs <- subset(data[subscripts,], flag == 1)
        panel.xyplot(obs$time, obs$resp, type="p", col=col[1], pch=16, cex=0.5, ...)
        panel.xyplot(fit$time, fit$fitted, type="p", col=col[1], pch=16, cex=0.5, ...)
        if (!any(grepl("fc", names(data)))) {
          panel.xyplot(data[subscripts,]$time, (data[subscripts,]$trend+data[subscripts,]$seasonal), type="l", col=col[2], lwd=1, ...)            
        } else {
          panel.xyplot(data[subscripts,]$time, (data[subscripts,]$data.seasonal+data[subscripts,]$fc.first+data[subscripts,]$fc.second), type="l", col=col[2], lwd=1, ...)
        }
      }
  )
  return(b)

}

plotEng.trend <- function(data, layout) {

  tmean <- ddply(
    .data = data,
    .variables = c("leaf","year"),
    .fun = function(r) {
      data.frame(mean = mean(r$trend + r$seasonal + r$remainder))
    }
  )
  tmv <- ddply(
    .data = arrange(tmean, year),
    .variables = "leaf",
    .fun = summarise,
    mv = as.numeric(filter(mean, rep(1/5,5), sides=2)),
    year = year
  )
  trendmean <- merge(x=data, y=tmv, by=c("leaf","year"), all.x=TRUE)

  b <- xyplot( trend ~ date | factor(leaf)
    , data = trendmean
    , xlab = list(label = "Month", cex = 1.5)
    , ylab = list(label = ylab, cex = 1.5)
    , layout = layout
    , key=list(
        text = list(label=c("trend","moving average of yearly mean")), 
        lines = list(pch=16, cex=1, lwd=2, type=c("l","p"), col=col[1:2]),
        columns=2
      )
    , prepanel = function(x,y,subscripts,...){ 
        v <- trendmean[subscripts,] 
        ylim <- range(v$mv, na.rm=TRUE) 
        ans <- prepanel.default.xyplot(v$date, v[, "trend"], ...) 
        ans$ylim <- range(c(ans$ylim, ylim), na.rm=TRUE) 
        ans 
      } 
    , scales = list(
        y = list(relation = 'free', cex=1.2, tick.number = 3), 
        x = list(at=c(0, 144, 288, 432, 576), relation = 'same', cex=1.2)
      )
    , panel = function(x, y, subscripts, ...) {
        v <- trendmean[subscripts,]
        panel.xyplot(v$date[seq(1,576,12)], v$mv[seq(1,576,12)], type="p", pch=16, cex=0.5, col=col[2], ...)
        panel.xyplot(x, y, type="l", lwd=2, col=col[1], ...)
      }
  )
  return(b)

}

plotEng.periodicseason <- function(data, layout) {
  
  mmean <- ddply(
    .data = data,
    .variable = c("leaf","month"),
    .fun = summarise,
    mean = mean(seasonal)
  )
  mmean <- arrange(mmean, match(month, month.abb))

  b <- xyplot( mean ~ match(month, month.abb) | factor(leaf)
    , data = mmean
    , xlab = list(label = "Month", cex = 1.5)
    , ylab = list(label = ylab, cex = 1.5)
    , layout = layout
    , scales = list(
        y = list(at=c(0, seq(-15,15,10)), cex=1.2), 
        x = list(at=c(1, 3, 5, 7, 9, 11), relation='same', cex=1.2)
      )
    , panel = function(x,y,...) {
        panel.xyplot(x,y,, type="b", cex=0.5, pch=16, ...)
      }
  )
  return(b)

}

plotEng.searemainder <- function(data, station, leaf) {

  smean <- ddply(
    .data = data,
    .variables = c("month"),
    .fun = function(r) {
      data.frame(mean = mean(r[, "seasonal"]))
    }
  )
  searmd <- arrange(merge(x=data, y=smean, by="month", all.x=TRUE), date)

  b <- xyplot( (seasonal - mean + remainder) ~ as.numeric(year) | factor(month, levels = month.abb),
    , data = searmd
    , xlab = list(label = "Year", cex = 1.5)
    , ylab = list(label = ylab, cex = 1.5)
    , sub = list(label=paste("Station ", station, "from cell", leaf), cex=1.2)
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
        panel.xyplot(as.numeric(searmd$year)[subscripts], (searmd$seasonal[subscripts]-searmd$mean[subscripts]), type="l", lwd=2, col=col[2], ...)
      }
  )
  return(b)

}

plotEng.QQremaider <- function(data, layout) {

  b <- qqmath( ~ remainder | factor(leaf)
    , data = data
    , distribution = qnorm
    , aspect = 1
    , layout = layout
    , xlab = list(label="Unit normal quantile", cex=1.5)
    , ylab = list(label=ylab, cex=1.5)
    , scales = list(x = list(cex=1.2), y = list(at=c(0,seq(-15,15,10)), cex=1.2))
    , prepanel = prepanel.qqmathline
    , panel = function(x, y,...) {
        panel.abline(h= c(seq(-15,15,by=10),0), v=c(-2,0,2), lty=1, lwd=1, col="lightgrey")
        panel.qqmathline(x, y=x)
        panel.qqmath(x, y, pch=16, cex=0.3, ...)
      }
  )
  return(b)

}

plotEng.remainderDate <- function(data, station, leaf) {
 
  b <- xyplot( remainder ~ date
    , data = data
    , xlab = list(label = "Month", cex=1.5)
    , ylab = list(label = ylab, cex=1.5)
    , sub = list(label=paste("Station ", station, "from cell", leaf), cex=1.2)
    , scale = list(x = list(at=seq(0, 576, 96)), cex=1.2)
    , key = list(
        text=list(label=c("remainder","loess smoothing")), 
        lines=list(pch=1, cex=1, lwd=2, type=c("p","l"), col=col[1:2]), 
        columns=2
      )
    , panel = function(x,y,...){
        panel.abline(h=0, color="black", lty=1)
        panel.xyplot(x,y, cex = 0.6, pch=1, ...)
        panel.loess(x,y, span=0.8, degree=1, lwd = 2, col=col[2],...)
      }
  )

}

plotEng.remainderMonth <- function(data, station, leaf) {

  b <- xyplot( remainder ~ as.numeric(year) | factor(month, levels = month.abb)
    , data = data
    , xlab = list(label = "Year", cex = 1.5)
    , ylab = list(label = ylab, cex = 1.5)
    , sub = list(label=paste("Station ", station, "from cell", leaf), cex=1.2)
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
        panel.loess(x,y, span=0.8, degree=1, lwd = 2, col=col[2],...)
      }
  )
  return(b)

}

plotEng.remainderMonth2 <- function(data, station, leaf) {

  b <- xyplot( remainder ~ as.numeric(year) | factor(month, levels = month.abb)
    , data = data
    , xlab = list(label = "Year", cex = 1.5)
    , ylab = list(label = ylab, cex = 1.5)
    , layout = c(2,6)
    , sub = list(label=paste("Station ", station, "from cell", leaf), cex=1.2)
    , scales = list(
        y = list(cex=1.2, at=seq(-10,10,5)), 
        x = list(tick.number=10, relation='same', cex=1.2)
      )
    , panel = function(x,y,...) {
        panel.abline(h=0, color="black", lty=1)
        panel.xyplot(x, y, type="b", pch=16, cex=0.5, ...)
      }
  )
  return(b)

}

plotEng.remainderACF <- function(data, layout) {

  alpha <- 0.95

  ACF <- ddply(
    .data = data,
    .variables = "leaf",
    .fun = function(r) {
      corr <- acf(r[, "remainder"], plot=FALSE)
      data.frame(
        correlation = corr$acf,
        lag = corr$lag 
      )
    }
  )
  clim <- qnorm((1 + alpha)/2)/sqrt(576)

  b <- xyplot( correlation ~ lag | factor(leaf)
    , data = ACF
    , subset = lag != 0
    , layout = layout
    , xlab = list(label = "Lag", cex = 1.5)
    , ylab = list(label = "Autocorrelation Function", cex = 1.5)
    , scales = list(
        y = list(cex=1.2), 
        x = list(relation='same', cex=1.2)
      )
    , panel = function(x,y,...) {
        panel.abline(h=0)
        panel.xyplot(x,y, type="h", lwd=1.5,...)
        panel.abline(h=c(-clim, clim), lty=2, col="red",...)
      }
  )
  return(b)

}


#############################################
##  Visualize the final residuals by month ##
#############################################
a1950.spafitVisualMon <- function(input, plotEng, vars, target) {

  output <- paste(input, "plot", sep=".")

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      year <- map.keys[[r]][1]
      month <- map.keys[[r]][2]
      pp <- plotEng(subset(map.values[[r]], !is.na(resp)), year, month, vars, target)
      index <- (as.numeric(year) - 1950)*12 + match(month, month.abb)
      rhcollect(index, serialize(pp, NULL))
    })
  })
  job$setup <- expression(
    map = {
      library(lattice)
    }
  )
  job$parameters <- list(
    plotEng = plotEng,
    col = col,
    ylab = ylab,
    target = target,
    vars = vars
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(output, type = "sequence")
  job$mapred <- list(
    mapred.reduce.tasks = 20,  #cdh3,4
    mapreduce.job.reduces = 20  #cdh5
  )
  job$readback <- TRUE
  job$jobname <- output
  job.mr <- do.call("rhwatch", job)
  
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("a1950", target, "vs", vars[1], vars[2], "ps", sep=".")),
    color = TRUE, 
    paper = "letter"
  )
  for(ii in order(unlist(lapply(job.mr, "[[", 1)))) {
    print(unserialize(job.mr[[ii]][[2]]))
  }
  dev.off()  

}

plotEng.RevsSpa <- function(data, year, month, vars, target) {

  if (target == "remainder"){
    ylab <- "Remainder of Stlplus"
  } else {
    ylab <- "Residual"
    data$spaResid <- with(data, remainder - spafit)
  }
  against <- vars[1, 1]
  condit <- vars[1, 2] 
  varname <- grep(condit, c("Latitude","Longitude","Elevation"), ignore.case=TRUE, value=TRUE)
  xlab <- grep(against, c("Latitude","Longitude","Elevation"), ignore.case=TRUE, value=TRUE)
  if (against == "elev") {
    against <- "elev2"
    xlab <- "Log Base 2 (Elevation + 128)"
  }
  if (condit == "elev"){
    data$elev <- 2^data$elev2 - 128
  }
  b <- xyplot( get(target) ~ get(against) | equal.count(get(condit), 20, overlap=0)
    , data = data
    , strip=strip.custom(var.name = varname, strip.levels=rep(FALSE, 2))
    , pch = 16
    , cex = 0.3
    , scale = list(
        y = list(relation = "free", tick.number = 5, cex=1.2),
        x = list(relation = "free", tick.number = 3, cex=1.2)
      )
    , layout = c(5,4)
    , xlab = list(label=xlab, cex=1.5)
    , ylab = list(label=ylab, cex=1.5)
    , sub = paste(year, month)
    , panel = function(x,y,...) {
        panel.abline(h=0, lwd=0.5, col="black")
        panel.xyplot(x,y,...)
        panel.loess(x,y, span=0.25, degree=1, col=col[2], evaluation=100,...)
    }
  )
  return(b)

}

plotEng.RevsFit <- function(data, year, month, vars=NULL, target=NULL) {

  b <- xyplot( remainder-spafit ~ spafit
    , data = data
    , pch = 16
    , cex = 0.5
    , scale = list(
        y = list(cex=1.2),
        x = list(cex=1.2)
      )
    , ylab = list(label="Residual", cex=1.5)
    , xlab = list(label="Spatial Fitted Value", cex=1.5)
    , sub = paste(year, month)
    , panel = function(x,y,...) {
        panel.abline(h=0, lwd=0.5, col="black")
        panel.xyplot(x,y,...)
        #panel.loess(x,y, span=0.25, degree=1, col=col[2], evaluation=100,...)
    }
  )
  return(b)

}


###############################################
##  Visualize the final residuals by station ##
###############################################
a1950.spafitVisualStat <- function(input, plotEng, name, sample = TRUE, multiple=NULL) {

  output <- paste(input, "plot", sep=".")

  job <- list()
  if (!is.null(multiple)) {
    job$map <- expression({
      lapply(seq_along(map.keys), function(r) {
        if (sample) {
          if (map.keys[[r]] %in% sample.a1950$station.id) {
            index <- as.numeric(sample.a1950$leaf[which(sample.a1950$station.id==map.keys[[r]])])
            num <- multiple[1] * multiple[2]
            map.values[[r]]$leaf <- index
            rhcollect(ceiling((index-0.5)/num), subset(map.values[[r]], flag==1))
          }
        } else {
          if (sum(map.values[[r]]$flag) == 576) {
            index <- as.numeric(which(a1950Nomiss==map.keys[[r]]))
            map.values[[r]]$station.id <- map.keys[[r]]
            num <- multiple[1] * multiple[2]
            rhcollect(ceiling((index-0.5)/num), map.values[[r]])
          }
        }        
      })
    })
  } else {
    job$map <- expression({
      lapply(seq_along(map.keys), function(r){
        if (sample) {
          if (map.keys[[r]] %in% sample.a1950$station.id) {
            index <- sample.a1950$leaf[which(sample.a1950$station.id==map.keys[[r]])]
            pp <- plotEng(map.values[[r]], map.keys[[r]], index)
            rhcollect(as.numeric(index), serialize(pp, NULL))
          }
        } else {
          if (sum(map.values[[r]]$flag) == 576) {
            pp <- plotEng(map.values[[r]], map.keys[[r]], NULL)
            rhcollect(map.keys[[r]], serialize(pp, NULL))
          }
        }
      })
    })
  }
  if (!is.null(multiple)) {
    job$reduce <- expression(
      pre = {
        combine <- data.frame()
      },
      reduce = {
        combine <- rbind(combine, do.call("rbind", reduce.values))
      },
      post = {
        pp <- plotEng(combine, multiple)
        rhcollect(reduce.key, serialize(pp, NULL))
      }
    )
  }
  if(sample) {
    job$shared <- c(file.path(rh.root, par$dataset, "a1950", "Rdata", "sample.a1950.RData"))
  } else {
    job$shared <- c(file.path(rh.root, par$dataset, "a1950", "Rdata", "nomiss.a1950.RData"))    
  } 
  job$setup <- expression(
    map = {
      if(sample) {load("sample.a1950.RData")} else {load("nomiss.a1950.RData")}
      library(lattice)
      library(plyr)
    },
    reduce = {
      library(plyr)
      library(lattice)
    }
  )
  job$parameters <- list(
    sample = sample,
    plotEng = plotEng,
    col = col,
    ylab = ylab,
    multiple = multiple
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(output, type = "sequence")
  job$mapred <- list(
    mapred.reduce.tasks = 20,  #cdh3,4
    mapreduce.job.reduces = 20  #cdh5
  )
  job$readback <- TRUE
  job$jobname <- output
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
plotEng.residualDate <- function(data, station, leaf) {
  
  set.seed(66)
  data$date <- sample(1:576, 576, replace=FALSE)
  b <- xyplot( remainder - spafit ~ date
    , data = data
    , group = flag
    , xlab = list(label = "Month", cex=1.5)
    , ylab = list(label = ylab, cex=1.5)
    , sub = list(label=paste("Station ", station), cex=1.2)
    , scale = list(x = list(at=seq(0, 576, 96)), cex=1.2)
    , key = list(
        text=list(label=c("remainder","loess smoothing")), 
        lines=list(pch=1, cex=1, lwd=2, type=c("p","l"), col=col[1:2]), 
        columns=2
      )
    , panel = function(x,y,...){
        panel.abline(h=0, color="black", lty=1)
        panel.xyplot(x,y, cex = 0.6, pch=1, ...)
        panel.loess(x,y, span=0.75, degree=1, lwd = 2, col=col[2],...)
      }
  )

}
plotEng.residACF <- function(data, layout) {

  alpha <- 0.95
  
  ACF <- ddply(
    .data = data,
    .variables = "station.id",
    .fun = function(r) {
      corr <- acf(r$remainder - r$spafit, plot=FALSE)
      data.frame(
        correlation = corr$acf,
        lag = corr$lag 
      )
    }
  )
  clim <- qnorm((1 + alpha)/2)/sqrt(576)

  b <- xyplot( correlation ~ lag | factor(station.id)
    , data = ACF
    , subset = lag != 0
    , layout = layout
    , xlab = list(label = "Lag", cex = 1.5)
    , ylab = list(label = "Autocorrelation Function", cex = 1.5)
    , scales = list(
        y = list(cex=1.2), 
        x = list(relation='same', cex=1.2)
      )
    , panel = function(x,y,...) {
        panel.abline(h=0)
        panel.xyplot(x,y, type="h", lwd=1.5,...)
        panel.abline(h=c(-clim, clim), lty=2, col="red",...)
      }
  )
  return(b)

}
plotEng.spafitDate <- function(data, station, leaf) {

  b <- xyplot( spafit ~ date
    , data = data
    , xlab = list(label = "Month", cex=1.5)
    , ylab = list(label = ylab, cex=1.5)
    , sub = list(label=paste("Station ", station, "from cell", leaf), cex=1.2)
    , scale = list(x = list(at=seq(0, 576, 96)), cex=1.2)
    , key = list(
        text=list(label=c("remainder spatial fit","loess smoothing")), 
        lines=list(pch=1, cex=1, lwd=2, type=c("p","l"), col=col[1:2]), 
        columns=2
      )
    , panel = function(x,y,...){
        panel.abline(h=0, color="black", lty=1)
        panel.xyplot(x,y, cex = 0.6, pch=1, ...)
        panel.loess(x,y, span=0.5, degree=1, lwd = 2, col=col[2],...)
      }
  )

}
plotEng.spafitACF <- function(data, layout) {

  alpha <- 0.95
  data <- arrange(data, date)
  ACF <- ddply(
    .data = data,
    .variables = "station.id",
    .fun = function(r) {
      corr <- acf(r$spafit, plot=FALSE)
      data.frame(
        correlation = corr$acf,
        lag = corr$lag 
      )
    }
  )
  clim <- qnorm((1 + alpha)/2)/sqrt(576)

  b <- xyplot( correlation ~ lag | factor(station.id)
    , data = ACF
    , subset = lag != 0
    , layout = layout
    , xlab = list(label = "Lag", cex = 1.5)
    , ylab = list(label = "Autocorrelation Function", cex = 1.5)
    , scales = list(
        y = list(cex=1.2), 
        x = list(relation='same', cex=1.2)
      )
    , panel = function(x,y,...) {
        panel.abline(h=0)
        panel.xyplot(x,y, type="h", lwd=1.5,...)
        panel.abline(h=c(-clim, clim), lty=2, col="red",...)
      }
  )
  return(b)

}  
plotEng.spafitDateMulti <- function(data, station, leaf) {

  data$factor <- factor(
    x = rep(paste("Period", 1:6), c(rep(108,5), 36)),
    levels = paste("Period", c(6:1))
  )
  data$time <- c(rep(0:107, times = 5), 0:35) 

  b <- xyplot( spafit ~ time | factor
    , data = data
    , xlab = list(label = "Month", cex = 1.5)
    , ylab = list(label = ylab, cex = 1.5)
    , sub = list(label=paste("Station ", station, "from cell", leaf), cex=1.2)
    , layout = c(1,6)
    , strip = FALSE,
    , xlim = c(0, 107)
    , ylim = c(
        min(c(data$spafit, data$spafitseasonal+data$spafittrend), na.rm=TRUE), 
        max(c(data$spafit, data$spafitseasonal+data$spafittrend), na.rm=TRUE)
      )
    , key=list(
        cex = 1.2,
        text = list(label=c("spatial smoothed value of remainder", "temporal fitted value")), 
        lines = list(pch=16, cex=0.7, lwd=1.5, type=c("p","l"), col=col[c(1:2)]),
        columns=2
      )
    , scales = list(
        y = list(tick.number=4, cex=1.2), 
        x = list(at=seq(0, 107, by=12), relation='same', cex=1.2)
      )
    , panel = function(x,y,subscripts,...) {
        panel.abline(v=seq(0,108, by=12), color="lightgrey", lty=3, lwd=0.5)
        panel.xyplot(x, y, type="p", col=col[1], pch=16, cex=0.5, ...)
        if (!any(grepl("fc", names(data)))) {
          panel.xyplot(data[subscripts,]$time, (data[subscripts,]$spafittrend+data[subscripts,]$spafitseasonal), type="l", col=col[2], lwd=1, ...)            
        } else {
          panel.xyplot(data[subscripts,]$time, (data[subscripts,]$data.seasonal+data[subscripts,]$fc.first+data[subscripts,]$fc.second), type="l", col=col[2], lwd=1, ...)
        }
      }
  )
  return(b)

}  

residSpaFitVisl <- function(i, j, bestStlplus) {

  FileInput <- file.path(rh.root, par$dataset, "a1950", "STL.bymonth", bestStlplus)

  spaPara <- data.frame(permutations(3, 2, c("lon","lat","elev")))
  para <- list(span=i, Edeg=j, degree=2, surf="direct")
  FileOutput <- file.path(
    rh.root, par$dataset, "a1950", "STL.bymth.remfit", 
    bestStlplus, "symmetric", para$surf, para$Edeg, paste("sp", para$span, sep="")
  )
  a1950.Spatialfit(input=FileInput, output=FileOutput, argumt=para)
  a1950.spafitVisualMon(input=FileOutput, plotEng.RevsFit, vars=NULL, target=NULL)
  for (k in 1:nrow(spaPara)) {
    vars <- spaPara[k, ]
    a1950.spafitVisualMon(input=FileOutput, plotEng.RevsSpa, vars, target="spaResid")
  }
  ## the residual against time for each station which does not have missing value
  FileInput <- FileOutput 
  FileOutput <- paste(FileInput, "bystation", sep=".")
  swapTostation(FileInput, FileOutput, elevFlag = FALSE)
  FileInput <- FileOutput
  a1950.spafitVisualStat(
    input=FileInput, plotEng.residualDate, 
    name="resid.vs.time", sample = FALSE, multiple=NULL
  )
  a1950.spafitVisualStat(
    input=FileInput, plotEng.residACF, 
    name="residACF", sample = FALSE, multiple=c(3, 3)
  )  
  ## the overall quantile plot of the residual for each month
  FileInput <- file.path(
    rh.root, par$dataset, "a1950", "STL.bymth.remfit", 
    bestStlplus, "symmetric", para$surf, para$Edeg, paste("sp", para$span, sep="")
  )
  df <- a1950.residQuant(
    input=FileInput, target="residual", by=NULL, 
    probs=seq(0.005, 0.995, 0.005), nBins = 10000, tails = 100
  )
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.residual.quant.ps"),
    color = TRUE, 
    paper = "letter"
  )
    b <- xyplot(q ~ fval
      , data = df
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Residual", cex=1.5)
      , scale = list(cex=1.2)
    )
    print(b)
  dev.off()      
  df <- a1950.residQuant(
    input=FileInput, target="residual", by=NULL, 
    probs=seq(0.005, 0.995, 0.005), nBins = 10000, tails = 0
  )
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.residual.centerquant.ps"),
    color = TRUE, 
    paper = "letter"
  )
    b <- xyplot(q ~ fval
      , data = df
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Residual", cex=1.5)
      , scale = list(cex=1.2, y=list(at=seq(-2,2,1)))
      , panel = function(x, y,...) {
          panel.abline(v=seq(0,1,0.2), h=seq(-2,2,1), col="lightgray", lwd=0.5)
          panel.xyplot(x,y,...)
      }
    )
    print(b)
  dev.off()  
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", "a1950.residual.QQ.ps"),
    color = TRUE, 
    paper = "letter"
  )
    b <- xyplot(q ~ qt(fval,3)
      , data = df
      , xlab = list(label="Quantiles of t-distribution", cex=1.5)
      , ylab = list(label="Residual", cex=1.5)
      , scale = list(cex=1.2, y=list(at=seq(-2,2,1)))
      , aspect = 1
      , panel = function(x, y, ...) {
          panel.qqmathline(y,y=y, distribution=function(p) qt(p, df=3),...)
          panel.xyplot(x, y, col = col[1], pch = 1, cex = 1, ...)
      }
    )
    print(b)
  dev.off() 

}
