source("~/Rhipe/rhinitial.R")
par <- list()
par$modified <- TRUE
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$loess <- "loess04" # loess04 is w/ elevation
par$family <- "symmetric"
par$degree <- 2
par$outer <- 5
par$type <- "same" # or "same", "decr"
par$parameters <- list(
  sw = "periodic",
  tw = 109,
  sd = 1,
  td = 2,
  inner = 1,
  outer = 1
) 
if (par$type == "same") {
  par$span <- rep(0.05, 5)
} else if (par$type == "incr") {
  par$span <- c(seq(0.03, 0.05, by=0.005))
} else {
  par$span <- c(seq(0.05, 0.03, by=-0.005))
}
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")


##############################################
## check the convergence of findal residuals##
##############################################
job4 <- list()
job4$map <- expression({
  lapply(seq_along(map.values), function(r) {
    file <- Sys.getenv("mapred.input.file")
    key <- substr(unlist(strsplit(tail(strsplit(file, "/")[[1]],3)[2], "[.]")), 8, 9)
    value <- with(map.values[[r]], (fitted-spatial)^2)
    rhcollect(as.numeric(key), value)
  })
})
job4$reduce <- expression(
  pre= {
    combine <- vector()
  },
  reduce = {
    combine <- c(combine, unlist(reduce.values))
  },
  post = {
    rhcollect(reduce.key, mean(combine, na.rm = TRUE))
  }
)
job4$input <- rhfmt(
  file.path(
    rh.datadir, par$dataset, "spatial", "a1950", par$family, 
    paste("sp", "0.05", sep=""), paste(par$loess, "loop", par$type, sep="."), paste("Spatial",1:length(par$span), sep="")
  ), 
  type = "sequence"
)
job4$output <- rhfmt(
  file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
    paste("sp", "0.05", sep=""), paste(par$loess, "loop", par$type, sep="."), "residuals"
  ), 
  type = "sequence"
)
job4$mapred <- list(mapred.reduce.tasks = 8)
job4$mon.sec <- 10
job4$jobname <- file.path(
  rh.datadir, par$dataset, "spatial", "a1950", par$family, 
  paste("sp", "0.05", sep=""), paste(par$loess, "loop", par$type, sep="."), "residuals"
)  
job4$readback <- FALSE
job.mr <- do.call("rhwatch", job4)

MSE <- function(type = "same") {
  
  rst <- rhread(
    file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
      paste("sp", "0.05", sep=""), paste(par$loess, "loop", type, sep="."), "residuals"
    )
  )

  result <- data.frame(
    iter = unlist(lapply(rst, "[[", 1)), 
    MSE  = unlist(lapply(rst, "[[", 2))
  )

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "MSE", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  b <- xyplot(MSE ~ iter
    , data = result[order(result$iter),]
    , aspect = 1
    , type = "b"
    , xlab = "Iteration time"
    , ylab = "Mean Squared Error"
    , scale = list(x = list(at=c(1,5,10,15,20)))
  )
  print(b)
  dev.off()

}


#########################################
##  Normality checking of the residual ##
#########################################
Qrst <- function(x) {
  x <- x[!is.na(x)]
  a <- sort(x)
  idx <- round(seq(1, length(x), length.out = 100))
  f.value <- (idx - 0.5) / length(a)
  qnorm <- qnorm(f.value)
  value <- data.frame(
    residual = a[idx[2:99]], 
    qnorm = qnorm[2:99], 
    fv = f.value[2:99]
  )
}

job5 <- list()
job5$map <- expression({
  lapply(seq_along(map.values), function(r) {
    file <- Sys.getenv("mapred.input.file")
    key <- substr(unlist(strsplit(tail(strsplit(file, "/")[[1]],3)[2], "[.]"))[1], 6, 6)
    value <- with(map.values[[r]], tmax - trend - seasonal - spatial)
    rhcollect(as.numeric(key), value)
  })
})
job5$reduce <- expression(
  pre= {
    combine <- vector()
  },
  reduce = {
    combine <- c(combine, unlist(reduce.values))
  },
  post = {
    value <- Qrst(combine)
    rhcollect(reduce.key, value)
  }
)
job5$parameters <- list(
  Qrst = Qrst
)
job5$input <- rhfmt(
  file.path(
    rh.datadir, par$dataset, "spatial", "a1950", par$family, 
    paste("sp", "0.05", sep=""), paste(par$loess, "loopnew", par$type, sep="."), paste("Outer",1:5, ".bystation", sep="")
  ), 
  type = "sequence"
)
job5$output <- rhfmt(
  file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
    paste("sp", "0.05", sep=""), paste(par$loess, "loopnew", par$type, sep="."), "residuals.qnorm"
  ), 
  type = "sequence"
)
job5$mapred <- list(mapred.reduce.tasks = 1)
job5$mon.sec <- 10
job5$jobname <- file.path(
  rh.datadir, par$dataset, "spatial", "a1950", par$family, 
  paste("sp", "0.05", sep=""), paste(par$loess, "loopnew", par$type, sep="."), "residuals.qnorm"
)  
job5$readback <- FALSE
job.mr <- do.call("rhwatch", job5)

QQ.resid <- function(type = "same") {
  
  rst <- rhread(  
    file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
      paste("sp", "0.05", sep=""), paste(par$loess, "loopnew", type, sep="."), "residuals.qnorm"
    )
  )  

  result <- do.call(rbind, lapply(rst, "[[", 2))
  result$iter <- rep(unlist(lapply(rst, "[[", 1)), each = 98)

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "QQ.resid", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
    b <- xyplot(residual ~ qnorm | factor(iter)
      , data = result
      , pch = 1
      , cex = 0.4
      , aspect = "xy"
      , layout = c(5,1)
      , xlab = "Unit normal quantile"
      , ylab = "Residuals"
      , panel = function(x,y,...){
          panel.qqmathline(y,y=y,...)
          panel.xyplot(x,y,...)
      }
    )
    print(b)
  dev.off()

}


#################################################################
##  Check the final residual from temporal and spatial aspect  ##
#################################################################
tempCheck <- function(type = "same", iter = 5) {

  library(plyr)
  library(magrittr)

  load(file.path(local.datadir, "samplestation.a1950.RData")) ##kd-tree is built in kdfindcells.R 
  
  rst <- rhread(  
    file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
      paste("sp", "0.05", sep=""), paste(par$loess, "loopnew", type, sep="."), paste("Outer", iter, ".bystation", sep="")
    )
  )  
  
  idx <- which(lapply(rst, "[[", 1) %in% tmax.sample.a1950[order(tmax.sample.a1950$station.id), 1])
  
  stationLeaf <- data.frame(
    station.id = do.call(rbind, lapply(idx, function(r){rst[[r]][[1]]})),
    idx = idx,
    stringsAsFactors = FALSE
  ) %>% merge(tmax.sample.a1950, by="station.id")

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "comps.by.stations", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(1:128, function(i) {
    r <- with(stationLeaf, stationLeaf[leaf == i, "idx"])
    data <- rst[[r]][[2]]
    data <- data.frame(cbind(
      do.call(c, with(data, list(trend, seasonal, spatial, tmax-trend-seasonal-spatial))), 
      rep(c("trend", "seasonal", "spatial", "residual"), each=576)
    ), stringsAsFactors = FALSE)
    names(data) <- c("value", "comp")
    data$value <- as.numeric(data$value)
    tmp <- ddply(
      .data = data,
      .variable = "comp",
      .fun = summarise,
      m = mean(value)
    )
    data$mean <- rep(tmp[c(4,2,3,1),2], each = 576)
    b <- qqmath( ~ (value-mean) | factor(comp)
      , data = data
      , distribution = qunif
      , pch = 1
      , cex = 0.4
      , aspect = 2
      , layout = c(4,1)
      , main = paste("Quantiles of Components")
      , sub = paste("Station from cell", tmax.sample.a1950[tmax.sample.a1950$station.id == rst[[r]][[1]], 2]) 
      , xlab = "f-value"
      , ylab = "Maximum Temperature"
      , panel = function(x,...){
          panel.abline(h=seq(-15,15,by=3),v=seq(0,1,by=0.2), lwd=0.5, lty=1, col="lightgray")
          panel.abline(h=0, col="black", lwd=0.5)
          panel.qqmath(x,...)
      }
    )
    #print(b)
  })  
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "resid.by.stations", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(1:128, function(i) {
    r <- with(stationLeaf, stationLeaf[leaf == i, "idx"])
    data <- rst[[r]][[2]]
    data <- data[order(data$time), ]
    data$time1 <- rep(0:95, 6)
    data$factor <- factor(
      rep(1:6, each = 96),
      levels = c(6:1)
    )
    b <- xyplot((tmax - trend - seasonal - spatial) ~ time1 | factor
      , data = data
      , pch = 1
      , strip = FALSE
      , cex = 0.4
      , aspect = "xy"
      , layout = c(1,6)
      , xlim = c(-1, 96)
      , scale = list(x = list(at=seq(0, 95, by=12)))
      , main = paste("Station from cell", tmax.sample.a1950[tmax.sample.a1950$station.id == rst[[r]][[1]], 2]) 
      , xlab = "Month"
      , ylab = "Residuals"
      , panel = function(x,y,...){
          panel.abline(h=0, col="black", lwd=0.5)
          panel.xyplot(x,y,...)
      }
    )
    #print(b)
  })
  dev.off()

  result <- do.call(rbind, lapply(stationLeaf$idx, function(r){
    subset(rst[[r]][[2]], select = c(tmax, seasonal, trend, spatial, station.id, time))
  })) %>% merge(stationLeaf, by = "station.id")
  
  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "QQ.by.stations", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  b <- qqmath(~ (tmax - trend - seasonal - spatial) | factor(leaf)
    , data = result
    , pch = 16
    , cex = 0.4
    , layout = c(7,2)
    , xlab = "Unit normal quantile"
    , ylab = "Residuals"
    , panel = function(x,...) {
        panel.qqmathline(x,...)
        panel.qqmath(x,...)
    }
  )
  print(b)
  dev.off()

  result <- do.call(rbind, lapply(stationLeaf$idx, function(r){
    subset(rst[[r]][[2]], select = c(tmax, seasonal, trend, spatial, station.id, year, month, time))
  })) %>% merge(stationLeaf, by = "station.id")

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "seasonal.by.stations", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  b <- xyplot( seasonal ~ factor(month, levels=c(
    "Jan","Feb","Mar","Apr","May","June",
    "July","Aug", "Sep", "Oct", "Nov", "Dec"
    )) | factor(leaf)
    , data = result[order(result$time),]
    , subset = year == "1950"
    , pch = 16
    , cex = 0.5
    , type = "b"
    , scale = list(x=list(at=c(1, 3, 5, 7, 9, 11), relation='same'))
    , aspect = "xy"
    , layout = c(5,4)
    , xlab = "Month"
    , ylab = "Seasonal Component"
    , panel = function(x,y,...) {
        panel.xyplot(x,y,...)
    }
  )
  print(b)
  dev.off()

  mav <- function(x,n=5){  

    data.frame(mv=matrix(filter(x$mean,rep(1/n,n), sides=2), ncol=1), stringsAsFactors=FALSE)  

  }  

  dr <- result %>% ddply(
    .variable = c("station.id", "year"),
    .fun = summarise,
    mean = mean(tmax, na.rm=TRUE)
  ) %>% ddply(
    .variable = "station.id",
    .fun= mav
  )  

  result <- result[order(result$station.id, result$time),]
  result <- dr[rep(row.names(dr), each=12), "mv", drop=FALSE] %>% cbind(result)

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "trend.by.stations", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  b <- xyplot( trend ~ time | factor(leaf) 
    , data = result
    , xlab = list(label = "Month")
    , ylab = list(label = "Maximum Temperature (degrees centigrade)")
    , xlim = c(0, 576)
    , pch = 16
    , aspect = "xy"
    , layout = c(3,3)
    , key=list(
        text = list(label=c("trend component","moving average of yearly mean")), 
        lines = list(pch=16, cex=0.7, lwd=1.5, type=c("l","p"), col=col[1:2]),
        columns=2
      )
    , scales = list(
        y = list(relation = 'free'), 
        x=list(at=seq(0, 576, by=120), relation = 'same')
      )
    , panel = function(x, y, subscripts, ...) {
        sub <- result[subscripts, ]
        panel.xyplot(
          x = sub$time[seq(1, 576, by=12)],
          y = sub$mv[seq(1, 576, by=12)],
          type="p", col=col[2], cex = 0.5, ...
        )
        panel.xyplot(x, y, type="l", col=col[1], ...)
      }
  )
  print(b)
  dev.off()
  
  monthLevel <- c(
    "Jan","Feb","Mar","Apr","May","June",
    "July","Aug", "Sep", "Oct", "Nov", "Dec"
  )

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "residmonth.by.stations", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
    for(i in 1:128){
      b <- xyplot( (tmax - trend - seasonal - spatial) ~ as.numeric(year) | factor(month, levels=monthLevel)
        , data = subset(result, leaf == i)
        , xlab = list(label = "Year")
        , ylab = list(label = "Maximum Temperature (degrees centigrade)")
        , main = paste("Residuals from station of cell", i)
        , pch = 16
        , cex = 0.5
        , layout = c(12,1)
        , strip = TRUE
        , scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same'))
        , panel = function(x,y,...){
            panel.abline(h=0, color="black", lty=1)
            panel.xyplot(x,y,...)
            panel.loess(x,y,span=3/4, degree=2, col=col[2],...)
          }
      )
      print(b)
    }
  dev.off()  

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "residmonth2.by.stations", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
    for(i in 1:128){
      b <- xyplot( (tmax - trend - seasonal - spatial)  ~ as.numeric(year) | factor(month, levels=monthLevel)
        , data = subset(result, leaf == i)
        , xlab = list(label = "Year")
        , ylab = list(label = "Maximum Temperature (degrees centigrade)")
        , main = paste("Residuals from station of cell", i)
        , type = "b"      
        , pch = 16
        , cex = 0.5
        , layout = c(2,6)
        , strip = TRUE
        , scales = list(y = list(relation = 'same', alternating=TRUE), x=list(tick.number=10, relation='same'))
        , panel = function(x,y,...) {
            panel.abline(h=0, color="black", lty=1)
            panel.xyplot(x,y,...)
          }
      )
      print(b)
    }
  dev.off()

  ACF <- ddply(
    .data = result,
    .variables = "leaf",
    .fun = summarise,
    correlation = c(acf((tmax - trend - seasonal - spatial), plot=FALSE)$acf),
    lag = c(acf((tmax - trend - seasonal - spatial), plot=FALSE)$lag) 
  )  

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "residACF.by.stations", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )     
    for(i in 1:128){
      b <- xyplot( correlation ~ lag
        , data = subset(ACF, leaf==i & lag!=0)
        , xlab = list(label = "Lag")
        , ylab = list(label = "ACF")
        , main = list(label = paste("Station from cell", i))
        , type = "h"
        , panel = function(x,y,...) {
            panel.abline(h=0)
            panel.xyplot(x,y,...)
          }
      )
      print(b)
    }
  dev.off()

}

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

spatialCheck <- function(type = "same", iter = 5) {
  
  library(plyr)
  library(magrittr)

  rst <- rhread(  
    file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
      paste("sp", "0.05", sep=""), paste(par$loess, "loopnew", type, sep="."), paste("Spatial", iter, sep="")
    )
  )  
  
  tmp <- do.call(rbind, lapply(rst,"[[",1)) %>%
    data.frame(stringsAsFactors=FALSE) 
  tmp$X2 <- tmp$X2 %>% substr(1,3) %>% match(month.abb)
  mod <- tmp %>% with(tmp[order(X1, X2),]) %>% row.names() %>% as.numeric()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "comps.by.month", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    data <- rst[[r]][[2]]
    data <- data.frame(cbind(
      do.call(c, with(data, list(trend, seasonal, spatial, tmax-trend-seasonal-spatial))), 
      rep(c("trend", "seasonal", "spatial", "residual"), each=7738)
    ), stringsAsFactors = FALSE)
    names(data) <- c("value", "comp")
    data$value <- as.numeric(data$value)
    tmp <- ddply(
      .data = data,
      .variable = "comp",
      .fun = summarise,
      m = mean(value, na.rm = TRUE)
    )
    data$mean <- rep(tmp[c(4,2,3,1),2], each = 7738)
    b <- qqmath( ~ (value-mean) | factor(comp)
      , data = data
      , subset = (value-mean) > -50
      , distribution = qunif
      , pch = 1
      , cex = 0.4
      , aspect = 2
      , layout = c(4,1)
      , main = paste("Quantiles of Components for", rst[[r]][[1]][1], rst[[r]][[1]][2])
      , xlab = "f-value"
      , ylab = "Maximum Temperature"
      , panel = function(x,...){
          panel.abline(h=seq(-15,15,by=3),v=seq(0,1,by=0.2), lwd=0.5, lty=1, col="lightgray")
          panel.abline(h=0, col="black", lwd=0.5)
          panel.qqmath(x,...)
      }
    )
    #print(b)
  })  
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "QQ.resid.by.month", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    data <- rst[[r]][[2]]
    data <- data.frame(cbind(
      do.call(c, with(data, list(tmax-trend-seasonal, tmax-trend-seasonal-spatial))), 
      rep(c("residual+spatial", "residual"), each=7738)
    ), stringsAsFactors = FALSE)
    names(data) <- c("value", "comp")
    data$value <- as.numeric(data$value)
    b <- qqmath( ~ value | factor(comp)
      , data = data
      , distribution = qnorm
      , pch = 1
      , cex = 0.4
      , aspect = 1
      , layout = c(2,1)
      , main = paste("Quantiles of Components for", rst[[r]][[1]][1], rst[[r]][[1]][2])
      , xlab = "Unit normal quantile"
      , ylab = "Maximum Temperature"
      , panel = function(x,...){
          panel.abline(h=0, col="black", lwd=0.5)
          panel.qqmathline(x,...)
          panel.qqmath(x,...)
      }
    )
    #print(b)
  })  
  dev.off()

  mainlab <- "Residuals vs. Latitude"

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "resid.vs.lat.lon", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( (tmax - trend - seasonal - spatial) ~ lat | equal.count(lon, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Longitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
           relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Latitude"
      , ylab = "Residual"
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2],evaluation=100,...)
      }
    )
    #print(b)
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "resid.vs.lat.elev", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( (tmax-trend-seasonal-spatial) ~ lat | equal.count(exp(elev2)-128, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Elevation", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Latitude"
      , ylab = "Residual"
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2],evaluation=100,...)
      }
    )
    print(b)
  })
  dev.off()

  mainlab <- "Residuals vs. Longitude"

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "resid.vs.lon.lat", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( (tmax-trend-seasonal-spatial) ~ lon | equal.count(lat, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Latitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Longitude"
      , ylab = "Residual"
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2], evaluation=100,...)
      }
    )
    #print(b)
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "resid.vs.lon.elev", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( tmax-trend-seasonal-spatial ~ lon | equal.count(exp(elev2)-128, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Elevation", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Longitude"
      , ylab = "Residual"
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2],evaluation=100,...)
      }
    )
    #print(b)
  })
  dev.off()
 
  mainlab <- "Residuals vs. Elevation"

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "resid.vs.elev.lat", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( tmax-trend-seasonal-spatial ~ elev2 | equal.count(lat, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Latitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Log (Elevation + 128) (log base 2 meter)"
      , ylab = "Residual"
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2], evaluation=100,...)
      }
    )
    #print(b)
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "resid.vs.elev.lon", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( tmax-trend-seasonal-spatial ~ elev2 | equal.count(lon, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Longitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Log (Elevation + 128) (log base 2 meter)"
      , ylab = "Residual"
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2],evaluation=100,...)
      }
    )
    #print(b)
  })
  dev.off()

}