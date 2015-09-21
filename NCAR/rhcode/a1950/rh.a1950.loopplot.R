##############################################
## check the convergence of findal residuals##
##############################################
Qrst <- function(x, n) {
  x <- x[!is.na(x)]
  a <- sort(x)
  idx <- round(seq(1, length(x), length.out = n))
  f.value <- (idx - 0.5) / length(a)
  qnorm <- qnorm(f.value)
  value <- data.frame(
    residual = a[idx[2:(n-1)]], 
    qnorm = qnorm[2:(n-1)], 
    fv = f.value[2:(n-1)]
  )
}


compare <- function(comp = "resid", family, type, degree, span, index) {
  
  rst <- rhread(file.path(
    rh.root, par$dataset, "a1950", "backfitting", family, type, degree,
    paste("sp", span[1], sep=""), index, paste(comp, "compare", sep="")
  ))

  trellis.device(
    device = postscript, 
    file = file.path(
      local.root, "output", paste(par$dataset, "backfitting", comp,"compare", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "letter"
  )
  for(i in 1:576) {
    tmp <- ddply(.data=rst[[i]][[2]], .vari="iter", .fun = function(r) {Qrst(r$target, n=1000)})

    b <- xyplot(residual ~ residual | iter 
      , data = subset(tmp, as.numeric(iter) <5)
      , xlab = list(label="Residual from i+1th iteration", cex=1.5)
      , ylab = list(label="Residual from ith iteration", cex=1.5)
      , scale = list(cex=1.2)
      , aspect = 1
      , sub = paste(rst[[i]][[1]][1], rst[[i]][[1]][2])
      , panel = function(x,y,subscripts,...) {
          index <- as.numeric(unique(tmp[subscripts, "iter"]))
          sub1 <- subset(tmp, as.numeric(iter) == index)
          sub2 <- subset(tmp, as.numeric(iter) == (index+1))
          panel.xyplot(y=sub1$residual, x=sub2$residual, pch=1)
          panel.abline(a=c(0,1))
      }
    )
    print(b)
  }
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
    key <- substr(unlist(strsplit(tail(strsplit(file, "/")[[1]],3)[2], "[.]")), 8, 9)
    value <- with(map.values[[r]], tmax - spatial - trend - seasonal)
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
    paste("sp", par$span[1], sep=""), paste(par$loess, par$loop, par$type, sep="."), paste("Spatial",1:length(par$span), sep="")
  ), 
  type = "sequence"
)
job5$output <- rhfmt(
  file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
    paste("sp", par$span[1], sep=""), paste(par$loess, par$loop, par$type, sep="."), "residuals.qnorm"
  ), 
  type = "sequence"
)
job5$mapred <- list(mapred.reduce.tasks = 1)
job5$mon.sec <- 10
job5$jobname <- file.path(
  rh.datadir, par$dataset, "spatial", "a1950", par$family, 
  paste("sp", par$span[1], sep=""), paste(par$loess, par$loop, par$type, sep="."), "residuals.qnorm"
)  
job5$readback <- FALSE
job.mr <- do.call("rhwatch", job5)

QQ.resid <- function(type = "same") {
  
  rst <- rhread(  
    file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
      paste("sp", par$span[1], sep=""), paste(par$loess, par$loop, type, sep="."), "residuals.qnorm"
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
    paper = "letter"
  )
    b <- xyplot(residual ~ qnorm | factor(iter)
      , data = result
      , pch = 1
      , cex = 0.4
      , aspect = "xy"
      , layout = c(5,2)
      , scale = list(cex=1.2)
      , xlab = list(label="Unit normal quantile", cex=1.5)
      , ylab = list(label="Residuals", cex=1.5)
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
    file.path(
      rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", par$span[1], sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), paste("Spatial", ".bystation", sep="")
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
    paper = "letter"
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
    b <- qqmath( ~ (value-mean) | factor(comp, levels= c("seasonal", "trend", "spatial", "residual"))
      , data = data
      , distribution = qunif
      , aspect = 2
      , layout = c(4,1)
      , sub = paste("Station from cell", tmax.sample.a1950[tmax.sample.a1950$station.id == rst[[r]][[1]], 2]) 
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Maximum Temperature (degrees centigrade)", cex=1.5)
      , scale = list(cex=1.2)
      , panel = function(x,...){
          panel.abline(h=seq(-15,15,by=3),v=seq(0,1,by=0.2), lwd=0.5, lty=1, col="lightgray")
          panel.abline(h=0, col="black", lwd=0.5)
          panel.qqmath(x,pch=1, cex=0.4, ...)
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
    paper = "letter"
  )
  lapply(1:128, function(i) {
    r <- with(stationLeaf, stationLeaf[leaf == i, "idx"])
    data <- rst[[r]][[2]]
    data <- arrange(data, time)
    data$time1 <- rep(0:95, 6)
    data$factor <- factor(
      rep(1:6, each = 96),
      levels = c(6:1)
    )
    b <- xyplot((tmax - trend - seasonal - spatial) ~ time1 | factor
      , data = data
      , strip = FALSE
      , layout = c(1,6)
      , xlim = c(-1, 96)
      , scale = list(
          x = list(at=seq(0, 95, by=12), cex=1.2), 
          y = list(cex=1.2, tick.number=3)
        )
      , key = list(
          text=list(label=c("remainder","loess smoothing")), 
          lines=list(pch=16, cex=1, lwd=2, type=c("p","l"), col=col[1:2]), 
          columns=2
        )
      , between = list(y=0.5)
      , sub = paste("Station from cell", tmax.sample.a1950[tmax.sample.a1950$station.id == rst[[r]][[1]], 2]) 
      , xlab = list(label="Month", cex=1.5)
      , ylab = list(label="Residuals", cex=1.5)
      , panel = function(x,y,...){
          panel.abline(h=0, col="black", lwd=0.5)
          panel.loess(x,y,degree=2,span=0.15, col=col[2], lwd = 2, evaluation=200,...)
          panel.xyplot(x,y, pch=16, cex=0.6,...)
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
    paper = "letter"
  )
  b <- qqmath(~ (tmax - trend - seasonal - spatial) | factor(leaf)
    , data = result
    , aspect = "xy"
    , layout = c(16,1)
    , xlab = list(label="Unit normal quantile", cex=1.5)
    , ylab = list(label="Residuals", cex=1.5)
    , scale = list(cex=1.2)
    , panel = function(x,...) {
        panel.qqmathline(x,...)
        panel.qqmath(x,pch=16, cex=0.4, ...)
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
    paper = "letter"
  )
  b <- xyplot( seasonal ~ factor(month, levels=c(
    "Jan","Feb","Mar","Apr","May","June",
    "July","Aug", "Sep", "Oct", "Nov", "Dec"
    )) | factor(leaf)
    , data = arrange(result, time)
    , subset = year == "1950"
    , scale = list(x=list(at=c(1, 3, 5, 7, 9, 11), cex=1.2), y=list(cex=1.2))
    , layout = c(4,4)
    , xlab = list(label="Month", cex=1.5)
    , ylab = list(label="Seasonal Component", cex=1.2)
    , panel = function(x,y,...) {
        panel.xyplot(x,y,pch=1, cex=0.5, type="b", ...)
    }
  )
  print(b)
  dev.off()

  mav <- function(x,n=10){  

    data.frame(mv=matrix(filter(x$mean,rep(1/n,n), sides=2), ncol=1), stringsAsFactors=FALSE)  

  }  

  dr <- ddply(
    .data = result,
    .variable = c("station.id", "year"),
    .fun = summarise,
    mean = mean(tmax, na.rm=TRUE)
  ) %>% ddply(
    .variable = "station.id",
    .fun= mav
  )  

  result <- arrange(result, station.id, time)
  result <- dr[rep(row.names(dr), each=12), "mv", drop=FALSE] %>% cbind(result)

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "trend.by.stations", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "letter"
  )
  b <- xyplot( mv ~ time | factor(leaf) 
    , data = result
    , xlab = list(label = "Month", cex=1.5)
    , ylab = list(label = "Maximum Temperature (degrees centigrade)", cex=1.5)
    , xlim = c(0, 576)
    , layout = c(4,3)
    , key=list(
        text = list(label=c("trend component","moving average of yearly mean")), 
        lines = list(pch=16, cex=1, lwd=2, type=c("l","p"), col=col[1:2]),
        columns=2
      )
    , prepanel = function(x,y,subscripts,...){ 
        v <- result[subscripts,] 
        ylim <- range(v$trend, na.rm=TRUE) 
        ans <- prepanel.default.xyplot(v$time, v$mv, ...) 
        ans$ylim <- range(c(ans$ylim, ylim), na.rm=TRUE) 
        ans 
      }
    , scales = list(
        y = list(relation = 'free', cex=1.2), 
        x=list(at=seq(0, 576, by=120), cex=1.2)
      )
    , panel = function(x, y, subscripts, ...) {
        sub <- result[subscripts, ]
        panel.xyplot(
          x = sub$time[seq(1, 576, by=12)], y = sub$mv[seq(1, 576, by=12)],
          type="p", col=col[2], cex = 0.5, pch=16,...
        )
        panel.xyplot(x = sub$time, y=sub$trend, type="l", col=col[1], lwd=2, ...)
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
    paper = "letter"
  )
    for(i in 1:128){
      b <- xyplot( (tmax - trend - seasonal - spatial) ~ as.numeric(year) | factor(month, levels=monthLevel)
        , data = result
        , subset = leaf == i
        , xlab = list(label = "Year", cex=1.5)
        , ylab = list(label = "Maximum Temperature (degrees centigrade)", cex=1.5)
        , sub = paste("Station from cell", i)
        , layout = c(12,1)
        , key=list(
            text = list(label=c("residual","loess smoothing")), 
            lines = list(pch=16, cex=1, lwd=2, type=c("p","l"), col=col[1:2]),
            columns=2
          )
        , scales = list(
            y = list(cex=1.2), 
            x = list(tick.number=4, cex=1.2)
          )
        , panel = function(x,y,...){
            panel.abline(h=0, color="black", lty=1)
            panel.xyplot(x,y, pch=16, cex=0.5,...)
            panel.loess(x,y,span=3/4, degree=1, lwd=2, col=col[2],...)
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
  #  for(i in 1:128){
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
  #  }
  dev.off()

  ACF <- ddply(
    .data = result,
    .variables = "leaf",
    .fun = function(r) {
      corr <- acf(with(r, tmax - trend - seasonal - spatial), plot=FALSE)
      data.frame(correlation = corr$acf, lag = corr$lag)
    }
  )

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "residACF.by.stations", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "letter"
  )     
    b <- xyplot( correlation ~ lag | factor(leaf)
      , data = ACF
      , subset = lag != 0
      , xlab = list(label = "Lag", cex=1.5)
      , ylab = list(label = "ACF", cex=1.5)
      , scale = list(cex=1.2)
      , layout = c(2,2)
      , panel = function(x,y,...) {
          panel.abline(h=0)
          panel.xyplot(x,y, type="h", ...)
        }
    )
    print(b)
  dev.off()

}

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

spatialCheck <- function(type = "same", iter = 5) {
  
  library(plyr)
  library(magrittr)

  rst <- rhread(  
    file.path(
      rh.datadir, par$dataset, "spatial", "a1950", par$family, 
      paste("sp", par$span[1], sep=""), paste(par$loess, par$loop, par$type, sep="."), paste("Spatial", length(par$span), sep="")
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
    paper = "letter"
  )
  lapply(mod, function(r) {
    data <- subset(rst[[r]][[2]], !is.na(tmax)&!is.na(fitted))
    data <- data.frame(cbind(
      do.call(c, with(data, list(trend, seasonal, spatial, tmax-trend-seasonal-spatial))), 
      rep(c("trend", "seasonal", "spatial", "residual"), each=nrow(data))
    ), stringsAsFactors = FALSE)
    names(data) <- c("value", "comp")
    data$value <- as.numeric(data$value)
    tmp <- ddply(
      .data = data,
      .variable = "comp",
      .fun = summarise,
      m = mean(value, na.rm = TRUE)
    )
    data$mean <- rep(tmp[c(4,2,3,1),2], each = nrow(data)/4)
    b <- qqmath( ~ (value-mean) | factor(comp, levels= c("seasonal", "trend", "spatial", "residual"))
      , data = data
      , distribution = qunif
      , aspect = 2
      , layout = c(4,1)
      , scale = list(cex=1.2)
      , sub = paste(rst[[r]][[1]][1], rst[[r]][[1]][2])
      , xlab = list(label="f-value", cex=1.5)
      , ylab = list(label="Maximum Temperature (degrees centigrade)", cex=1.5)
      , panel = function(x,...){
          panel.abline(h=seq(-15,15,by=5),v=seq(0,1,by=0.2), lwd=0.5, lty=1, col="lightgray")
          panel.abline(h=0, col="black", lwd=0.5)
          panel.qqmath(x,pch=1, cex=0.4, ...)
      }
    )
  })  
  dev.off()

  
  data <- ldply(1:576, function(r) {
    Qrst(with(rst[[r]][[2]], tmax-trend-seasonal-spatial))
  })
  key <- data.frame(
    matrix(unlist(lapply(rst, "[[", 1)), ncol=2, byrow=TRUE), stringsAsFactors=FALSE
  )
  names(key) <- c("year","month")
  data <- cbind(
    data, key[rep(1:nrow(key), each=98),]
  )

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "QQ.resid.by.month", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "letter"
  )
    b <- xyplot( residual ~ qnorm | factor(month)*factor(year)
      , data = data
      , pch = 1
      , cex = 0.4
      , aspect = 1
      , layout = c(4,3)
      , scale = list(cex=1.2)
      , xlab = list(label="Unit normal quantile", cex=1.5)
      , ylab = list("Maximum Temperature (degrees centigrade)", cex=1.5)
      , panel = function(x,y,...){
          b1 <- quantile(y, probs = 0.25)
          a1 <- qnorm(0.25)
          b2 <- quantile(y, probs = 0.75)
          a2 <- qnorm(0.75)
          panel.abline(b = (b2-b1)/(a2-a1), a = (b1 - a1*(b2-b1)/(a2-a1)), col="black", lwd=0.5)
          panel.xyplot(x,y,...)
      }
    )
    print(b)
  dev.off()


  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "resid.vs.lat.lon", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "letter"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( (tmax - trend - seasonal - spatial) ~ lat | equal.count(lon, 20, overlap=0)
      , data = rst[[r]][[2]]
      , strip=strip.custom(var.name = "Longitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(relation = "free", cex=1.2),
          x = list(relation = "free", tick.number = 3, cex=1.2)
        )
      , layout = c(5,4)
      , xlab = list(label="Latitude", cex=1.5)
      , ylab = list(label="Residual", cex=1.5)
      , sub = paste(rst[[r]][[1]][1], rst[[r]][[1]][2])
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
    paper = "letter"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( (tmax-trend-seasonal-spatial) ~ lat | equal.count(exp(elev2)-128, 20, overlap=0)
      , data = rst[[r]][[2]]
      , strip=strip.custom(var.name = "Elevation", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(relation = "free", cex=1.2),
          x = list(relation = "free", tick.number = 3, cex=1.2)
        )
      , layout = c(5,4)
      , xlab = list(label="Latitude", cex=1.5)
      , ylab = list(label="Residual", cex=1.5)
      , sub = paste(rst[[r]][[1]][1], rst[[r]][[1]][2])
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
      local.output, paste(par$dataset, "backfitting", type, "resid.vs.lon.lat", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "letter"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( (tmax-trend-seasonal-spatial) ~ lon | equal.count(lat, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Latitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(relation = "free", cex=1.2),
          x = list(relation = "free", tick.number = 3, cex=1.2)
        )
      , layout = c(5,4)
      , xlab = list(label="Longitude", cex=1.5)
      , ylab = list(label="Residual", cex=1.5)
      , sub = paste(data[[r]][[1]][1], data[[r]][[1]][2])
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
    paper = "letter"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( tmax-trend-seasonal-spatial ~ lon | equal.count(exp(elev2)-128, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Elevation", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(relation = "free", cex=1.2),
          x = list(relation = "free", tick.number = 3, cex=1.2)
        )
      , layout = c(5,4)
      , xlab = list(label="Longitude", cex=1.5)
      , ylab = list(label="Residual", cex=1.5)
      , sub = paste(data[[r]][[1]][1], data[[r]][[1]][2])
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
      local.output, paste(par$dataset, "backfitting", type, "resid.vs.elev.lat", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "letter"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( tmax-trend-seasonal-spatial ~ elev2 | equal.count(lat, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Latitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(relation = "free", cex=1.2),
          x = list(relation = "free", tick.number = 3, cex=1.2)
        )
      , layout = c(5,4)
      , xlab = list(label="Log (Elevation + 128) (log base 2 meter)", cex=1.5)
      , ylab = list(label="Residual", cex=1.5)
      , sub = paste(data[[r]][[1]][1], data[[r]][[1]][2])
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
    paper = "letter"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( tmax-trend-seasonal-spatial ~ elev2 | equal.count(lon, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Longitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(relation = "free", cex=1.2),
          x = list(relation = "free", tick.number = 3, cex=1.2)
        )
      , layout = c(5,4)
      , xlab = list(label="Log (Elevation + 128) (log base 2 meter)", cex=1.5)
      , ylab = list(label="Residual", cex=1.5)
      , sub = paste(data[[r]][[1]][1], data[[r]][[1]][2])
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