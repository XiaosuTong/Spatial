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
    key <- substr(unlist(strsplit(tail(strsplit(file, "/")[[1]],3)[2], "[.]")), 8, 9)
    value <- with(map.values[[r]], fitted-spatial)
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
    paste("sp", "0.05", sep=""), paste(par$loess, "loop", par$type, sep="."), paste("Spatial",1:length(par$span), sep="")
  ), 
  type = "sequence"
)
job5$output <- rhfmt(
  file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
    paste("sp", "0.05", sep=""), paste(par$loess, "loop", par$type, sep="."), "residuals.qnorm"
  ), 
  type = "sequence"
)
job5$mapred <- list(mapred.reduce.tasks = 1)
job5$mon.sec <- 10
job5$jobname <- file.path(
  rh.datadir, par$dataset, "spatial", "a1950", par$family, 
  paste("sp", "0.05", sep=""), paste(par$loess, "loop", par$type, sep="."), "residuals.qnorm"
)  
job5$readback <- FALSE
job.mr <- do.call("rhwatch", job5)

QQ.resid <- function(type = "same") {
  
  rst <- rhread(  
    file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
      paste("sp", "0.05", sep=""), paste(par$loess, "loop", type, sep="."), "residuals.qnorm"
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
      , layout = c(5,2)
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
tempCheck <- function(type = "same", iter = 20) {

  library(plyr)
  library(magrittr)

  load(file.path(local.datadir, "samplestation.a1950.RData")) ##kd-tree is built in kdfindcells.R 
  
  rst <- rhread(  
    file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
      paste("sp", "0.05", sep=""), paste(par$loess, "loop", type, sep="."), paste("Spatial", iter, ".bystation", sep="")
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
    b <- xyplot((fitted - spatial) ~ time1 | factor
      , data = data
      , pch = 1
      , strip = FALSE
      , cex = 0.4
      , aspect = "xy"
      , layout = c(1,6)
      , xlim = c(-1, 96)
      , scale = list(x = list(at=seq(0, 95, by=12)))
      , main = paste("Residuals from iteration", iter)
      , sub = paste("Station from cell", tmax.sample.a1950[tmax.sample.a1950$station.id == rst[[r]][[1]], 2]) 
      , xlab = "Month"
      , ylab = "Residuals"
      , panel = function(x,y,...){
          panel.abline(h=0, col="black", lwd=0.5)
          panel.xyplot(x,y,...)
      }
    )
    print(b)
  })
  dev.off()

  result <- do.call(rbind, lapply(stationLeaf$idx, function(r){
    subset(rst[[r]][[2]], select = c(fitted, spatial, station.id, time))
  })) %>% merge(stationLeaf, by = "station.id")
  
  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "backfitting", type, "QQ.by.stations", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  b <- qqmath(~ (fitted - spatial) | factor(leaf)
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

}

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

spatialCheck <- function(type = "same", iter = 20) {
  
  library(plyr)
  library(magrittr)

  rst <- rhread(  
    file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, 
      paste("sp", "0.05", sep=""), paste(par$loess, "loop", type, sep="."), paste("Spatial", iter, sep="")
    )
  )  
  
  tmp <- do.call(rbind, lapply(rst,"[[",1)) %>%
    data.frame(stringsAsFactors=FALSE) 
  tmp$X2 <- tmp$X2 %>% substr(1,3) %>% match(month.abb)
  mod <- tmp %>% with(tmp[order(X1, X2),]) %>% row.names() %>% as.numeric()

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
    b <- xyplot( (fitted - spatial) ~ lat | equal.count(lon, 20, overlap=0)
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
    print(b)
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
    b <- xyplot( (fitted-spatial) ~ lat | equal.count(exp(elev2)-128, 20, overlap=0)
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
    b <- xyplot( (fitted-spatial) ~ lon | equal.count(lat, 20, overlap=0)
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
    print(b)
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
    b <- xyplot( fitted-spatial ~ lon | equal.count(exp(elev2)-128, 20, overlap=0)
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
    print(b)
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
    b <- xyplot( fitted-spatial ~ elev2 | equal.count(lat, 20, overlap=0)
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
    print(b)
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
    b <- xyplot( fitted-spatial ~ elev2 | equal.count(lon, 20, overlap=0)
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
    print(b)
  })
  dev.off()

}