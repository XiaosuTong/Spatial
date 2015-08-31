source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")
source("~/Projects/Spatial/NCAR/myloess/kdtree.R")

par$N <- 576
par$family <- "symmetric"
par$all <- TRUE
par$degree <- 2
par$span <- seq(0.01,0.1,0.005)
par$drop <- c(TRUE, FALSE) # drop for the qudratic term of elevation

sampleStation <- function(df, seed, rep) {
  
  data <- subset(df, !is.na(resp))

  row.names(data) <- 1:nrow(data)

  rst <- cppkdtree(as.matrix(data[,c("lat", "lon")]), 200)
  
  idx <- ddply(
	  .data = rst,
	  .variable = "leaf",
	  .fun = function(r) {
		  set.seed(seed)
		  data.frame(rowID = sample(r$idx, rep))
	  }
  )

  predStation <- as.character(data$station.id[idx$rowID])

  predStation

}

  
cv <- function(a,b) {
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      
      v <- map.values[[r]]
      v$elev2 <- log2(v$elev + 128)
      seed <- set.seed(map.keys[[r]])
    
      stations <- sampleStation(v, seed, 4)
      if (par$drop[a]) {
        lo.fit <- my.loess2( resp ~ lon + lat + elev2, 
          data       = subset(v, !(station.id %in% stations)), 
          degree     = par$degree, 
          span       = par$span[b],
          parametric = "elev2",
          drop.square = "elev2",
          family     = par$family
        )
      } else {
        lo.fit <- my.loess2( resp ~ lon + lat + elev2, 
          data       = subset(v, !(station.id %in% stations)), 
          degree     = par$degree, 
          span       = par$span[b],
          parametric = "elev2",
          family     = par$family
        )
      }
      
      value <- subset(v, station.id %in% stations)  
      fit <- my.predict.loess(
        object = lo.fit, 
        newdata = data.frame(
          lon = value$lon, 
          lat = value$lat,
          elev2 = value$elev2
        )
      )
      value$fitted <- fit
      value$station.id <- as.character(value$station.id)
      rhcollect(map.keys[[r]], value)
    })
  })
  job$setup <- expression(
    map = {
      library(plyr)
      system("chmod 777 myloess2.so")
      dyn.load("myloess2.so")
      system("chmod 777 cppkdtree.so")
      dyn.load("cppkdtree.so")
    }
  )
  job$shared <- c(
    file.path(rh.root, par$dataset, "shareRLib", "myloess2.so"),
    file.path(rh.root, par$dataset, "shareRLib", "cppkdtree.so")
  )
  job$parameters <- list(
    a = a,
    b = b,
    par = par,
    sampleStation = sampleStation,
    my.loess2 = my.loess2,
    my.simple2 = my.simple2,
    my.predict.loess = my.predict.loess,
    my.predLoess = my.predLoess,
    cppkdtree = cppkdtree
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bymonth"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(
      rh.root, par$dataset, "a1950", "spatial", par$family, "spatial.cv", 
      paste("elevd", as.numeric(par$drop[a]), "sp", par$span[b], sep="")
    ), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 20,  #cdh3,4
    mapreduce.job.reduces = 20  #cdh5 
    #rhipe_reduce_buff_size = 10000
  )
  job$mon.sec <- 10
  job$jobname <- file.path( 
    paste("elevd", as.numeric(par$drop[a]), "sp", par$span[b], sep="")
  )
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)

}


job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
    file <- Sys.getenv("mapred.input.file")
    drop <- substr(tail(strsplit(file, "/")[[1]],3)[2], 6, 6)
    span <- substr(tail(strsplit(file, "/")[[1]],3)[2], 9, 13)
    residual <- sum(with(map.values[[r]], (resp - fitted)^2), na.rm=TRUE)
    rhcollect(c(drop, span), residual)
  })
})
job$reduce <- expression(
  pre= {
    combine <- 0
  },
  reduce = {
    combine <- sum(combine, unlist(reduce.values))
  },
  post = {
    rhcollect(reduce.key, combine/(1024*576))
  }
)
job$input <- rhfmt(
  file.path(rh.root, par$dataset, "a1950", "spatial", par$family, "spatial.cv"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.root, par$dataset, "a1950", "spatial", par$family, "spatial.cvComp"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 10,  #cdh3,4
  mapreduce.job.reduces = 10  #cdh5 
  #rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 10
job$jobname <- file.path("spatial.cvComp")  
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


rst <- rhread("/wsc/tongx/Spatial/tmp/tmax/a1950/spatial/symmetric/spatial.cvComp")
result <- data.frame(
  matrix(as.numeric(unlist(lapply(rst, "[[", 1))), ncol=2, byrow=T)
)
result$MSE <- unlist(lapply(rst, "[[", 2))
names(result)[1:2] <- c("drop","span")

trellis.device(
  device = postscript, 
  file = file.path(local.root, "output", "crossValidMean.ps"),
  color = TRUE, 
  paper = "letter"
)
  xyplot(MSE ~ span
    , data = arrange(result, span)
    , group = drop
    , key = list(
        text = list(label=c("qudratic elevation","linear elevation")),
        lines = list(pch=1, cex=1, type=c("p","p"), col=col[1:2]), 
        columns = 2
      )
    , xlab = list(label="Span", cex=1.5)
    , ylab = list(label="MSE", cex=1.5)
    , scale = list(
        x = list(at=seq(0.01,0.1, by=0.01), cex=1.2),
        y = list(cex=1.2)
      )
    , panel = function(x,y, ...) {
        panel.xyplot(x,y,pch=1, type="b",cex=1,...)
        panel.abline(h=min(y), v=x[which.min(y)], lty=1, lwd=0.5, col="black")
    }
  )
dev.off()


