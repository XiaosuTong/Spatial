##################
##
##  loopnew is the current version of backfitting, inner=5, outer=5
##
##  loopout1 is inner=40, outer=1
##
##  loop.periodic is inner= 40, outer=1, sw="periodic", tw=109, td=2, elev.degree=2
##
##  looptest is change the spatial loess to be tmax-trend ~ lon + lat + elev WRONG!!!!!backfitting is substract out all other
##  component when fit one
##  sw = "periodic", tw = 109, td = 2
##
##  loop.nonperiodic is inner = 40, outer=1, sw = 37, sd =1, tw = 109, td=1
##
##  loop.longtrend is inner = 20, outer=1, sw = 37, sd = 1, td=1, tw=425
##
##################
par$modified <- TRUE
par$loess <- "loess04" # loess04 is w/ elevation
par$family <- "symmetric"
par$degree <- 2
par$outer <- 1
par$loop <- "loop.longtrend" 
par$type <- "same" # or "same", "decr"
par$parameters <- list(
	sw = 37,
	tw = 425,
	sd = 1,
	td = 1,
	inner = 1,
	outer = 1
) 
if (par$type == "same") {
  par$span <- rep(0.035, 40)
} else if (par$type == "incr") {
  par$span <- c(seq(0.03, 0.05, by=0.005))
} else {
  par$span <- c(seq(0.05, 0.03, by=-0.005))
}
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")

# The first input file is the loess04.bystation.all which imputed all NA
# Key is the station.id

backfitSwap.station <- function() {

    # job3 is just changing the key from month and year to station.id from job22
    job3 <- list()
    job3$map <- expression({
      lapply(seq_along(map.values), function(r) {
        map.values[[r]]$year <- map.keys[[r]][1]
        map.values[[r]]$month <- map.keys[[r]][2]
        lapply(1:dim(map.values[[r]])[1], function(i){
          value <- map.values[[r]][i, ]
          rhcollect(as.character(value$station.id), value)
        })
      })
    })
    job3$reduce <- expression(
      pre = {
        combine <- data.frame()
      },
      reduce = {
        combine <- rbind(combine, do.call(rbind, reduce.values))
      },
      post = {
        rhcollect(reduce.key, combine)
      }
    )
    job3$combiner <- TRUE
    job3$input <- rhfmt(FileInput , type = "sequence")
    job3$output <- rhfmt(FileOutput, type = "sequence")
    job3$mapred <- list(mapred.reduce.tasks = 72)
    job3$mon.sec <- 10
    job3$jobname <- FileOutput  
    job3$readback <- FALSE  

    job.mr <- do.call("rhwatch", job3)

}

backfitSwap.month <- function() {

  job21 <- list()
  job21$map <- expression({
    lapply(seq_along(map.values), function(r) {
      d_ply(
        .data = map.values[[r]],
        .variable = c("year","month"),
        .fun = function(k, station = map.keys[[r]]) {
          key <- c(unique(k$year), unique(as.character(k$month)))
          value <- subset(k, select = -c(year, month))
          value$station.id <- station
          value$lon <- as.numeric(attributes(map.values[[r]])$location$lon)
          value$lat <- as.numeric(attributes(map.values[[r]])$location$lat)
          value$elev2 <- as.numeric(attributes(map.values[[r]])$location$elev2)
          rhcollect(key, value)
      })
    })
  })
  job21$reduce <- expression(
    pre = {
      combine <- data.frame()
    },
    reduce = {
      combine <- rbind(combine, do.call(rbind, reduce.values))
    },
    post = {
      rhcollect(reduce.key, combine)
    }
  )
  job21$setup <- expression(
    map = {
      library(plyr)
    }
  )
  job21$mapred <- list(
    mapred.reduce.tasks = 72,
    mapred.tasktimeout = 0,
    rhipe_reduce_buff_size = 10000
  )
  job21$combiner <- TRUE
  job21$input <- rhfmt(FileInput, type="sequence")
  job21$output <- rhfmt(FileOutput, type="sequence")
  job21$mon.sec <- 10
  job21$jobname <- FileOutput
  job21$readback <- FALSE  

  job.mr <- do.call("rhwatch", job21)  

}

backfitSpatial <- function() {

  job22 <- list()
  job22$map <- expression({
    lapply(seq_along(map.values), function(r) {
      value <- map.values[[r]]
      lo.fit <- my.loess2( (resp - trend - seasonal) ~ lon + lat + elev2, 
        data    = subset(value, !is.na(fitted)), 
        degree  = argumt$degree, 
        span    = argumt$span,
        weights = subset(value, !is.na(fitted))$weights,
        family  = argumt$family,
        parametric = "elev2",
        control = loess.control(surface = "direct")
      )
      fit <- my.predict.loess(
        object = lo.fit, 
        newdata = data.frame(
          lon = value$lon, 
          lat = value$lat,
          elev2 = value$elev2
        )
      )
      value$spatial <- fit
      rhcollect(map.keys[[r]], value)
    })
  })
  job22$parameters <- list(
    argumt = arg,
    my.loess2 = my.loess2,
    my.simple2 = my.simple2,
    my.predict.loess = my.predict.loess,
    my.predLoess = my.predLoess
  )
  job22$setup <- expression(
    map = {
      library(maps)
      system("chmod 777 myloess2.so")
      dyn.load("myloess2.so")
    }
  )
  job22$shared <- c(
    file.path(file.path(rh.root, par$dataset, "shareRLib", "myloess2.so"))
  )
  job22$mapred <- list(
    mapred.reduce.tasks = 72,
    mapred.tasktimeout = 0
  )
  job22$input <- rhfmt(FileInput, type="sequence")
  job22$output <- rhfmt(FileOutput, type="sequence")
  job22$mon.sec <- 10
  job22$jobname <- FileOutput
  job22$readback <- FALSE  

  job.mr <- do.call("rhwatch", job22)  

}

backfitSTL <- function() {

  # In each loop, the first job is STL+ at each station, output key is station.id
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      v <- arrange(map.values[[r]], time)
      if (first) {
        Index <- which(is.na(v$resp))
        v$resp[Index] <- v$fitted[Index]
        v$spatial <- 0
        v$trend <- 0
        v$weight <- 1
      } 
      v.stl <- stl3(
        x = with(v, resp - spatial), 
        t = v$time, 
        trend = v$trend,
        weight = v$weight, 
        n.p = 12, 
        s.window = parameters$sw, 
        s.degree = parameters$sd, 
        t.window = parameters$tw, 
        t.degree = parameters$td, 
        inner = parameters$inner, 
        outer = parameters$outer
      )$data
      if (first) {
        value <- cbind(
          subset(v, select = -c(station.id, elev, lon, lat, elev2, trend)), 
          subset(v.stl, select = -c(weights, sub.labels, raw, remainder))
        )
      } else {
        value <- cbind(
          subset(v, select = -c(station.id, lon, lat, elev2, trend, seasonal)), 
          subset(v.stl, select = -c(weights, sub.labels, raw, remainder))
        )        
      }
      attr(value, "location") <- c(v[1, c("lon", "lat", "elev2", "station.id")])
      rhcollect(map.keys[[r]], value)
    })
  })
  job$parameters <- list(
    parameters = par$parameters,
    dataset = par$dataset,
    first = i == 1 && o == 1
  )
  job$setup <- expression(
    map = {
      library(lattice)
      library(yaImpute)
      library(plyr)
      library(stl3)
    }
  )
  job$input <- rhfmt(FileInput , type = "sequence")
  job$output <- rhfmt(FileOutput, type = "sequence")
  job$mapred <- list(mapred.reduce.tasks = 72)
  job$mon.sec <- 5
  job$jobname <- FileOutput  
  job$readback <- FALSE  

  job.mr <- do.call("rhwatch", job)

}

backfitWeights <- function() {
  
    FileOutput <- file.path(
      rh.root, par$dataset, "a1950", "spatial", par$family, paste("sp", par$span[1], sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), "residual"
    )  

    job4 <- list()
    job4$map <- expression({
      lapply(seq_along(map.values), function(r) {
        R.abs <- abs(with(map.values[[r]], resp - spatial - seasonal - trend))
        rhcollect(1, R.abs)
      })
    })
    job4$reduce <- expression(
      pre = {
        combine <- vector()
      },
      reduce = {
        combine <- c(combine, unlist(reduce.values))
      },
      post = {
        n <- length(combine)
        mid1 <- floor(n/2+1)
        mid2 <- n - mid1+1
        h <- 3 * sum(sort(combine)[mid1:mid2])
        h9 <- .999 * h
        h1 <- .001 * h
        rhcollect(reduce.key, c(h, h9, h1))
      }
    )
    job4$input <- rhfmt(FileInput, type = "sequence")
    job4$output <- rhfmt(FileOutput, type = "sequence")
    job4$mapred <- list(mapred.reduce.tasks = 1, rhipe_reduce_buff_size = 10000)
    job4$mon.sec <- 10
    job4$jobname <- FileOutput  
    job4$readback <- TRUE 

  weight <- do.call("rhwatch", job4)[[1]][[2]]
  return(weight)

}

backfitRobust <- function(weight) {

    job5 <- list()
    job5$map <- expression({
      lapply(seq_along(map.values), function(r) {
        R.abs <- abs(with(map.values[[r]], resp - spatial - seasonal - trend))
        w <- (1 - (R.abs / h)^2)^2
        w[R.abs <= h1] <- 1
        w[R.abs >= h9] <- 0
        w[w == 0] <- 1e-6
        w[is.na(w)] <- 1
        map.values[[r]]$weight <- w
        rhcollect(map.keys[[r]], map.values[[r]])
      })
    })
    job5$parameters <- list(
      h  = weight[1],
      h9 = weight[2],
      h1 = weight[3]
    )
    job5$input <- rhfmt(FileInput, type = "sequence")
    job5$output <- rhfmt(FileOutput, type = "sequence")
    job5$mapred <- list(mapred.reduce.tasks = 72)
    job5$mon.sec <- 10
    job5$jobname <- FileOutput  
    job5$readback <- FALSE 

    job.mr <- do.call("rhwatch", job5)

}

backfitAll <- function(span, family, type, parameter) {

  FileInput <- file.path(
    rh.root, par$dataset, "a1950", "bymonth.fit", family, type, degree, paste("sp", span[1], sep="")
  )


for(o in 1:par$outer) {
  
  for(i in 1:length(par$span)) {
    
    FileOutput <- file.path(
      rh.root, par$dataset, "a1950", "spatial", par$family, paste("sp", par$span[1], sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), "STL"
    )  

    backfitSTL()
    
    # The output from job is the input to job2
    FileInput <- FileOutput
    
    # Second job output is the spatial fitting, output key is month and year
    FileOutput <- file.path(
      rh.root, par$dataset, "a1950", "spatial", par$family, paste("sp", par$span[1], sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), "byMonth"
    )  

    backfitSwap.month()
 
    # The output from job21 is the input to job22
    FileInput <- FileOutput  

    FileOutput <- file.path(
      rh.root, par$dataset, "a1950", "spatial", par$family, paste("sp", par$span[1], sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), paste("Spatial", i, sep="")
    )  

    arg <- list(degree = par$degree, span = par$span[i], family = "gaussian")



    # The output from job22 is the input to job3
    FileInput <- FileOutput  

    FileOutput <- file.path(
      rh.root, par$dataset, "a1950", "spatial", par$family, paste("sp", par$span[1], sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), paste("Spatial", ".bystation", sep="")
    )    

    backfitSwap.station()

    FileInput <- FileOutput  

  }

  #outer that calculate the weights
  
  if(par$outer > 1) {
    
    w <- backfitWeights()

    FileOutput <- file.path(
      rh.root, par$dataset, "a1950", "spatial", par$family, paste("sp", par$span[1], sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), paste("Outer", o, ".bystation", sep="")
    )
    
    backfitRobust(weight = w)
    
    FileInput <- FileOutput

  }

}

}