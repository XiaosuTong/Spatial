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
parameter <- data.frame(
	sw = "periodic",
	tw = 425,
	sd = 1,
	td = 1,
	fcw = 425,
	fcd = 1,
	scw = 214,
	scd = 2,
  fc.flag = TRUE,
  stringsAsFactors=FALSE
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

backfitSwap.station <- function(input, output) {

    # job3 is just changing the key from month and year to station.id from job22
    job3 <- list()
    job3$map <- expression({
      lapply(seq_along(map.values), function(r) {
        map.values[[r]]$year <- map.keys[[r]][1]
        map.values[[r]]$month <- map.keys[[r]][2]
        lapply(1:nrow(map.values[[r]]), function(i){
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
    job3$input <- rhfmt(input , type = "sequence")
    job3$output <- rhfmt(output, type = "sequence")
    job3$mapred <- list(mapred.reduce.tasks = 72)
    job3$mon.sec <- 10
    job3$jobname <- output  
    job3$readback <- FALSE  

    job.mr <- do.call("rhwatch", job3)

}

backfitSwap.month <- function(input, output) {

  job <- list()
  job$map <- expression({
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
  job$reduce <- expression(
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
  job$setup <- expression(
    map = {
      library(plyr)
    }
  )
  job$mapred <- list(
    mapred.reduce.tasks = 72,
    mapred.tasktimeout = 0,
    rhipe_reduce_buff_size = 10000
  )
  job$combiner <- TRUE
  job$input <- rhfmt(input, type="sequence")
  job$output <- rhfmt(output, type="sequence")
  job$mon.sec <- 10
  job$jobname <- output
  job$readback <- FALSE  

  job.mr <- do.call("rhwatch", job)  

}

backfitSpatial <- function(input, output, argumt, fc.flag) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      value <- map.values[[r]]
      if(!fc.flag) {
        lo.fit <- my.loess2( (resp - trend - seasonal) ~ lon + lat + elev2, 
          data    = subset(value, !is.na(fitted)), 
          degree  = argumt$degree, 
          span    = argumt$span,
          weights = subset(value, !is.na(fitted))$weights,
          family  = argumt$family,
          parametric = "elev2"
        )
      } else {
        lo.fit <- my.loess2( (resp - data.seasonal - fc.first - fc.second) ~ lon + lat + elev2, 
          data    = subset(value, !is.na(fitted)), 
          degree  = argumt$degree, 
          span    = argumt$span,
          weights = subset(value, !is.na(fitted))$weights,
          family  = argumt$family,
          parametric = "elev2"
        )  
      }
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
  job$parameters <- list(
    argumt = arg,
    fc.flag = fc.flag,
    my.loess2 = my.loess2,
    my.simple2 = my.simple2,
    my.predict.loess = my.predict.loess,
    my.predLoess = my.predLoess
  )
  job$setup <- expression(
    map = {
      system("chmod 777 myloess2.so")
      dyn.load("myloess2.so")
    }
  )
  job$shared <- c(
    file.path(file.path(rh.root, par$dataset, "shareRLib", "myloess2.so"))
  )
  job$mapred <- list(
    mapred.reduce.tasks = 72,
    mapred.tasktimeout = 0
  )
  job$input <- rhfmt(input, type="sequence")
  job$output <- rhfmt(output, type="sequence")
  job$mon.sec <- 10
  job$jobname <- output
  job$readback <- FALSE  

  job.mr <- do.call("rhwatch", job)  

}

backfitSTL <- function(input, output, parameter, iiter, oiter) {

  # In each loop, the first job is STL+ at each station, output key is station.id
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      v <- map.values[[r]]
      v <- arrange(v, year, match(month, month.abb))
      if (first) {
        v$time <- 1:576
        Index <- which(is.na(v$resp))
        v$resp[Index] <- v$fitted[Index]
        v$spatial <- 0
        v$trend <- 0
        v$fc.first <- 0
        v$fc.second <- 0
        v$weight <- 1
      } 
      if (!fc.flag) {
        v.stl <- stl3(
          x = with(v, resp - spatial), t = v$time, 
          trend = v$trend, weight = v$weight, n.p = 12, 
          s.window = parameter$sw, s.degree = parameter$sd, t.window = parameter$tw, t.degree = parameter$td, inner = 1, outer = 1
        )$data
        v.stl <- subset(v.stl, select = -c(weights, sub.labels, raw, remainder))
        v <- subset(v, select = -c(trend))
      } else {
        v.stl <- do.call("cbind", stl3(
          x = with(v, resp - spatial), t = v$time, 
          fc.first = v$fc.first, fc.second = v$fc.second, weight = v$weight, n.p = 12, 
          s.window = parameter$sw, s.degree = parameter$sd, t.window = parameter$tw, t.degree = parameter$td, 
          fc.window = c(parameter$fcw, parameter$scw), fc.degree = c(parameter$fcd, parameter$scd), inner = 1, outer = 1
        )[c("data","fc")])
        names(v.stl)[grep("fc.fc", names(v.stl))] <- c("fc.first", "fc.second")
        v.stl <- subset(v.stl, select = -c(data.weights, data.trend, data.sub.labels, data.raw, data.remainder, fc.remainder))
        v <- subset(v, select = -c(fc.first, fc.second))
      } 

      if (!first) {
        if(!fc.flag) {
          v <- subset(v, select = -c(seasonal))
        } else {
          v <- subset(v, select = -c(data.seasonal))
        }
      } else {
        v <- subset(v, select = -c(elev))  
      }
      value <- cbind(v, v.stl)
      attr(value, "location") <- c(v[1, c("lon", "lat", "elev2", "station.id")])
      value <- subset(value, select = -c(station.id, lon, lat, elev2))
      rhcollect(map.keys[[r]], value)
    })
  })
  job$parameters <- list(
    sw = parameter$sw, tw = parameter$tw, sd = parameter$sd, td = parameter$td, fcw = parameter$fcw, 
    fcd = parameter$fcd, scw = parameter$scw, scd = parameter$scd,
    fc.flag = parameter$fc.flag, first = iiter == 1 && oiter == 1
  )
  job$setup <- expression(
    map = {
      library(lattice)
      library(yaImpute)
      library(plyr)
      library(stl3)
    }
  )
  job$input <- rhfmt(input , type = "sequence")
  job$output <- rhfmt(output, type = "sequence")
  job$mapred <- list(mapred.reduce.tasks = 72)
  job$mon.sec <- 10
  job$jobname <- output  
  job$readback <- FALSE  

  job.mr <- do.call("rhwatch", job)

}

backfitWeights <- function(input, output, fc.flag) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      if(!fc.flag) {
        R.abs <- abs(with(map.values[[r]], resp - spatial - seasonal - trend))
      } else {
        R.abs <- abs(with(map.values[[r]], resp - spatial - data.seasonal - fc.first - fc.second))
      }
      rhcollect(1, R.abs)
    })
  })
  job$reduce <- expression(
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
  job$parameters <- list(fc.flag = fc.flag)
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(output, type = "sequence")
  job$mapred <- list(mapred.reduce.tasks = 1, rhipe_reduce_buff_size = 10000)
  job$mon.sec <- 10
  job$jobname <- output  
  job$readback <- TRUE 

  weight <- do.call("rhwatch", job)[[1]][[2]]
  return(weight)

}

backfitRobust <- function(input, output, weight, fc.flag) {

    job <- list()
    job$map <- expression({
      lapply(seq_along(map.values), function(r) {
        if(!fc.flag) {
          R.abs <- abs(with(map.values[[r]], resp - spatial - seasonal - trend))
        } else {
          R.abs <- abs(with(map.values[[r]], resp - spatial - data.seasonal - fc.first - fc.second))
        }
        w <- (1 - (R.abs / h)^2)^2
        w[R.abs <= h1] <- 1
        w[R.abs >= h9] <- 0
        w[w == 0] <- 1e-6
        w[is.na(w)] <- 1
        map.values[[r]]$weight <- w
        rhcollect(map.keys[[r]], map.values[[r]])
      })
    })
    job$parameters <- list(
      h  = weight[1],
      h9 = weight[2],
      h1 = weight[3],
      fc.flag = fc.flag
    )
    job$input <- rhfmt(input, type = "sequence")
    job$output <- rhfmt(output, type = "sequence")
    job$mapred <- list(mapred.reduce.tasks = 72)
    job$mon.sec <- 10
    job$jobname <- output  
    job$readback <- FALSE 

    job.mr <- do.call("rhwatch", job)

}

backfitAll <- function(span, family, type, parameter, index, degree) {

  FileInput <- file.path(
    rh.root, par$dataset, "a1950", "bymonth.fit", family, type, degree, paste("sp", span[1], sep="")
  )

  FileOutput <- file.path(
    rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), "startbymonth"
  )
  
  try(backfitSwap.station(input=FileInput, output=FileOutput))

  FileInput <- FileOutput

  for(o in 1:1) {
    for(i in 1:3) {
      
      FileOutput <- file.path(
        rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), index, "STL"
      )  
  
      try(backfitSTL(input=FileInput, output=FileOutput, parameter=parameter, iiter=i, oiter=o))
      
      # The output from job is the input to job2
      FileInput <- FileOutput
      
      # Second job output is the spatial fitting, output key is month and year
      FileOutput <- file.path(
        rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), index, "bymonth" 
      )  
  
      try(backfitSwap.month(input=FileInput, output=FileOutput))
   
      # The output from job21 is the input to job22
      FileInput <- FileOutput  
  
      FileOutput <- file.path(
        rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), index, paste("spatial", i, sep="")
      )  
  
      arg <- list(degree = 2, span = span[i], family = "symmetric")
  
      try(backfitSpatial(input=FileInput, output=FileOutput, argumt=arg, fc.flag=TRUE))
  
      # The output from job22 is the input to job3
      FileInput <- FileOutput  
  
      FileOutput <- file.path(
        rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), index, "spatial.bystation"
      )    
  
      try(backfitSwap.station(input=FileInput, output=FileOutput))
  
      FileInput <- FileOutput  
  
    }
  
    #outer that calculate the weights
    
    if(par$outer > 0) {
      
      FileOutput <- file.path(
        rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), "residual"
      )  
      try(w <- backfitWeights(input = FileInput, output = FileOutput, fc.flag=TRUE))
  
      FileOutput <- file.path(
        rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), 
        paste("Outer", o, ".bystation", sep="")
      )
      
      try(backfitRobust(input = FileInput, output = FileOutput, weight = w, fc.flag=TRUE))
      
      FileInput <- FileOutput
  
    }
  
  }
}