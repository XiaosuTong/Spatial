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
par$loop <- "loopelevno2" # loopnew is the current version of backfitting, looptest is degree=1 for elevation(drop.square)
par$type <- "same" # or "same", "decr"
par$parameters <- list(
	sw = 37,
	tw = 109,
	sd = 2,
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

# The first input file is the loess04.bystation.all which imputed all NA
# Key is the station.id

FileInput <- file.path(
  rh.datadir, par$dataset, "spatial", "a1950", par$family, 
  paste("sp", "0.05", sep=""), paste(par$loess, "bystation.all", sep=".")
)  


for(o in 1:par$outer) {
  
  for(i in 1:length(par$span)) {
    
    FileOutput <- file.path(
      rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", "0.05", sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), paste("STL", i, sep="")
    )  

    # In each loop, the first job is STL+ at each station, output key is station.id
    job1 <- list()
    job1$map <- expression({
      lapply(seq_along(map.values), function(r) {
        map.values[[r]] <- map.values[[r]][order(map.values[[r]]$time), ]
        if (first) {
          Index <- which(is.na(map.values[[r]]$tmax))
          map.values[[r]]$tmax[Index] <- map.values[[r]]$fitted[Index]
          map.values[[r]]$spatial <- 0
          map.values[[r]]$trend <- 0
          map.values[[r]]$weight <- 1
        } 
        v.stl <- stl3(
          x = with(map.values[[r]], tmax - spatial), 
          t = map.values[[r]]$time, 
          trend = map.values[[r]]$trend,
          weight = map.values[[r]]$weight, 
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
            subset(map.values[[r]], select = -c(station.id, elev, lon, lat, elev2, trend)), 
            subset(v.stl, select = -c(weights, sub.labels, raw, remainder))
          )
        } else {
          value <- cbind(
            subset(map.values[[r]], select = -c(station.id, lon, lat, elev2, trend, seasonal)), 
            subset(v.stl, select = -c(weights, sub.labels, raw, remainder))
          )        
        }
        attr(value, "location") <- c(map.values[[r]][1, c("lon", "lat", "elev2", "station.id")])
        rhcollect(map.keys[[r]], value)
      })
    })
    job1$parameters <- list(
      parameters = par$parameters,
      dataset = par$dataset,
      first = i == 1 && o == 1
    )
    job1$setup <- expression(
      map = {
        library(lattice)
        library(yaImpute, lib.loc = lib.loc)
        library(stl3, lib.loc = lib.loc)
      }
    )
    job1$input <- rhfmt(FileInput , type = "sequence")
    job1$output <- rhfmt(FileOutput, type = "sequence")
    job1$mapred <- list(mapred.reduce.tasks = 72)
    job1$mon.sec <- 5
    job1$jobname <- FileOutput  
    job1$readback <- FALSE  

    job.mr <- do.call("rhwatch", job1)
    
    arg <- list(degree = par$degree, span = par$span[i], family = "gaussian")
    
    # The output from job1 is the input to job2
    FileInput <- FileOutput
    
    # Second job output is the spatial fitting, output key is month and year
    FileOutput <- file.path(
      rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", "0.05", sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), paste("Spatial", i, sep="")
    )  

    job2 <- list()
    job2$map <- expression({
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
    job2$reduce <- expression(
      pre = {
        combine <- data.frame()
      },
      reduce = {
        combine <- rbind(combine, do.call(rbind, reduce.values))
      },
      post = {
        lo.fit <- my.loess2( (tmax - seasonal - trend) ~ lon + lat + elev2, 
          data    = subset(combine, !is.na(fitted)), 
          degree  = argumt$degree, 
          span    = argumt$span,
          weights = subset(combine, !is.na(fitted))$weights,
          family  = argumt$family,
          parametric = "elev2",
          drop.square = "elev2",
          control = loess.control(surface = "direct")
        )
        fit <- my.predict.loess(
          object = lo.fit, 
          newdata = data.frame(
            lon = combine$lon, 
            lat = combine$lat,
            elev2 = combine$elev2
          )
        )
        combine$spatial <- fit 
        rhcollect(reduce.key, combine)
      }
    )
    job2$parameters <- list(
      argumt = arg,
      my.loess2 = my.loess2,
      my.simple2 = my.simple2,
      my.predict.loess = my.predict.loess,
      my.predLoess = my.predLoess
    )
    job2$setup <- expression(
      map = {
        library(plyr)
      },
      reduce = {
        library(maps, lib.loc = lib.loc)
        system("chmod 777 myloess2.so")
        dyn.load("myloess2.so")
      }
    )
    job2$shared <- c(
      file.path(rh.datadir, par$dataset, "shareRLib", "myloess2.so")
    )
    job2$mapred <- list(
      mapred.reduce.tasks = 72,
      mapred.tasktimeout = 0,
      rhipe_reduce_buff_size = 2500
    )
    job2$input <- rhfmt(FileInput, type="sequence")
    job2$output <- rhfmt(FileOutput, type="sequence")
    job2$mon.sec <- 10
    job2$jobname <- FileOutput
    job2$readback <- FALSE  

    job.mr <- do.call("rhwatch", job2)  

    # The output from job2 is the input to job3
    FileInput <- FileOutput  

    FileOutput <- file.path(
      rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", "0.05", sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), paste("Spatial", i, ".bystation", sep="")
    )  

    # job3 is just changing the key from month and year to station.id from job2
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

    FileInput <- FileOutput  

  }

  #outer that calculate the weights
  
  if(par$outer > 1) {

    FileOutput <- file.path(
      rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", "0.05", sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), "residual"
    )  

    job4 <- list()
    job4$map <- expression({
      lapply(seq_along(map.values), function(r) {
        R.abs <- abs(with(map.values[[r]], tmax - spatial - seasonal - trend))
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

    FileOutput <- file.path(
      rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", "0.05", sep=""), 
      paste(par$loess, par$loop, par$type, sep="."), paste("Outer", o, ".bystation", sep="")
    )

    job5 <- list()
    job5$map <- expression({
      lapply(seq_along(map.values), function(r) {
        R.abs <- abs(with(map.values[[r]], tmax - spatial - seasonal - trend))
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

    FileInput <- FileOutput

  }

}

