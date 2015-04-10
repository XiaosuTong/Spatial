source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$loess <- "loess04" # loess04 is w/ elevation
par$family <- "symmetric"
par$degree <- 2
par$type <- "same" # or "same", "decr"
par$parameters <- list(
	sw = "periodic",
	tw = 109,
	sd = 1,
	td = 2,
	inner = 2,
	outer = 0
) 
if (par$type == "same") {
  par$span <- rep(0.05, 10)
} else if (par$type == "incr") {
  par$span <- c(seq(0.01, 0.05, by=0.005), rep(0.05, 8))
} else {
  par$span <- c(rep(0.05,3), seq(0.045, 0.01, by=-0.005), rep(0.01, 5))
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

for(i in 1:length(par$span)) {
  
  FileOutput <- file.path(
    rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", "0.05", sep=""), 
    paste(par$loess, "loop", par$type, sep="."), paste("STL", i, sep="")
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
      }
      v.stl <- stl2(
        x = with(map.values[[r]], tmax - spatial), 
        t = map.values[[r]]$time, 
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
          subset(map.values[[r]], select = -c(station.id, elev, lon, lat, elev2)), 
          subset(v.stl, select = -c(weights, sub.labels, raw))
        )
      } else {
        value <- cbind(
          subset(map.values[[r]], select = -c(station.id, lon, lat, elev2)), 
          subset(v.stl, select = -c(weights, sub.labels, raw))
        )        
      }
      attr(value, "location") <- c(map.values[[r]][1, c("lon", "lat", "elev2", "station.id")])
      rhcollect(map.keys[[r]], value)
    })
  })
  job1$parameters <- list(
    parameters = par$parameters,
    dataset = par$dataset,
    first = i == 1
  )
  job1$setup <- expression(
    map = {
      library(lattice)
      library(yaImpute, lib.loc = lib.loc)
      library(stl2, lib.loc = lib.loc)
    }
  )
  job1$input <- rhfmt(FileInput , type = "sequence")
  job1$output <- rhfmt(FileOutput, type = "sequence")
  job1$mapred <- list(mapred.reduce.tasks = 72)
  job1$mon.sec <- 5
  job1$jobname <- FileOutput  
  job1$readback <- FALSE

  job.mr <- do.call("rhwatch", job1)
  
  arg <- list(degree = par$degree, span = par$span[i], family = par$family)
  
  # The output from job1 is the input to job2
  FileInput <- FileOutput
  
  # Second job output is the spatial fitting, output key is month and year
  FileOutput <- file.path(
    rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", "0.05", sep=""), 
    paste(par$loess, "loop", par$type, sep="."), paste("Spatial", i, sep="")
  )

  job2 <- list()
  job2$map <- expression({
    lapply(seq_along(map.values), function(r) {
      d_ply(
        .data = map.values[[r]],
        .variable = c("year","month"),
        .fun = function(k, station = map.keys[[r]]) {
          key <- c(unique(k$year), unique(as.character(k$month)))
          k$fitted <- with(k, tmax - trend - seasonal)
          value <- subset(k, select = -c(trend, seasonal, year, month))
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
      lo.fit <- my.loess2( fitted ~ lon + lat + elev2, 
        data    = subset(combine, !is.na(fitted)), 
        degree  = argumt$degree, 
        span    = argumt$span,
        family  = argumt$family,
        parametric = "elev2",
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
      combine$spatial <- fit # remainder-spatial will be useful later for checking the convergence issue
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
    paste(par$loess, "loop", par$type, sep="."), paste("Spatial", i, ".bystation", sep="")
  )

  # job3 is just changing the key from month and year to station.id from job2
  job3 <- list()
  job3$map <- expression({
    lapply(seq_along(map.values), function(r) {
      map.values[[r]]$year <- map.keys[[r]][1]
      map.values[[r]]$month <- map.keys[[r]][2]
      lapply(1:dim(map.values[[r]])[1], function(i){
        value <- subset(map.values[[r]][i, ], select = -c(remainder))
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

##############################################
## check the convergence of findal residuals##
##############################################
job4 <- list()
job4$map <- expression({
  lapply(seq_along(map.values), function(r) {
    file <- Sys.getenv("mapred.input.file")
    key <- substr(unlist(strsplit(tail(strsplit(file, "/")[[1]],3)[2], "[.]")), 8, 9)
    value <- with(map.values[[r]], (remainder-spatial)^2)
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

rst <- rhread("/ln/tongx/Spatial/tmp/tmax/spatial/a1950/symmetric/sp0.05/loess04.loop.incr/residuals")
result <- data.frame(idx=unlist(lapply(rst, "[[", 1)), MSE = unlist(lapply(rst, "[[", 2)))
xyplot(MSE ~ idx, data=result)

