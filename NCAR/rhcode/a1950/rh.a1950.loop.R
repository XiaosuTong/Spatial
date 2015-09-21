

backfitSwap.station <- function(input, output) {

    # job is just changing the key from month and year to station.id from job22
    job <- list()
    job$map <- expression({
      lapply(seq_along(map.values), function(r) {
        map.values[[r]]$year <- map.keys[[r]][1]
        map.values[[r]]$month <- map.keys[[r]][2]
        lapply(1:nrow(map.values[[r]]), function(i){
          value <- map.values[[r]][i, ]
          rhcollect(as.character(value$station.id), value)
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
    job$combiner <- TRUE
    job$input <- rhfmt(input , type = "sequence")
    job$output <- rhfmt(output, type = "sequence")
    job$mapred <- list(mapred.reduce.tasks = 72)
    job$mon.sec <- 10
    job$jobname <- output  
    job$readback <- FALSE  

    job.mr <- do.call("rhwatch", job)

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
        if(!fc.flag) {
          v$trend <- 0
        } else {
          v$fc.first <- 0
          v$fc.second <- 0
        }
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
      location <- c(value[1, c("lon", "lat", "elev2", "station.id")])
      value <- subset(value, select = -c(station.id, lon, lat, elev2))
      attr(value, "location") <- location
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

backfitStart <- function(span, family, type, degree) {

  FileInput <- file.path(
    rh.root, par$dataset, "a1950", "bymonth.fit", family, type, degree, paste("sp", span[1], sep="")
  )

  FileOutput <- file.path(
    rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), "startbystation"
  )
  
  try(backfitSwap.station(input=FileInput, output=FileOutput))

}

backfitAll <- function(span, family, type, parameter, index, degree, inner, outer) {

  FileInput <- file.path(
    rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), "startbystation"
  )

  for(o in 1:outer) {
    for(i in 1:inner) {
      
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
  
      # The output from job22 is the input to job
      FileInput <- FileOutput  
  
      FileOutput <- file.path(
        rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), index, paste("spatial",i,".bystation",sep="")
      )    
  
      try(backfitSwap.station(input=FileInput, output=FileOutput))
  
      FileInput <- FileOutput  
  
    }
  
    #outer that calculate the weights
    
    if(outer > 0) {
      
      FileOutput <- file.path(
        rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), index, "residual"
      )  

      try(w <- backfitWeights(input = FileInput, output = FileOutput, fc.flag=TRUE))
  
      FileOutput <- file.path(
        rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), index,
        paste("Outer", o, ".bystation", sep="")
      )
      
      try(backfitRobust(input = FileInput, output = FileOutput, weight = w, fc.flag=TRUE))
      
      FileInput <- FileOutput
  
    }
  
  }
}

backfitConverge <- function(span, family, type, index, degree, fc.flag) {

  input <- rhls(
    file.path(rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), index)
  )$file
  
  files <- input[grep("spatial[0-9]", input)]

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      file <- Sys.getenv("mapred.input.file")
      key <- substr(unlist(strsplit(tail(strsplit(file, "/")[[1]],3)[2], "[.]")), 8, 9)
      if(fc.flag) {
        value <- with(map.values[[r]], (resp - spatial - fc.first - fc.second - data.seasonal)^2)
      } else {
        value <- with(map.values[[r]], (resp - spatial - trend - seasonal)^2)
      }
      rhcollect(as.numeric(key), value)
    })
  })
  job$reduce <- expression(
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
  job$parameters <- list(fc.flag=fc.flag)
  job$input <- rhfmt(files, type = "sequence")
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "backfitting", family, type, degree,
      paste("sp", span[1], sep=""), index, "converge"
    ), 
    type = "sequence"
  )
  job$mapred <- list(mapred.reduce.tasks = 8)
  job$mon.sec <- 10
  job$jobname <- file.path(
    rh.root, par$dataset, "a1950", "backfitting", family, type, degree,
    paste("sp", span[1], sep=""), index, "converge"
  )
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)

}

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

backfitResidcomp <- function(by, family, type, degree, span, index, fc.flag, comp) {

  input <- rhls(
    file.path(rh.root, par$dataset, "a1950", "backfitting", family, type, degree, paste("sp", span[1], sep=""), index)
  )$file

  if(by == "station") {
    files <- input[grep("spatial[0-9].", unlist(lapply(strsplit(input, "/"), "[[",14)))]
  } else {
    files <- input[grep("^spatial[0-9]$", unlist(lapply(strsplit(input, "/"), "[[",14)))]
  }

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      file <- Sys.getenv("mapred.input.file")
      index <- substr(unlist(strsplit(tail(strsplit(file, "/")[[1]],3)[2], "[.]")), 8, 9)
      if (comp == "resid") {
        if (fc.flag) {
          value <- data.frame(target=with(map.values[[r]], resp - fc.first - fc.second - data.seasonal))
        } else {
          value <- data.frame(target=with(map.values[[r]], resp - seasonal - trend))
        }
      } else if(comp == "residfit") {
        if (fc.flag) {
          value <- data.frame(
            fit=with(map.values[[r]], fc.first + fc.second + data.seasonal), 
            resid=with(map.values[[r]], resp - fc.first - fc.second - data.seasonal)
          )
        } else {
          value <- data.frame(
            fit=with(map.values[[r]], seasonal + trend),
            resid=with(map.values[[r]], resp - trend - seasonal)
          )
        }
      } else {
        value <- data.frame(target=with(map.values[[r]], get(comp)))
      }
      rhcollect(index, value)
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
  job$input <- rhfmt(files, type = "sequence")
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "backfitting", family, type, degree,
      paste("sp", span[1], sep=""), index, paste(comp,"compare", sep="")
    ), 
    type = "sequence"
  ) 
  job$parameters <- list(fc.flag=fc.flag, comp=comp)
  job$mapred <- list(mapred.reduce.tasks = 10)
  job$mon.sec <- 10
  job$jobname <- file.path(
    rh.root, par$dataset, "a1950", "backfitting", family, type, degree,
    paste("sp", span[1], sep=""), index, paste(comp,"compare", sep="")
  )
  job$readback <- FALSE
  job$combiner <- TRUE
  job.mr <- do.call("rhwatch", job)

}