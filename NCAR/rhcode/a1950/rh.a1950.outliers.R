outliersStatAll <- function(input, output, lim=5) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      value <- subset(map.values[[r]], flag == 1)
      value$resid <- with(value, remainder - spafit)
      value <- subset(value, resid <= -lim | resid >= lim)
      rhcollect(1, unique(value$station.id))
    })
  })
  job$reduce <- expression(
    pre = {
      combine <- character()
    },
    reduce = {
      combine <- c(combine, unlist(reduce.values))
    },
    post = {
      rhcollect(reduce.key, unique(combine))
    }
  )
  job$setup <- expression(
    map = {
      library(plyr, lib.loc=lib.loc)
    }
  )
  job$mapred <- list(
    mapred.reduce.tasks = 1,
    mapred.tasktimeout = 0
  )
  job$parameters <- list(
    lim =lim
  )
  job$input <- rhfmt(input, type="sequence")
  job$output <- rhfmt(output, type="sequence")
  job$mon.sec <- 10
  job$jobname <- output
  job$readback <- TRUE
  job.mr <- do.call("rhwatch", job)    

  outliers.a1950.stations <- job.mr[[1]][[2]]
  rhsave(outliers.a1950.stations, file=file.path(rh.root, par$dataset, "a1950","Rdata","outliers.a1950.RData"))

}

outliersTotal <- function(input, output, lim=2) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      value <- subset(map.values[[r]], flag == 1)
      value$resid <- with(value, remainder - spafit)
      value <- subset(value, resid <= -lim | resid >= lim)
      value <- subset(value, select = c(station.id, resp, remainder, lon, lat, elev2, spafit))
      value$month <- map.keys[[r]][2]
      value$year <- as.numeric(map.keys[[r]][1])
      rhcollect(2, data.frame(count=sum(map.values[[r]]$flag)))
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
      rhcollect(reduce.key, unique(combine))
    }
  )
  job$parameters <- list(
    lim =lim
  )
  job$setup <- expression(
    map = {
      library(plyr, lib.loc=lib.loc)
    }
  )
  job$mapred <- list(
    mapred.reduce.tasks = 1,
    mapred.tasktimeout = 0
  )
  job$input <- rhfmt(input, type="sequence")
  job$output <- rhfmt(output, type="sequence")
  job$mon.sec <- 10
  job$jobname <- output
  job$readback <- TRUE
  job.mr <- do.call("rhwatch", job)    

  print(sum(job.mr[[1]][[2]]$count))
  return(job.mr[[2]][[2]])

}


outliersCount <- function(input, output, lim=2, by) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      resid <- with(map.values[[r]], remainder - spafit)
      outlier <- sum(resid <=-lim | resid >= lim)
      impute <- sum(map.values[[r]]$flag == 0)
      if (by == "station") {
        value <- data.frame(station.id = map.keys[[r]], out = outlier, imp = impute)
      } else {
        value <- data.frame(year = map.keys[[r]][1], month = map.keys[[r]][2], out = outlier, imp = impute)
      }
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
      rhcollect(reduce.key, unique(combine))
    }
  )
  job$parameters <- list(
    lim =lim,
    by = by
  )
  job$setup <- expression(
    map = {
      library(plyr, lib.loc=lib.loc)
    }
  )
  job$mapred <- list(
    mapred.reduce.tasks = 1,
    mapred.tasktimeout = 0
  )
  job$input <- rhfmt(input, type="sequence")
  job$output <- rhfmt(output, type="sequence")
  job$mon.sec <- 10
  job$jobname <- FileOutput
  job$readback <- TRUE
  job.mr <- do.call("rhwatch", job)  

  return(job.mr[[1]][[2]])
  
}


localRadius <- function() {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      value <- map.values[[r]]
      value$outFlag <- with(value, remainder - spafit) <= -5 | with(value, remainder - spafit) >= 5
      x1 <- subset(value, outFlag, select=c(lon, lat))
      x2 <- subset(value, select=c(lon, lat))
      if(nrow(x1) > 0) {
        dist <- as.data.frame(rdist.earth(x2, x1, miles = FALSE, R = NULL))
        R <- as.numeric(apply(dist, 2, function(r) {sort(round(r, 3))[floor(0.015*7738)]}))
        v <- subset(value, outFlag, select=c(lon, lat, elev2, resp, fitted, seasonal, trend, remainder, spafit))
        v$radius <- R
        v$month <- map.keys[[r]][2]
        v$year <- as.numeric(map.keys[[r]][1])
        rhcollect(1, v)
      }    
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
  job$setup <- expression(
    map = {
      library(plyr, lib.loc=lib.loc)
      suppressMessages(library(spam, lib.loc=lib.loc))
      suppressMessages(library(grid, lib.loc=lib.loc))
      suppressMessages(library(fields, lib.loc=lib.loc))
    }
  )
  job$mapred <- list(
    mapred.reduce.tasks = 1,
    mapred.tasktimeout = 0,
    #mapreduce.task.timeout = 0
  )
  job$input <- rhfmt(FileInput, type="sequence")
  job$output <- rhfmt(FileOutput, type="sequence")
  job$mon.sec <- 10
  job$jobname <- FileOutput
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)  

}


outliersStatTop <- function(input, output, lim=2, top=2^4) {
  
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      resid <- with(map.values[[r]], remainder - spafit)
      outlier <- sum(resid <=-lim | resid >= lim)
      if (outlier > top) {
        rhcollect(1, map.keys[[r]])
      }
    })
  })
  job$reduce <- expression(
    pre = {
      combine <- character()
    },
    reduce = {
      combine <- c(combine, unlist(reduce.values))
    },
    post = {
      rhcollect(reduce.key, combine)
    }
  )
  job$parameters <- list(
    lim = lim,
    top = top
  )
  job$mapred <- list(
    mapred.reduce.tasks = 1,
    mapred.tasktimeout = 0
  )
  job$input <- rhfmt(input, type="sequence")
  job$output <- rhfmt(output, type="sequence")
  job$mon.sec <- 10
  job$jobname <- FileOutput
  job$readback <- TRUE
  job.mr <- do.call("rhwatch", job)  

  return(job.mr[[1]][[2]])
  
}
