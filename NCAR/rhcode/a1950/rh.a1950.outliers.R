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
      value <- subset(map.values[[r]], flag == 1)
      resid <- with(value, remainder - spafit)
      outlier <- sum(resid <=-lim | resid >= lim)
      if (by == "station") {
        value <- data.frame(station.id = map.keys[[r]], out = outlier)
      } else {
        value <- data.frame(year = map.keys[[r]][1], month = map.keys[[r]][2], out = outlier)
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
  job$jobname <- output
  job$readback <- TRUE
  job.mr <- do.call("rhwatch", job)  

  return(job.mr[[1]][[2]])
  
}

## MissCount is about to count the outliers which
## remainder is small but spafit is very large
outliersMissCount <- function(input, output, lim=2) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      value <- subset(map.values[[r]], flag == 1)
      value$resid <- with(value, remainder - spafit)
      value <- subset(value, resid <= -lim | resid >= lim)
      value <- subset(value, remainder < 1 & remainder > -1)
      if (nrow(value) > 0) {
        value$station.id  <- map.keys[[r]]
        value$lon <- attributes(map.values[[r]])$loc[1]
        value$lat <- attributes(map.values[[r]])$loc[2]
        value$elev2 <- attributes(map.values[[r]])$loc[3]      
        rhcollect(1, value)
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
      rhcollect(reduce.key, unique(combine))
    }
  )
  job$parameters <- list(
    lim =lim
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

  return(job.mr[[1]][[2]])
  
}

##  outliersTop is trying to find either a month or a station
##  which contains outliers more than top, or contains a outliers larger
##  than ORtop.
outliersTop <- function(input, output, lim=2, top=2^4, ORtop=0, by) {
  
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      value <- subset(map.values[[r]], flag == 1)
      resid <- with(value, remainder - spafit)
      outlier <- sum(resid <= -lim | resid >= lim)
      if (outlier >= top) {
        rhcollect(1, map.keys[[r]])
        rhcounter("outliers", "toomuch", 1)
      }
      if (max(resid) >= ORtop | min(resid) <= -ORtop) {
        rhcollect(1, map.keys[[r]])
        rhcounter("outliers","toolarge", 1)
      }
    })
  })
  if(by == "station") {
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
  } else {
    job$reduce <- expression(
      pre = {
        combine <- data.frame()
      },
      reduce = {
        combine <- rbind(combine, do.call("rbind", reduce.values))
      },
      post = {
        names(combine) <- c("year","month")
        rhcollect(reduce.key, combine)
      }
    )
  } 
  job$parameters <- list(
    lim = lim,
    top = top,
    ORtop = ORtop
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

  if(by=="station") { 
    outliers.a1950.stations <- (job.mr[[1]][[2]])
    rhsave(outliers.a1950.stations, file="/ln/tongx/Spatial/tmp/tmax/a1950/Rdata/outliersStatTop.a1950.RData")
  } else {
    outliers.a1950.stations <- (job.mr[[1]][[2]])
    rhsave(outliers.a1950.stations, file="/ln/tongx/Spatial/tmp/tmax/a1950/Rdata/outliersMonthTop.a1950.RData")    
  }
  
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