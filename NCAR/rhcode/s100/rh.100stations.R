##########################################
##  create key-values for 100 stations  ##
##########################################
#######################################################
##  Aggregate 100 stations subset for visualization  ##
#######################################################
s100 <- function() {

  rhload(
    file.path(rh.root, par$dataset, "100stations", "Rdata", "100stations.RData")
  )
  
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      if(map.keys[[r]] %in% target) {
        value <- arrange(map.values[[r]], year, match(month, month.abb))
        rhcollect(map.keys[[r]], value)
      }
    })
  })
  job$parameters <- list(
    target = stations100
  )
  job$setup <- expression(
    map = {
      library(lattice)
      library(plyr)
    }
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "All", "bystation"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "bystation"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 50,  #cdh3,4
    mapreduce.job.reduces = 50  #cdh5
  )
  job$readback <- FALSE
  job$jobname <- file.path(rh.root, par$dataset, "100stations", "bystation")
  job.mr <- do.call("rhwatch", job)
  
  
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      value <- map.values[[r]]
      value$station.id <- map.keys[[r]]
      value$lat <- as.numeric(attributes(map.values[[r]])$location["lat"])
      value$lon <- as.numeric(attributes(map.values[[r]])$location["lon"])
      value$elev <- as.numeric(attributes(map.values[[r]])$location["elev"])
      rhcollect(1, value)
    })
  })
  job$reduce <- expression(
    pre = {
      combined <- data.frame()
    },
    reduce = {
      combined <- rbind(combined, do.call("rbind", reduce.values))
    },
    post = { 
      rhcollect(reduce.key, combined)
    }
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "bystation"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "aggregated"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 1,  #cdh3,4
    mapreduce.job.reduces = 1  #cdh5
  )
  job$readback <- FALSE
  job$copyFiles <- TRUE
  job$jobname <- file.path(rh.root, par$dataset, "100stations", "aggregated")
  job.mr <- do.call("rhwatch", job)

}


s100.STLfit <- function(sw, sd, tw, td, fcw=NULL, fcd=NULL) {

  tuning <- list(sw=sw, sd=sd, tw=tw, td=td, fcw=fcw, fcd=fcd)

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      value <- map.values[[r]]
      value$station.id <- map.keys[[r]]
      value$date <- 1:1236
      if (is.null(par$fcw)) {
        
        fit <- stlplus(
          x=value$resp, t=value$date, n.p=12, s.window=par$sw, s.degree=par$sd, 
          t.window=par$tw, t.degree=par$td, inner=10, outer=0
        )$data
        value <- cbind(value, subset(fit, select = c(seasonal, trend, remainder)))

      } else {

        fit <- do.call("cbind", stlplus(
          x=value$resp, t=value$date, n.p=12, s.window=par$sw, s.degree=par$sd, t.window=par$tw, 
          t.degree=par$td, fc.window=c(par$tw,par$fcw), fc.degree=c(par$td,par$fcd), inner=10, outer=0
        )[c("data","fc")])
        value <- cbind(value, subset(fit, select = -c(data.raw, data.trend, data.remainder, data.weights, data.sub.labels)))
        names(value)[grep("fc.fc", names(value))] <- c("fc.first", "fc.second")

      }
      rhcollect(1, value)
    })
  })
  job$reduce <- expression(
    pre = {
      combined <- data.frame()
    },
    reduce = {
      combined <- rbind(combined, do.call("rbind", reduce.values))
    },
    post = { 
      rhcollect(reduce.key, combined)
    }
  )
  job$parameters <- list(
    par = tuning
  )
  job$setup <- expression(
    map = {
      suppressMessages(library(stlplus, lib.loc=lib.loc))
    }
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "bystation"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "STL", paste("t",tuning$tw, "td", tuning$td, "_s", tuning$sw, "sd", tuning$sd, "_f", tuning$fcw, "fd", tuning$fcd, sep="")), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 1,  #cdh3,4
    mapreduce.job.reduces = 1  #cdh5
  )
  job$readback <- FALSE
  job$jobname <- file.path(rh.root, par$dataset, "100stations", "STL", paste("t",tuning$tw, "td", tuning$td, "_s", tuning$sw, "sd", tuning$sd, "_f", tuning$fcw, "fd", tuning$fcd, sep=""))
  
  job.mr <- do.call("rhwatch", job)

}