#######################################################################
##Order the stations based on the deviation from the normal distribution
##For each station, given lap, given group, the distribution of 601 
##prediction error is compared with normal distribution. The sum of abs
##deviation over all lap for one station under given group is calculated.
#######################################################################
datadir <- "~/Projects/Spatial/NCAR/RData/"
dataset <- "tmax"
par <- list()
par$machine <- "gacrux"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Rhipe/rhinitial.R")

#load(paste(datadir, dataset, "div.stations.RData", sep = ""))
for(index in c("E5")){
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r){
      v <- map.values[[r]]
      v$group <- map.keys[[r]][1]
      myfun <- function(data){
        yy <- quantile(data$residual, c(0.25, 0.75))
        xx <- qnorm(c(0.25, 0.75))
        r <- diff(yy)/diff(xx)
        x <- qnorm(ppoints(length(data$residual)))
        y <- r*x + yy[1] - xx[1]*r
        div <- sum(abs(sort(data$residual) - y))
        rhcollect(c(unique(data$group) ,unique(data$station.id)), div)
      }
      tmp <- d_ply(
        .data = v,
        .variables = "station.id",
        .fun = myfun
      )
    })
  })
  job$reduce <- expression(
    pre = {
      total <- 0
    },
    reduce = {
      total <- total + sum(unlist(reduce.values))
    },
    post = {
      rhcollect(reduce.key, total)
    }
  )
  job$setup <- expression(
    map = {
      library(plyr)
    }
  )
  job$setup <- expression(
    map = {
      library(plyr)
    }
  )
  job$input <- rhfmt(
    file.path(rh.datadir, dataset,"100stations","sharepredict",index,"36.lap.station"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.datadir, dataset,"100stations","sharepredict",index,"orderstations"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 72
  )
  job$jobname <- paste(dataset, "order stations")
  job$readback <- FALSE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r){
      rhcollect(map.keys[[r]][1], c(station.id = map.keys[[r]][2], div = map.values[[r]]))
    })
  })
  job$reduce <- expression(
    pre = {
      combined <- data.frame()
    },
    reduce = {
      combined <- rbind(combined, do.call(rbind, reduce.values)) 
    },
    post = {
      combined$div <- as.numeric(as.character(combined$div))
      combined <- combined[order(combined$div),]
      rhcollect(reduce.key, combined)
    }
  )
  job$input <- rhfmt(
    file.path(rh.datadir, dataset,"100stations","sharepredict",index,"orderstations"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.datadir, dataset,"100stations","sharepredict",index,"group.orderstations"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 72
  )
  job$jobname <- paste(dataset, "order stations")
  job$readback <- FALSE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)
}


