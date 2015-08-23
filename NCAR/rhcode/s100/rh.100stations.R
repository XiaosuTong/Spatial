source("~/Rhipe/ross.initial.R")

Machine <- "rossmann"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

par <- list()
par$dataset <- "tmax"

if (Machine == "adhara") {
  root <- "/ln/tongx/Spatial/tmp"
} else if (Machine == "rossmann") {
  root <- "/wsc/tongx/Spatial/tmp"
}
if(par$dataset == "precip") {
  Nstations <- 11918
} else {
  Nstations <- 8125
}


rhload(
	file.path(root, par$dataset, "100stations", "Rdata", "100stations.RData")
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
  file.path(root, par$dataset, "All", "bystation"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(root, par$dataset, "100stations", "bystation"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 50,  #cdh3,4
  mapreduce.job.reduces = 50  #cdh5
)
job$readback <- FALSE
job$jobname <- file.path(root, par$dataset, "100stations", "bystation")
job.mr <- do.call("rhwatch", job)



job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
    value <- map.values[[r]]
    value$factor <- factor(
      x = rep(paste("Period", 1:9), c(rep(144,8),84)), 
      levels = paste("Period", c(9:1))
    )
    value$time <- c(rep(0:143,8), 0:83)
    value$station.id <- map.keys[[r]] 
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
  file.path(root, par$dataset, "100stations", "bystation"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(root, par$dataset, "100stations", "aggregated"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 1,  #cdh3,4
  mapreduce.job.reduces = 1  #cdh5
)
job$readback <- FALSE
job$copyFiles <- TRUE
job$jobname <- file.path(root, par$dataset, "100stations", "aggregated")
job.mr <- do.call("rhwatch", job)