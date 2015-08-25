rhload(
  file.path(root, par$dataset, "100stations", "Rdata", "100stations.RData")
)
##########################################
##  create key-values for 100 stations  ##
##########################################
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

#######################################################
##  Aggregate 100 stations subset for visualization  ##
#######################################################
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


