noMissing <- function() {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      if(sum(is.na(map.values[[r]]$resp)) == 0) {
        value <- arrange(map.values[[r]], year, match(month, month.abb))
        rhcollect(map.keys[[r]], value)
      }
    })
  })
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
    file.path(rh.root, par$dataset, "nomissing", "bystation"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 50,  #cdh3,4
    mapreduce.job.reduces = 50  #cdh5
  )
  job$readback <- FALSE
  job$jobname <- file.path(rh.root, par$dataset, "nomissing", "bystation")
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
    file.path(rh.root, par$dataset, "nomissing", "bystation"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "nomissing", "bymonth"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 50,  #cdh3,4
    mapreduce.job.reduces = 50  #cdh5
  )
  job$readback <- FALSE
  job$copyFiles <- TRUE
  job$jobname <- file.path(rh.root, par$dataset, "nomissing", "bymonth")
  job.mr <- do.call("rhwatch", job)

}
