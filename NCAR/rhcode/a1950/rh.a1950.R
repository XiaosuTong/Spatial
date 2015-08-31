job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
    sub <- subset(map.values[[r]], year >= 1950)
    if(sum(is.na(sub$resp)) != 576) {
      attr(sub, "location") <- attributes(map.values[[r]])$location
      rhcollect(map.keys[[r]], sub)
    }
  })
})
job$input <- rhfmt(
  file.path(rh.root, par$dataset, "All", "bystation"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.root, par$dataset, "a1950", "bystation"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 100,  #cdh3,4
  mapreduce.job.reduces = 100  #cdh5
)
job$readback <- FALSE
job$combiner <- TRUE
job$jobname <- file.path(rh.root, par$dataset, "a1950", "bystation")
job.mr <- do.call("rhwatch", job)


job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
    v <- map.values[[r]]
    v$station.id <- as.character(map.keys[[r]])
    v$elev <- as.numeric(attributes(v)$location["elev"])
    v$lon <- as.numeric(attributes(v)$location["lon"])
    v$lat <- as.numeric(attributes(v)$location["lat"])
    lapply(1:nrow(v), function(k){
      value <- subset(v[k,], select = -c(year, month))
      rhcollect(c(v$year[k], v$month[k]), value)
    })
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
    row.names(combined) <- NULL
    rhcollect(reduce.key, combined)
  }
)
job$input <- rhfmt(
  file.path(rh.root, par$dataset, "a1950", "bystation"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.root, par$dataset, "a1950", "bymonth"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 100,  #cdh3,4
  mapreduce.job.reduces = 100  #cdh5
)
job$readback <- FALSE
job$combiner <- TRUE
job$jobname <- file.path(rh.root, par$dataset, "a1950", "bymonth")
job.mr <- do.call("rhwatch", job)

