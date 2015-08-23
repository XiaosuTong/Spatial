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

job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
    value <- data.frame(
      count = sum(!is.na(map.values[[r]]$resp)),
      station.id = map.keys[[r]],
      elev = as.numeric(attributes(map.values[[r]])$location["elev"]),
      lon = as.numeric(attributes(map.values[[r]])$location["lon"]),
      lat = as.numeric(attributes(map.values[[r]])$location["lat"]),
      stringsAsFactors = FALSE
    )
    rhcollect(1, value)
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
  file.path(root, par$dataset, "All", "bystation"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(root, par$dataset, "All", "stationcount"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 1,  #cdh3,4
  mapreduce.job.reduces = 1  #cdh5
)
job$readback <- FALSE
job$combiner <- TRUE
job$jobname <- file.path(root, par$dataset, "All", "stationcount")
job.mr <- do.call("rhwatch", job)