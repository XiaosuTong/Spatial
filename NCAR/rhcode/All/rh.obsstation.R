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

rst <- rhread(file.path(root, par$dataset, "All", "stationcount"))[[1]][[2]]

library(plyr)
library(lattice)


trellis.device(
  device = postscript, 
  file = file.path(local.output, "obs.station.ps"), 
  color = TRUE, 
  paper = "letter"
)
  b <- qqmath( ~ log2(count)
    , data = rst
    , distr = qunif
    , scales = list(
        y = list(at = seq(4,10,1), cex = 1.2),
        x = list(cex = 1.2)
      )
    , xlab = list(label = "f-value", cex = 1.5)
    , ylab = list(label = "Log Base 2 Number of Observation", cex = 1.5) 
    , panel = function(x,...) {
        panel.abline(v = seq(0,1,0.2), h = seq(4,10,1), lwd=0.5, col = "lightgray")
        panel.abline(h = log2(1236), lwd = 0.5, lty = 2, col = "black")
        panel.text(x=0.2, y= log2(1236), "Full Obs")
        panel.qqmath(x, pch = 16, cex = 0.3,...)
      }
  )
  print(b)
dev.off()

trellis.device(
  device = postscript, 
  file = file.path(local.output, "obsrate.station.ps"), 
  color = TRUE, 
  paper = "letter"
)
  b <- qqmath( ~ count/1236
    , data = rst
    , distr = qunif
    , scales = list(
        y = list(cex = 1.2),
        x = list(cex = 1.2)
      )
    , xlab = list(label = "f-value", cex = 1.5)
    , ylab = list(label = "Rate of Valid Observation", cex = 1.5) 
    , panel = function(x,...) {
        panel.abline(v = seq(0,1,0.2), h = seq(0,1,0.2), lwd=0.5, col = "lightgray")
        #panel.abline(h = log2(1236), lwd = 0.5, lty = 2, col = "black")
        #panel.text(x=0.2, y= log2(1236), "Full Obs")
        panel.qqmath(x, pch = 16, cex = 0.3,...)
      }
  )
  print(b)
dev.off()

