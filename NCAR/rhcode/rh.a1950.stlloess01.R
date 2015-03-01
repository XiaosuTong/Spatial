########################################################
## spatial loess fit of stl remainder at stations after 1950 
##
## - Find all locations that can do stl2
########################################################
source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$family <- "symmetric"
par$sw <- "periodic"
par$sd <- 1
par$tw <- 241
par$td <- 1
par$inner <- 10
par$outer <- 0
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

# load the stations vector after 1950, there are 7,738 stations
load(file.path(local.datadir, "stations.a1950.RData"))
month.or <- c(
		"Jan","Feb","Mar","Apr","May","June",
		"July","Aug", "Sep", "Oct", "Nov", "Dec"
)
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")

job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
    if (map.keys[[r]] %in% stations) {         # only includes 7,738 stations
      rhcounter("stations", "_a1950_", 1)
      value <- subset(map.values[[r]], year >= 1950)
      if (sum(!is.na(value$resp)) >= 300) {   # only includes 4,978 stations each of which has 300+ obs
        rhcounter("stations", "_300_", 1)
        value$month <- factor(value$month, levels = month.or)
        value <- value[with(value, order(year, month)),]
        if (nrow(value) != 576) rhcounter("stations", "_576_", 1)
        value$time <- 1:nrow(value)
        fit <- stl2(
          x = value$resp, 
          t = value$time, 
          n.p = 12, 
          s.window = sw, 
          s.degree = sd, 
          t.window = tw, 
          t.degree = td, 
          inner = inner, 
          outer = outer
        )$data
        value <- cbind(
        	subset(value,  select = c(month, year, resp, time)),
          subset(fit, select = -c(weights, sub.labels, raw))
        )
        rhcollect(map.keys[[r]], value)
      }
    }
  })
})
job$setup <- expression(
	map = {
    library(lattice)
    library(yaImpute, lib.loc = lib.loc)
    library(stl2, lib.loc = lib.loc)
	}
)
job$parameters <- list(
  stations = stations.a1950.tmax,
  sw = par$sw,
  sd = par$sd,
  tw = par$tw,
  td = par$td,
  inner = par$inner,
  outer = par$outer
)
job$input <- rhfmt(
  file.path(rh.datadir, par$dataset, "data", "bystation"),
  type = "sequence"
) 
job$output <- rhfmt(
	file.path(rh.datadir, par$dataset, "stl", "a1950", "E1"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72
)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "stl", "a1950", "E1")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)

job <- list()
job$map <- expression({

})
job$reduce <- expression(
	pre = {

	},
	reduce = {

	},
	post = {

	}
)
job$input <- c()
job$output <- rhfmt("/ln/tongx/userhipe/partition", type="sequence")
job$mapred <- list(
	mapred.reduce.tasks = 72
)
job$mon.sec <- 5
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)
