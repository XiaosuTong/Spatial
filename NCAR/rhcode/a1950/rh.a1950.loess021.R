########################################################
## spatial loess fit at all stations after 1950 
## includes all stations 7,738. All stations will be used
## for plotting condtional on year*month
##
##  - my.loess02 is the loess function that calculate kd-tree 
##    nodes, and all necessary information for interpolation
##
##  - loess02:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess02"
##      is key-value pairs that key is c(year, month), value is loess 
##      fit with family, degree, and span w/o elevation at that month.
##
##  - loess04:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess04"
##      is key-value pairs that key is c(year, month), value is loess 
##      fit with family, degree, and span w/ elevation at that month.
##      
##  - loess02/loess04.bystation.all:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess02.bystation.all"
##      is key-value pairs that key is station.id, value is data.frame
##      for the station. All stations included.
##
##  - loess02/loess04.bystation.all.10pc:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess02.bystation.all.10pc"
##      merge loess02.bystation.all to 10 key-value pairs, key is meaningless
##      faster to combine all 10 pieces for ploting.
##
###################################################
source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$loess <- "loess04"
#par$family <- "gaussian"
#par$span <- 0.025
par$family <- "symmetric"
par$span <- 0.05
par$degree <- 2
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")

## change the key from c(year, month) to station.id, 
## includes all stations
job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		m <- match(substr(map.keys[[r]][2],1,3), month.abb)
		map.values[[r]]$time <- (as.numeric(map.keys[[r]][1]) - 1950)*12 + m
		map.values[[r]]$year <- map.keys[[r]][1]
		map.values[[r]]$month <- map.keys[[r]][2]
		lapply(1:dim(map.values[[r]])[1], function(i){
			value <- map.values[[r]][i, ]
			rhcollect(as.character(value$station.id), value)
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
		rhcollect(reduce.key, combined)
	}
)
job$input <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), par$loess
	), 
	type = "sequence"
)
job$output <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), paste(par$loess, "bystation.all", sep=".")
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72
)
job$mon.sec <- 5
job$jobname <- file.path(
	rh.datadir, par$dataset, "spatial", "a1950",
	par$family, paste("sp", par$span, sep=""), paste(par$loess, "bystation.all", sep=".")
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


## group all stations to 10 groups
job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r){
		rhcollect(sample(1:10,1), map.values[[r]])
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
		rhcollect(reduce.key, combined)
	}
)
job$input <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), 
		paste(par$loess, "bystation.all", sep=".")
	), 
	type = "sequence"
)
job$output <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), 
		paste(par$loess, "bystation.all", "10pc", sep=".")
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 10
)
job$mon.sec <- 5
job$jobname <- file.path(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), 
		paste(par$loess, "bystation.all", "10pc", sep=".")
	)
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)
