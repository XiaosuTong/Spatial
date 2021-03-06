########################################################
## spatial loess fit at all stations after 1950 only use
## longitude and latitude.
##
##  - my.loess02: 
##      is the loess function that calculate kd-tree 
##      nodes, and all necessary information for interpolation
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
##  - loess02.bystation:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess02.bystation"
##      is key-value pairs that key is station.id, value is data.frame
##      for the station. Only stations have over 300 obs got kept.
##
##  - loess02.bystation.10pc:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess02.bystation.10pc"
##      merge loess02.bystation to 10 key-value pairs, key is meaningless
##      faster to combine all 10 pieces for ploting.
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
par$loess <- "loess02"
#par$family <- "gaussian"
#par$span <- 0.025
par$family <- "symmetric"
par$span <- 0.05
par$degree <- 2
par$all <- TRUE
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")

job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
	i <- ceiling(map.keys[[r]]/12)
	j <- map.keys[[r]] - (i - 1)*12
	month <- c(
		"Jan","Feb","Mar","Apr","May","June",
		"July","Aug", "Sep", "Oct", "Nov", "Dec"
	)
	m <- month[j]
	y <- i + 1949
	v <- subset(
		get(paste(par$dataset, "a1950", sep=".")), 
		year == y & month == m
	)[, c("station.id", "elev", "lon", "lat", par$dataset)]
	lo.fit <- my.loess2( get(par$dataset) ~ lon + lat, 
		data    = v, 
		degree  = par$degree, 
		span    = par$span,
		family  = par$family,
		control = loess.control(iterations = 10)
	)
	fit <- my.predict.loess(
		object = lo.fit, 
    newdata = data.frame(
    	lon = v$lon, 
    	lat = v$lat
    )
	)
	v$fitted <- fit
	v$station.id <- as.character(v$station.id)
	rhcollect(c(y, m), v)
  })
})
job$setup <- expression(
	map = {
		system("chmod 777 myloess2.so")
		dyn.load("myloess2.so")
	  load(paste(par$dataset, "a1950", "RData", sep="."))
	  library(maps, lib.loc = lib.loc)
	}
)
job$shared <- c(
	file.path(
		rh.datadir, par$dataset, "shareRLib", "myloess2.so"
	),
	file.path(
		rh.datadir, par$dataset, "a1950", "Rdata", 
		paste(par$dataset, "a1950", "RData", sep=".")
	)
)
job$parameters <- list(
	par = par,
	my.loess2 = my.loess2,
	my.simple2 = my.simple2,
	my.predict.loess = my.predict.loess,
	my.predLoess = my.predLoess
)
job$input <- c(par$N, 100) 
job$output <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), par$loess
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72, 
	rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 5
job$jobname <- file.path(
	rh.datadir, par$dataset, "spatial", "a1950",
	par$family, paste("sp", par$span, sep=""), par$loess
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)

## change the key from c(year, month) to station.id, 
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
if (!par$all) {
  job$reduce <- expression(
  	pre = {
	  	combined <- data.frame()
	  },
	  reduce = {
		  combined <- rbind(combined, do.call(rbind, reduce.values))
	  },
	  post = {
		  if(sum(!is.na(combined$tmax)) >= 300){
			  rhcollect(reduce.key, combined)
		  } 
	  }
  )
  file <- "bystation"
} else {
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
  file <- "bystation.all"
}
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
		par$family, paste("sp", par$span, sep=""), paste(par$loess, file, sep=".")
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72
)
job$mon.sec <- 5
job$jobname <- file.path(
	rh.datadir, par$dataset, "spatial", "a1950",
	par$family, paste("sp", par$span, sep=""), paste(par$loess, file, sep=".")
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
		paste(par$loess, file, sep=".")
	), 
	type = "sequence"
)
job$output <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), 
		paste(par$loess, file, "10pc", sep=".")
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
		paste(par$loess, file, "10pc", sep=".")
	)
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)
