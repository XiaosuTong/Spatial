########################################################
## spatial loess fit at all stations after 1950 with 
## elevation as conditional parameter 
########################################################

# my.loess02 is the loess function that calculate kd-tree 
# nodes, and all necessary information for interpolation

#"/ln/tongx/Spatial/tmp/tmax/spatial/a1950/loess04"
#is the span=0.05, degree=2 with residuals only for stations
#with elevation in the model as condtional parametric

source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$span <- 0.05
par$degree <- 2
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
	lo.fit <- my.loess2( get(par$dataset) ~ lon + lat + elev, 
		data       = v, 
		degree     = par$degree, 
		span       = par$span,
		parametric = "elev"
	)
	fit <- my.predict.loess(
		object = lo.fit, 
    newdata = data.frame(
    	lon = v$lon, 
    	lat = v$lat,
    	elev = v$elev
    )
	)
	v$fitted <- fit
	v$year <- y
	v$month <- m
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
	file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess04"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72, 
	rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess04")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)

## change the key from c(year, month) to station.id
job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		m <- match(substr(map.keys[[r]][2],1,3), month.abb)
		map.values[[r]]$time <- (as.numeric(map.keys[[r]][1]) - 1950)*12 + m
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
		if(sum(!is.na(combined$tmax)) >= 50){
			rhcollect(reduce.key, combined)
		}
	}
)
job$input <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess04"), 
	type = "sequence"
) 
job$output <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess04.bystation"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72
)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess04.bystation")
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
		rh.datadir, par$dataset, "spatial", "a1950", "loess04.bystation"
	), 
	type = "sequence"
) 
job$output <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", "loess04.bystation.10pc"
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 10
)
job$mon.sec <- 5
job$jobname <- file.path(
	rh.datadir, par$dataset, "spatial", "a1950", "loess04.bystation.10pc"
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)
