########################################################
## spatial loess fit at all stations after 1950 
## evaluate the residuals at grid points
########################################################

# my.loess02 is the loess function that calculate kd-tree 
# nodes, and all necessary information for interpolation

#"/ln/tongx/Spatial/tmp/tmax/spatial/a1950/loess01"
#is the span=0.05, degree=2 of residuals for grid points

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
	lo.fit <- my.loess2( get(par$dataset) ~ lon + lat, 
		data    = v, 
		degree  = par$degree, 
		span    = par$span
	)
	fit <- my.predict.loess(
		object = lo.fit, 
    newdata = data.frame(
    	lon = v$lon, 
    	lat = v$lat
    )
	)
	v$fitted <- fit
	v <- v[!is.na(v[[par$dataset]]), ] ## remove all NA for residual calculation
	resid.fit <- my.loess2( (get(par$dataset)-fitted) ~ lon + lat,
		data = v,
		degree  = par$degree, 
		span    = par$span
	)
	new.grid <- expand.grid(
		lon = seq(-124, -67, by = 0.5),
		lat = seq(25, 49, by = 0.5)
	)
	grid.fit <- my.predict.loess(
		object = resid.fit,
		newdata = data.frame(
			lon = new.grid$lon,
			lat = new.grid$lat
		)
	)
	new.grid$resid.fit <- grid.fit
	instate <- !is.na(map.where("state", new.grid$lon, new.grid$lat))
	value <- new.grid[instate, ]
	rhcollect(c(y, m), value)
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
	file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess01"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72, 
	rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess01")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)