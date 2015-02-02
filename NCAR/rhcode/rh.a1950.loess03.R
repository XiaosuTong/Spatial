########################################################
## spatial loess fit at all stations after 1950 with 
## elevation as conditional parameter 
## evaluate the residuals at grid points lon*lat
########################################################

# my.loess02 is the loess function that calculate kd-tree 
# nodes, and all necessary information for interpolation

#"/ln/tongx/Spatial/tmp/tmax/spatial/a1950/loess03"
#is the span=0.05, degree=2 with residuals for grid points on lon*lat
#model is tmax ~ lat + lon + elev which elev is conditional parametric

source("~/Rhipe/rhinitial.R")
par <- list()
library(maps)
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$span <- 0.05
par$degree <- 2
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")

new.grid <- expand.grid(
	lon = seq(-124, -67, by = 0.25),
	lat = seq(25, 49, by = 0.25)
)
instate <- !is.na(map.where("state", new.grid$lon, new.grid$lat))
new.grid <- new.grid[instate, ]

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
		v <- v[!is.na(v[[par$dataset]]), ] ## remove all NA for residual calculation
		resid.fit <- my.loess2( (get(par$dataset)-fitted) ~ lon + lat,
			data    = v,
			degree  = par$degree, 
			span    = par$span
		)
		grid.fit <- my.predict.loess(
			object  = resid.fit,
			newdata = data.frame(
				lon = new.grid$lon,
				lat = new.grid$lat
			)
		)
		new.grid$resid.fit <- grid.fit
		rhcollect(c(y, m), new.grid)
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
	my.predLoess = my.predLoess,
	new.grid = new.grid
)
job$input <- c(par$N, 100) 
job$output <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess03"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72, 
	rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 5
job$jobname <- file.path(
	rh.datadir, par$dataset, "spatial", "a1950", "loess03"
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)
