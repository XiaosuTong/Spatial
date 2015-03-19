########################################################
## spatial loess fit at all stations after 1950 with 
## elevation as conditional parameter 
## evaluate the residuals at grid points lon*lat
##
##  - model:
##      tmax ~ lat + lon + log2(elev+128) which elev is conditional parametric
##
##  - my.loess02: 
##      is the loess function that calculate kd-tree 
##      nodes, and all necessary information for interpolation
##
##  - loess03:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess03"
##      is key-value pairs that key is c(year, month), value is loess 
##      smoothing with family, degree, and span of the residuals at 
##      grid points in that month.
##
##  - loess03.rob:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess03.rob"
##      is key-value pairs that key is c(year, month), value is robust loess 
##      smoothing with family, degree, and span of the residuals at 
##      grid points in that month.
##
########################################################
source("~/Rhipe/rhinitial.R")
par <- list()
library(maps)
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$degree <- 2
par$loess <- "loess03.rob"
#par$family <- "gaussian"
#par$span <- 0.05
par$span <- 0.05
par$family <- "symmetric"
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
		v$elev2 <- log2(v$elev + 128)
		lo.fit <- my.loess2( get(par$dataset) ~ lon + lat + elev2, 
			data       = v, 
			degree     = par$degree, 
			span       = par$span,
			parametric = "elev2",
			family     = par$family,
			control = loess.control(iterations = 10)
		)
		fit <- my.predict.loess(
			object = lo.fit, 
    	newdata = data.frame(
    		lon = v$lon, 
    		lat = v$lat,
    		elev2 = v$elev2
    	)
		)
		v$fitted <- fit
		v <- v[!is.na(v[[par$dataset]]), ] ## remove all NA for residual calculation
		resid.fit <- my.loess2( (get(par$dataset)-fitted) ~ lon + lat,
			data    = v,
			degree  = par$degree, 
			span    = par$span,
			control = loess.control(iterations = 10)
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
