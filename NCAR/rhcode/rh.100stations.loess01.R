source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 1236
par$span <- 0.06
par$degree <- 2
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")

#load stations.RData, have to check this part
load(file.path(local.datadir, "stations.RData"))
stations.100 <- get(grep(par$dataset, ls(), value=T))
par$stations.100 <- stations.100

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
	y <- i + 1894
	v <- subset(
		get(par$dataset), 
		year == y & month == m
	)[, c("station.id", "elev", "lon", "lat", par$dataset)]
	lo.fit <- my.loess2( get(par$dataset) ~ lon + lat, 
		data    = v, 
		degree  = par$degree, 
		span    = par$span
	)
	new <- subset(v, station.id %in% par$stations.100)
	fit <- my.predict.loess(
		object = lo.fit, 
    newdata = data.frame(
    	lon = new$lon, 
    	lat = new$lat
    )
	)
	new$fitted <- fit
	new$year <- rep(y, nrow(new))
	new$month <- rep(m, nrow(new))
	lapply(1:nrow(new), function(k) {
		key <- as.character(v$station.id[k])
		rhcollect(key, new[k, ])
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
		combined$month <- factor(
			combined$month, 
			levels = c(
				"Jan", "Feb", "Mar", "Apr", "May", "June", 
				"July", "Aug", "Sep", "Oct", "Nov", "Dec"
			)
		)
		combined <- combined[with(combined, order(year, month)), ]
		combined$time <- 0:1235
		rhcollect(reduce.key, combined)
	}	
)
job$setup <- expression(
	map = {
		dyn.load("/home/shaula/u16/tongx/Projects/Spatial/NCAR/myloess/shareLib/myloess2.so")
	  load(paste(par$dataset, "RData", sep="."))
	}
)
job$shared <- c(
	file.path(rh.datadir, par$dataset, "Rdata", paste(par$dataset, "RData", sep="."))
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
	file.path(rh.datadir, par$dataset, "spatial", "100stations", "loess02"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72, 
	rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "100stations", "loess02")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)