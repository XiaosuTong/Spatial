#initialize the rhipe and setup the directories
source("~/Rhipe/rhinitial.R")
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
dataset <- "tmax"
#get the loess fit at each non-NA location of the using the whole dataset. 
#The output key is the c(year, month), and value is the fitted value of locations.
N <- 1236
span <- 0.05
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
		get(dataset), 
		year == y & month == m
	)[, c("station.id", "elev", "lon", "lat", dataset)]
	lo.fit <- loess( get(dataset) ~ lon + lat, 
		data    = v, 
		degree  = 1, 
		span    = span, 
		control = loess.control(surface = "direct")
	)
	#pred <- predict(lo.fit, data.frame(lon=v$lon, lat=v$lat))
	v <- na.omit(v)
	value <- cbind(v, lo.fit$fitted)
	names(value)[6] <- "predict"
	value$resid <- value[,dataset] - value[,"predict"]
	rhcounter("LOESS","_ALL_",1)
	rhcollect(c(y, m), value)
  })
})
job$reduce <- expression(
	pre={
		MSE <- 0
	},
	reduce={
    	MSE <- MSE + unlist(lapply(reduce.values, function(r) sum(r$resid^2)))
	},
	post={
    	rhcollect(reduce.key, MSE)
   	}
)
job$setup <- expression(
	map = {
	    load(paste(dataset, "RData", sep="."))
	}
)
job$shared <- c(
	file.path(rh.datadir, dataset, "Rdata", paste(dataset, "RData", sep="."))
)
job$parameters <- list(
	dataset = dataset, 
	span = span
)
job$input <- c(N, 242) 
job$output <- rhfmt(
	file.path(rh.datadir, dataset, "spatial", "loess03"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72, 
	rhipe_reduce_full_size = 10000
)
job$jobname <- file.path(rh.datadir, dataset, "spatial", "loess03")
job$mon.sec <- 10
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)

