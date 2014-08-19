#initialize the rhipe and setup the directories
source("~/Rhipe/rhinitial.R")
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
dataset <- "tmax"
#get the loess fit at each vertix of the kd tree(different tree for each month) using the whole dataset. 
#The output key is the c(year, month), and value is the fitted value of all vertixes for that time.
N <- 1236
source("~/Projects/Spatial/NCAR/code/spatial/kdtree.R")
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
	v <- na.omit(v)
	kd <- kdtree(
		df = v[,c("station.id","lon","lat")], 
		alpha = 0.05
	)
	lo.fit <- loess( get(dataset) ~ lon + lat, 
		data    = v, 
		degree  = 1, 
		span    = 0.1, 
		control = loess.control(surface = "direct")
	)
	value <- predict(
		lo.fit, 
		data.frame(lon = kd$vertix$lon, lat = kd$vertix$lat)
	)
	value <- cbind(kd$vertix, value)
	names(value) <- c("lon","lat","fitted")
	value$fac <- 1:nrow(kd$vertix)
	rhcollect(c(y, m), value)	
  })
})
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
	kdtree = kdtree
)
job$input <- c(N, 242) 
job$output <- rhfmt(
	file.path(rh.datadir, dataset, "spatial", "loess02"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72, 
	rhipe_reduce_full_size = 10000
)
job$jobname <- file.path(rh.datadir, dataset, "spatial", "loess02")
job$mon.sec <- 10
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


