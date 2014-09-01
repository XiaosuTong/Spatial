#initialize the rhipe and setup the directories
source("~/Rhipe/rhinitial.R")
par <- list()
par$dataset <- "tmax"
par$N <- 1236
par$span <- 0.2
par$degree <- 2
par$machine <- "gacrux"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/rhcode/my.loess02.R")
#get the loess fit at each vertix of the kd tree(different tree for each month) using the whole dataset. 
#The output key is the c(year, month), and value is the fitted value of all vertixes for that time.

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
	v <- na.omit(v)
	lo.fit <- my.loess2( get(par$dataset) ~ lon + lat, 
		data    = v, 
		degree  = par$degree, 
		span    = par$span, 
	)
	#value1 is the fitted value and first derivatives of each vertices
	value1 <- setNames(
		data.frame(matrix(lo.fit$kd$vval, byrow = TRUE, ncol = 3)),
		c("fitted", "dlon", "dlat")
	)
	#value2 is the location of each vertices
	value2 <- setNames(
		lo.fit$kd$vert2,
		c("lon", "lat")
	)
	value <- cbind(value2, value1)
	rhcollect(c(y, m), value)	
  })
})
job$setup <- expression(
	map = {
		dyn.load("myloess2.so")
	    load(paste(par$dataset, "RData", sep="."))
	}
)
job$shared <- c(
	file.path(rh.datadir, par$dataset, "shareRLib", "myloess2.so"),
	file.path(rh.datadir, par$dataset, "Rdata", paste(par$dataset, "RData", sep="."))
)
job$parameters <- list(
	par = par,
	my.loess = my.loess2, 
	my.simple = my.simple2
)
job$input <- c(par$N, 100) 
job$output <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "loess02"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72, 
	rhipe_reduce_full_size = 10000
)
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "loess02")
job$mon.sec <- 10
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


