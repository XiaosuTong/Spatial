source("~/Rhipe/rhinitial.R")
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
par <- list()
par$dataset <- "tmax"
par$N <- 1236
par$span <- 0.2
par$degree <- 2
if(par$dataset == "precip"){
	load(file.path(local.datadir, "USpinfo.RData"))
    info <- USpinfo
}else{
	load(file.path(local.datadir, "UStinfo.RData"))
    info <- UStinfo
}
rm(list=grep("US", ls(), value=T))
#load stations.RData, have to check this part
load(file.path(datadir, "stations.RData"))
stations.100 <- get(grep(par$dataset, ls(), value=T))
par$stations.100 <- stations.100$station.id

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
		year == y & month == m & station.id %in% par$stations.100
	)[, c("station.id", "elev", "lon", "lat", par$dataset)]
	lo.fit <- loess( get(par$dataset) ~ lon + lat, 
		data    = v, 
		degree  = par$degree, 
		span    = par$span
	)
	v$fitted <- lo.fit$fitted
	v$year <- rep(y, nrow(value))
	v$month <- rep(m, nrow(value))
	lapply(1:nrow(v), function(r) {
		key <- v$station.id[k]
		rhcollect(key, v[k , !(names(v) %in% c("station.id"))])
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
		if(nrow(combined) != 1236) {
			rhcounter("BUG", "1236", 1)
		}
		stopifnot(nrow(combined) != 1236)
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
	    load(paste(par$dataset, "RData", sep="."))
	}
)
job$shared <- c(
	file.path(rh.datadir, par$dataset, "Rdata", paste(par$dataset, "RData", sep="."))
)
job$parameters <- list(
	par = par,
)
job$input <- c(par$N, 242) 
job$output <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "100stations", "loess01"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 66, 
	rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "100stations", "loess01")
job$readback <- FALSE
job$noeval <- TRUE
job.mr <- do.call("rhwatch", job)
z <- rhex(job.mr, async = TRUE)