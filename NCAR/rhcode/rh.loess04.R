#initialize the Rhipe package and the setup directories
source("~/Rhipe/rhinitial.R")
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
dataset <- "tmax"

N <- 1236
span <- 0.1
job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
	i <- ceiling(map.keys[[r]]/12)
	j <- map.keys[[r]] - (i - 1)*12
	month <- c("Jan","Feb","Mar","Apr","May","June","July","Aug", "Sep", "Oct", "Nov", "Dec")
	m <- month[j]
	y <- i + 1894
	v <- subset(get(dataset), year == y & month == m)
	lap <- max(unlist(lapply(kd001$cells, nrow)))
	tmp <- llply(1:lap, function(i){
   	  do.call(rbind, llply(.data = kd001$cells, 
		.fun = function(r){
           	if(nrow(r) >= i)
           	subset(r, sub == i)[,1:(ncol(r) - 1)]
   	  }))
	})
	lapply(1, function(r){
		key <- c(y, m)
		v.tmp <- v[with(v, which(station.id %in% tmp[[r]]$station.id)), c("station.id", "elev", "lon", "lat", dataset)]
		lo.fit <- loess( get(dataset) ~ lon + lat, 
			data    = v.tmp, 
			degree  = 1,
			span    = span,
			control = loess.control(surface = "direct")
		)
	    value <- predict(lo.fit, data.frame(lon = v.tmp$lon, lat = v.tmp$lat))		
        value <- cbind(v.tmp, fitted = value)
		rhcollect(key, value)
	})
  })
})
job$setup <- expression(
	map = {
		library(plyr)
		load("kd001.RData")
	    load(paste(dataset, "RData", sep="."))
	}
)
job$shared <- c(
	file.path(rh.datadir, dataset, "Rdata", "kd001.RData"),
	file.path(rh.datadir, dataset, "Rdata", paste(dataset, "RData", sep="."))
)
job$parameters <- list(dataset = dataset, span = span)
job$input <- c(N, 242) 
job$output <- rhfmt(file.path(rh.datadir, dataset, "spatial", "loess04"), type="sequence")
job$mapred <- list(mapred.reduce.tasks = 72, rhipe_reduce_full_size=10000)
job$mon.sec <- 10
job$jobname <- file.path(rh.datadir, dataset, "spatial", "loess04")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)
