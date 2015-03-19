#initialize the Rhipe package and the setup directories
source("~/Rhipe/rhinitial.R")
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
dataset <- "tmax"
#for each divided subset, predict the loess fitting at vertix of kdtree, 
#then recombined fitted value through all subsets. subsets are conducted 
#by using Near Exact Replicate idea. Randomly sampled one observation from 
#each cell to get one subset.
N <- 1236
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
	v <- subset(get(dataset), year == y & month == m)
	lap <- max(unlist(lapply(kd003$cells, nrow)))
	tmp <- llply(1:lap, function(i){
   	  do.call(rbind, llply(
   	  	.data = kd003$cells, 
		.fun = function(r){
           	if(nrow(r) >= i)
           	subset(r, sub == i)[,1:(ncol(r) - 1)]
   	  	})
   	  )
	})
	tmp[[16]] <- do.call(rbind, tmp[16:18])
	length(tmp) <- 16 
	lapply(1:length(tmp), function(r){
		key <- c(y, m)
		v.tmp <- v[
			with(v, which(station.id %in% tmp[[r]]$station.id)), 
			c("station.id", "elev", "lon", "lat", dataset)
		]
		lo.fit <- loess( get(dataset) ~ lon + lat, 
			data    = v.tmp, 
			degree  = 1,
			control = loess.control(surface = "direct")
		)
        value <- predict(
        	lo.fit, 
        	data.frame(lon = kd003$vertix$lon, lat = kd003$vertix$lat)
        )
        value <- cbind(kd003$vertix, value)
        names(value) <- c("lon","lat","fitted")		
		value$fac <- 1:1026
		rhcollect(key, value)
	})
  })
})
job$reduce <- expression(
	pre={
		combined <- data.frame()
	},
	reduce={
		combined <- rbind(combined, do.call(rbind, reduce.values))
		stopifnot(nrow(combined) == 1026*16)
		value <- ddply(
			.data = combined,
			.variable = "fac",
			.fun = summarise,
			lon = unique(lon),
			lat = unique(lat),
			mean = mean(fitted)
		)
		value <- value[, c("lon","lat","mean")]
	},
	post={
		rhcollect(reduce.key, value)
	}	
)
job$setup <- expression(
	map = {
		library(plyr)
		load("kd003.RData")
	    load(paste(dataset, "RData", sep="."))
	},
	reduce = {
		library(plyr)
	}
)
job$shared <- c(
	file.path(rh.datadir, dataset, "Rdata", "kd003.RData"),
	file.path(rh.datadir, dataset, "Rdata", paste(dataset, "RData", sep="."))
)
job$parameters <- list(
	dataset = dataset
)
job$input <- c(N, 242) 
job$output <- rhfmt(
	file.path(rh.datadir, dataset, "spatial", "DR.loess"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72, 
	rhipe_reduce_full_size = 10000
)
job$mon.sec <- 10
job$jobname <- file.path(rh.datadir, dataset, "spatial", "DR.loess")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)
