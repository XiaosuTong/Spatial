###############################################
##
##Get the loess fit at each vertix of the kd tree using the whole 
##dataset. The output key is the vertix id, and value is the 1236 
##fitted value for that vertix.
##
##initialize the rhipe and setup the directories
source("~/Rhipe/rhinitial.R")
##set up the parameters
par <- list()
par$dataset <- "tmax"
par$machine <- "gacrux"
par$N <- 1236
par$span <- 0.125
par$degree <- 2
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess01.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")

load(file.path(local.datadir, "info.RData"))
if(par$dataset == "precip"){
	info <- USpinfo
}else{
	info <- UStinfo
}
info$fack <- rep(1, nrow(info))
dyn.load("~/Projects/Spatial/NCAR/myloess/shareLib/myloess2.so")
par$kdwhole <- my.loess2(
	fack ~ lon + lat, 
	data    = info, 
	degree  = par$degree, 
	span    = par$span,
	normalize = FALSE
)$kd$vert2
names(par$kdwhole) <- c("lon", "lat")
##loess fit at each vertices
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
	lo.fit <- my.loess1( get(par$dataset) ~ lon + lat, 
		data    = v, 
		degree  = par$degree, 
		span    = par$span,
		control = loess.control(surface = "direct")
	)
	pred <- data.frame(
		fitted = my.predict.loess(
			lo.fit, 
      data.frame(lon = par$kdwhole$lon, lat = par$kdwhole$lat)
		)
	)
	value <- cbind(par$kdwhole, pred)
	value$year <- rep(y, nrow(value))
	value$month <- rep(m, nrow(value))
	value$fac <- seq_len(nrow(value))
	lapply(1:nrow(value), function(k){
		key <- value$fac[k]
		rhcollect(key, value[k, ])
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
	  load(paste(par$dataset, "RData", sep="."))
	  dyn.load("/home/shaula/u16/tongx/Projects/Spatial/NCAR/myloess/shareLib/myloess1.so")
	}
)
job$shared <- c(
	file.path(rh.datadir, par$dataset, "Rdata", paste(par$dataset, "RData", sep="."))
)
job$parameters <- list(
	par = par, 
	my.loess1 = my.loess1, 
	my.simple1 = my.simple1,
	my.predict.loess = my.predict.loess,
	my.predLoess = my.predLoess
)
job$input <- c(par$N, 100) 
job$output <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "loess01"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72, 
	rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "loess01")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


########################################
## STL fitting for each vertex
##
#a <- list(sw="periodic", sd=1, tw=241, td=1, inner=10, outer=0, flag=FALSE)
a <- list(
	sw="periodic", 
	sd=1, 
	tw=1855, 
	td=1, 
	fcw=c(1855,241), 
	fcd=c(1,1), 
	inner=10, 
	outer=0, 
	flag=TRUE
) 
par <- c(par, a)
job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
	tmp <- map.values[[r]]
	if(par$flag){
    	dr <- do.call("cbind", stl2(tmp$fitted, tmp$time,
    		n.p=12,
    		s.window = par$sw,
    		s.degree = par$sd,
    		t.window = par$tw,
   			t.degree = par$td,
   			fc.window = par$fcw,
   			fc.degree = par$fcd,
   			inner = par$inner,
    		outer = par$outer
		)[c("data","fc")])
		tmp <- cbind(
			tmp, 
			dr[, c(!(names(dr) %in% c("station.id", "data.raw", "data.sub.labels")))]
		)
		names(tmp)[grep("fc.fc", names(tmp))] <- c("fc.trend", "fc.second")
	}else{
    	dr <- stl2(tmp$fitted, tmp$time,
    		n.p = 12,
    		s.window = par$sw,
    		s.degree = par$sd,
    		t.window = par$tw,
    		t.degree = par$td,
    		inner = par$inner,
    		outer = par$outer)$data
		tmp <- cbind(
			tmp, 
			dr[, c(!(names(dr) %in% c("station.id", "raw", "sub.labels")))]
		)
	}
	rhcollect(map.keys[[r]], tmp)
  })
})
job$setup <- expression(
    map = {
		library(lattice)
   	library(yaImpute, lib.loc = lib.loc)
		library(stl2, lib.loc = lib.loc)
    },
)
job$parameters <- list(
	par = par
)
job$input <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "loess01"), 
	type = "sequence"
)
job$output <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "loess01.stl"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72, 
	rhipe_reduce_buff_size = 10000, 
	mapred.max.split.size = 200000
)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "loess01.stl")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


