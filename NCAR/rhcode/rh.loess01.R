###############################################
##
##get the loess fit at each vertix of the kd tree using the whole 
##dataset. The output key is the vertix id, and value is the 1236 
##fitted value for that vertix.
##
##initialize the rhipe and setup the directories
source("~/Rhipe/rhinitial.R")
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/rhcode/my.loess.R")
##set up the parameters
par <- list()
par$dataset <- "tmax"
par$N <- 1236
par$span <- 0.2
par$degree <- 2
##construct the kd tree for given parameters of span and degree
##station information is in USpinfo and UStinfo object.
if(par$dataset == "precip"){
	load(file.path(local.datadir, "USpinfo.RData"))
    info <- USpinfo
}else{
	load(file.path(local.datadir, "UStinfo.RData"))
    info <- UStinfo
}
rm(list=grep("US", ls(), value=T))
source("~/Projects/Spatial/NCAR/code/spatial/kdtree.R")
kd <- kdtree(info[, c("station.id","lon","lat")], alpha = par$span/5)
rhsave(
	list = ("kd"), 
	file = file.path(rh.datadir, par$dataset, "Rdata", paste("kd", ".RData", sep=""))
)

##loess fit at each vertices
job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
	i <- ceiling(map.keys[[r]]/12)
	j <- map.keys[[r]] - (i - 1)*12
	month <- c("Jan","Feb","Mar","Apr","May","June","July","Aug", "Sep", "Oct", "Nov", "Dec")
	m <- month[j]
	y <- i + 1894
	v <- subset(get(par$dataset), year == y & month == m)[, c("station.id", "elev", "lon", "lat", par$dataset)]
	lo.fit <- my.loess( get(par$dataset) ~ lon + lat, 
		data    = v, 
		degree  = par$degree, 
		span    = par$span, 
		control = loess.control(surface = "direct")
	)
	value <- predict(lo.fit, data.frame(lon = kd$vertix$lon, lat = kd$vertix$lat))
	value <- cbind(kd$vertix, value)
	names(value) <- c("lon","lat","fitted")
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
	pre={
		combined <- data.frame()
	},
	reduce={
		combined <- rbind(combined, do.call(rbind, reduce.values))
	},
	post={
		if(nrow(combined) == 1236) {
			rhcounter("BUG", "1236", 1)
		}
		stopifnot(nrow(combined) == 1236)
		combined$month <- factor(combined$month, 
			levels = c("Jan", "Feb", "Mar", "Apr", "May", "June", 
    		"July", "Aug", "Sep", "Oct", "Nov", "Dec")
		)
		combined <- combined[with(combined, order(year, month)), ]
		combined$time <- 0:1235
		rhcollect(reduce.key, combined)
	}	
)
job$setup <- expression(
	map = {
		load("kd.RData")
	    load(paste(par$dataset, "RData", sep="."))
	}
)
job$shared <- c(
	file.path(rh.datadir, par$dataset, "Rdata", "kd.RData"),
	file.path(rh.datadir, par$dataset, "Rdata", paste(par$dataset, "RData", sep="."))
)
job$parameters <- list(par = par, my.loess = my.loess, my.simple = my.simple)
job$input <- c(par$N, 242) 
job$output <- rhfmt(file.path(rh.datadir, par$dataset, "spatial", "loess01"), type="sequence")
job$mapred <- list(mapred.reduce.tasks = 66, rhipe_reduce_buff_size=10000)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "loess01")
job$readback <- FALSE
job$noeval <- TRUE
job.mr <- do.call("rhwatch", job)
z <- rhex(job.mr, async = TRUE)


########################################
## STL fitting for each vertex
##
a <- list(sw="periodic", sd=1, tw=241, td=1, inner=10, outer=0, flag=FALSE)
#a <- list(sw="periodic", sd=1, tw=1855, td=1, fcw=c(1855,241), fcd=c(1,1), inner=10, outer=0, flag=TRUE) 
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
		tmp <- cbind(tmp, dr[, c(!(names(dr) %in% c("station.id", "data.raw", "data.sub.labels")))])
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
		tmp <- cbind(tmp, dr[, c(!(names(dr) %in% c("station.id", "raw", "sub.labels")))])
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
job$parameters <- list(par = par)
job$input <- rhfmt(file.path(rh.datadir, par$dataset, "spatial", "loess01"), type = "sequence")
job$output <- rhfmt(file.path(rh.datadir, par$dataset, "spatial", "loess01.stl"), type="sequence")
job$mapred <- list(
	mapred.reduce.tasks = 66, 
	rhipe_reduce_buff_size = 10000, 
	mapred.max.split.size = 200000 # the size of the data frame for one vertex is about this size
)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "loess01.stl")
job$readback <- FALSE
job$noeval <- TRUE
job.mr <- do.call("rhwatch", job)
z <- rhex(job.mr, async = TRUE)

