########################################################
## spatial loess fit at all stations after 1950 
## includes only stations  4,978
##  - my.loess02: 
##      is the loess function that calculate kd-tree 
##      nodes, and all necessary information for interpolation
##
##  - loess02:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess02"
##      is key-value pairs that key is c(year, month), value is loess 
##      fit with family, degree, and span w/o elevation at that month.
##
##  - loess02.bystation:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess02.bystation"
##      is key-value pairs that key is station.id, value is data.frame
##      for the station. Only stations have over 300 obs got kept.
##
##  - loess02.bystation.10pc:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess02.bystation.10pc"
##      merge loess02.bystation to 10 key-value pairs, key is meaningless
##      faster to combine all 10 pieces for ploting.
##
###################################################
source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$loess <- "loess02"
#par$family <- "gaussian"
#par$span <- 0.025
par$family <- "symmetric"
par$span <- 0.05
par$degree <- 2
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")

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
	lo.fit <- my.loess2( get(par$dataset) ~ lon + lat, 
		data    = v, 
		degree  = par$degree, 
		span    = par$span,
		family  = par$family,
		control = loess.control(iterations = 10)
	)
	fit <- my.predict.loess(
		object = lo.fit, 
    newdata = data.frame(
    	lon = v$lon, 
    	lat = v$lat
    )
	)
	v$fitted <- fit
	v$station.id <- as.character(v$station.id)
	rhcollect(c(y, m), v)
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
	my.predLoess = my.predLoess
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

## change the key from c(year, month) to station.id, 
## only includes stations that have over 300 obs
job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		m <- match(substr(map.keys[[r]][2],1,3), month.abb)
		map.values[[r]]$time <- (as.numeric(map.keys[[r]][1]) - 1950)*12 + m
		map.values[[r]]$year <- map.keys[[r]][1]
		map.values[[r]]$month <- map.keys[[r]][2]
		lapply(1:dim(map.values[[r]])[1], function(i){
			value <- map.values[[r]][i, ]
			rhcollect(as.character(value$station.id), value)
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
		if(sum(!is.na(combined$tmax)) >= 300){
			rhcollect(reduce.key, combined)
		}
	}
)
job$input <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), par$loess
	), 
	type = "sequence"
)
job$output <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), paste(par$loess, "bystation", sep=".")
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72
)
job$mon.sec <- 5
job$jobname <- file.path(
	rh.datadir, par$dataset, "spatial", "a1950",
	par$family, paste("sp", par$span, sep=""), paste(par$loess, "bystation", sep=".")
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


## calculate the stl2 fitting for each station
parameter <- list(
	sw = "periodic",
	tw = 865,
	sd = 1,
	td = 1,
	fcw = 865,
	fcd = 1,
	ssw = 111,
	ssd = 1
) 
job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		map.values[[r]] <- map.values[[r]][order(map.values[[r]]$time), ]
		cycleSubIndices <- rep(c(1:12), ceiling(nrow(map.values[[r]])/12))[1:nrow(map.values[[r]])]

		if(sum(!is.na(map.values[[r]][[dataset]])) > 50 &
			all(by(map.values[[r]][[dataset]], list(cycleSubIndices), function(x) any(!is.na(x))))
		) {
			v.stl <- do.call("cbind", 
      	stl2(
          x = map.values[[r]][[dataset]], 
          t = map.values[[r]]$time, 
          n.p = 12, 
          s.window = sw, 
          s.degree = sd, 
          t.window = tw, 
          t.degree = td, 
          fc.window = c(fcw, ssw), 
          fc.degree = c(fcd, ssd), 
          inner = inner, 
          outer = outer
      	)[c("data","fc")]
    	)
			names(v.stl)[grep("fc.fc", names(v.stl))] <- c("fc.trend", "fc.seasonal")
			v.stl$fit <- v.stl$data.seasonal + v.stl$fc.trend + v.stl$fc.seasonal
			value <- map.values[[r]]
			value$stlfit<- v.stl$fit
			value <- value[!is.na(value[[dataset]]), ]
			rhcollect(sample(1:100, 1), value)
		}
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
		rhcollect(reduce.key, combined)
	}
)
job$parameters <- list(
  sw = parameter[["sw"]], 
  tw = parameter[["tw"]], 
  sd = parameter[["sd"]], 
  td = parameter[["td"]],
  fcw = parameter[["fcw"]], 
  fcd = parameter[["fcd"]], 
  ssw = parameter[["ssw"]], 
  ssd = parameter[["ssd"]], 
  inner = 10, 
  outer = 0, 
  dataset = par$dataset
)
job$setup <- expression(
  map = {
    library(lattice)
    library(yaImpute, lib.loc = lib.loc)
    library(stl2, lib.loc = lib.loc)
  }
)
job$input <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess02.bystation"), 
	type = "sequence"
) 
job$output <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess02.bystation.stl"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72
)
job$mon.sec <- 5
job$jobname <- file.path(rh.datadir, par$dataset, "spatial", "a1950", "loess02.bystation.stl")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


## group all stations to 10 groups
job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r){
		rhcollect(sample(1:10,1), map.values[[r]])
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
		rhcollect(reduce.key, combined)
	}
)
job$input <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), 
		paste(par$loess, "bystation", sep=".")
	), 
	type = "sequence"
)
job$output <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), 
		paste(par$loess, "bystation", "10pc", sep=".")
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 10
)
job$mon.sec <- 5
job$jobname <- file.path(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), 
		paste(par$loess, "bystation", "10pc", sep=".")
	)
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)
