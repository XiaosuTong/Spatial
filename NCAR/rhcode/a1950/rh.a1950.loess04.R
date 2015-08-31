########################################################
## spatial loess fit at all stations after 1950 using
## longitude, latitude, and elevation.
##
##  - my.loess02: 
##      is the loess function that calculate kd-tree 
##      nodes, and all necessary information for interpolation
##
##  - loess02:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess02"
##      is key-value pairs that key is c(year, month), value is loess 
##      fit with family, degree, and span w/o elevation at that month.
##
##  - loess04:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess04"
##      is key-value pairs that key is c(year, month), value is loess 
##      fit with family, degree, and span w/ elevation at that month.
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
##  - loess02/loess04.bystation.all:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess02.bystation.all"
##      is key-value pairs that key is station.id, value is data.frame
##      for the station. All stations included.
##
##  - loess02/loess04.bystation.all.10pc:
##      "/ln/tongx/Spatial/tmp/tmax/spatial/a1950/family/span/loess02.bystation.all.10pc"
##      merge loess02.bystation.all to 10 key-value pairs, key is meaningless
##      faster to combine all 10 pieces for ploting.
##
###################################################

par$loess <- "loess04"
par$span <- 0.035
par$family <- "symmetric"
par$degree <- 2
par$all <- TRUE
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")

job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
    v <- map.values[[r]]
    v$elev2 <- log2(v$elev + 128)
    lo.fit <- my.loess2( resp ~ lon + lat + elev2, 
      data       = v, 
      degree     = par$degree, 
      span       = par$span,
      parametric = "elev2",
      family     = par$family,
      control = loess.control(iterations = 10)
    )
  	fit <- my.predict.loess(
  		object = lo.fit, 
      newdata = data.frame(
      	lon = v$lon, 
      	lat = v$lat,
      	elev2 = v$elev2
      )
  	)
  	v$fitted <- fit
  	v$station.id <- as.character(v$station.id)
  	rhcollect(map.keys[[r]], v)
  })
})
job$setup <- expression(
  map = {
    system("chmod 777 myloess2.so")
    dyn.load("myloess2.so")
    library(maps)
  }
)
job$shared <- c(
	file.path(file.path(rh.root, par$dataset, "shareRLib", "myloess2.so"))
)
job$parameters <- list(
	par = par,
	my.loess2 = my.loess2,
	my.simple2 = my.simple2,
	my.predict.loess = my.predict.loess,
	my.predLoess = my.predLoess
)
job$input <- rhfmt(
  file.path(rh.root, par$dataset, "a1950", "bymonth"), 
  type = "sequence"
)
job$output <- rhfmt(
	file.path(
		rh.root, par$dataset, "a1950", "spatial",
		par$family, paste("sp", par$span, sep=""), par$loess
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72
)
job$mon.sec <- 10
job$jobname <- 	file.path(
  rh.root, par$dataset, "a1950", "spatial",
  par$family, paste("sp", par$span, sep=""), par$loess
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)

## change the key from c(year, month) to station.id
job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		m <- match(map.keys[[r]][2], month.abb)
		map.values[[r]]$time <- (as.numeric(map.keys[[r]][1]) - 1950)*12 + m
		map.values[[r]]$year <- map.keys[[r]][1]
		map.values[[r]]$month <- map.keys[[r]][2]
		lapply(1:dim(map.values[[r]])[1], function(i){
			value <- map.values[[r]][i, ]
			rhcollect(as.character(value$station.id), value)
		})
	})
})
if (!par$all) {
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
  file <- "bystation"
} else {
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
  file <- "bystation.all"
}
job$input <- rhfmt(
  file.path(
    rh.root, par$dataset, "a1950", "spatial", 
    par$family, paste("sp", par$span, sep=""), par$loess
  ), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(
    rh.root, par$dataset, "a1950", "spatial", 
    par$family, paste("sp", par$span, sep=""), paste(par$loess, file, sep=".")
  ), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 72
)
job$mon.sec <- 5
job$jobname <- file.path(
  rh.root, par$dataset, "a1950", "spatial", 
  par$family, paste("sp", par$span, sep=""), paste(par$loess, file, sep=".")
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)
