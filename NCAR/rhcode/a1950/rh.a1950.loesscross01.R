source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$family <- "symmetric"
par$all <- TRUE
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")
source("~/Projects/Spatial/NCAR/myloess/kdtree.R")

par$degree <- 2
par$span <- seq(0.1,0.5,0.005)
par$drop <- c(TRUE, FALSE) # drop for the qudratic term of elevation

sampleStation <- function(df, seed) {
  
  data <- subset(df, !is.na(tmax))

  row.names(data) <- 1:nrow(data)

  rst <- cppkdtree(as.matrix(data[,c("lat", "lon", "elev")]), 200)
  
  idx <- ddply(
	  .data = rst,
	  .variable = "leaf",
	  .fun = function(r) {
		  set.seed(seed)
		  data.frame(rowID = sample(r$idx, 3), rep = c(1,2,3))
	  }
  )

  predStation <- as.character(data$station.id[idx$rowID])

  predStation

}

for(a in 1:length(par$drop)) {
	for(b in 1:length(par$span)) {

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
    	)[, c("station.id", "elev", "lon", "lat", "tmax")]
    	v$elev2 <- log2(v$elev + 128)
      seed <- set.seed(map.keys[[r]])
      
      stations <- sampleStation(v, seed)

    	if (par$drop[a]) {
        lo.fit <- my.loess2( tmax ~ lon + lat + elev2, 
          data       = subset(v, !(station.id %in% stations)), 
          degree     = par$degree, 
          span       = par$span[b],
          parametric = "elev2",
          drop.square = "elev2",
          family     = par$family
        )
      } else {
        lo.fit <- my.loess2( tmax ~ lon + lat + elev2, 
          data       = subset(v, !(station.id %in% stations)), 
          degree     = par$degree, 
          span       = par$span[b],
          parametric = "elev2",
          family     = par$family
        )
      }
      
      value <- subset(v, station.id %in% stations)

    	fit <- my.predict.loess(
    		object = lo.fit, 
        newdata = data.frame(
        	lon = value$lon, 
        	lat = value$lat,
        	elev2 = value$elev2
        )
    	)
    	value$fitted <- fit
    	value$station.id <- as.character(value$station.id)
    	rhcollect(c(y, m), value)
      })
    })
    job$setup <- expression(
    	map = {
        library(plyr, lib.loc = lib.loc)
    		system("chmod 777 myloess2.so")
    		dyn.load("myloess2.so")
    		system("chmod 777 cppkdtree.so")
        dyn.load("cppkdtree.so")
    	  load(paste(par$dataset, "a1950", "RData", sep="."))
    	}
    )
    job$shared <- c(
    	file.path(
    		rh.datadir, par$dataset, "shareRLib", "myloess2.so"
    	),
    	file.path(
    		rh.datadir, par$dataset, "shareRLib", "cppkdtree.so"
    	),
    	file.path(
    		rh.datadir, par$dataset, "a1950", "Rdata", 
    		paste(par$dataset, "a1950", "RData", sep=".")
    	)
    )
    job$parameters <- list(
    	a = a,
    	b = b,
    	par = par,
    	sampleStation = sampleStation,
    	my.loess2 = my.loess2,
    	my.simple2 = my.simple2,
    	my.predict.loess = my.predict.loess,
    	my.predLoess = my.predLoess,
      cppkdtree = cppkdtree
    )
    job$input <- c(par$N, 100) 
    job$output <- rhfmt(
    	file.path(
    		rh.datadir, par$dataset, "spatial", "a1950", 
    		par$family, "spatial.cv", paste("elevd", as.numeric(par$drop[a]), "sp", par$span[b], sep="")
    	), 
    	type = "sequence"
    )
    job$mapred <- list(
    	mapred.reduce.tasks = 72, 
    	rhipe_reduce_buff_size = 10000
    )
    job$mon.sec <- 10
    job$jobname <- file.path( 
    	paste("elevd", as.numeric(par$drop[a]), "sp", par$span[b], sep="")
    )
    job$readback <- FALSE
    job.mr <- do.call("rhwatch", job)

  }
}