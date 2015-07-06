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
par$span <- seq(0.01,0.05,0.001)
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


job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
    file <- Sys.getenv("mapred.input.file")
    drop <- substr(tail(strsplit(file, "/")[[1]],3)[2], 6, 6)
    span <- substr(tail(strsplit(file, "/")[[1]],3)[2], 9, 13)
    residual <- mean(with(map.values[[r]], (tmax - fitted)^2), na.rm=TRUE)
    value <- data.frame(
      year = map.keys[[r]][1],
      month = map.keys[[r]][2],
      drop = as.numeric(drop),
      span = as.numeric(span),
      resid = residual,
      stringsAsFactors = FALSE
    )
    rhcollect(1, value)
  })
})
job$reduce <- expression(
  pre= {
    combine <- data.frame()
  },
  reduce = {
    combine <- rbind(combine, do.call(rbind, reduce.values))
  },
  post = {
    rhcollect(reduce.key, combine)
  }
)
job$input <- rhfmt(
  file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, "spatial.cv"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, "spatial.cvComp"), 
  type = "sequence"
)
job$mapred <- list(mapred.reduce.tasks = 1)
job$mon.sec <- 10
job$jobname <- file.path("spatial.cvComp")  
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)



rst <- rhread(
  file.path(rh.datadir, par$dataset, "spatial", "a1950", par$family, "spatial.cvComp")
)[[1]][[2]]

monOrder <- c(
    "Jan","Feb","Mar","Apr","May","June",
    "July","Aug", "Sep", "Oct", "Nov", "Dec"
)

rst$month <- factor(rst$month, levels = monOrder)

trellis.device(
  device = postscript, 
  file = file.path(local.output, "crossValid.ps"),
  color = TRUE, 
  paper = "legal"
)
xyplot(resid ~ as.numeric(as.character(span)) | month*year
  , data = rst
  , subset = resid < 5
  , group = drop
  , cex = 0.5
  , key = list(
      text = list(label=c(
        "qudratic elevation",
        "linear elevation"
      )),
      lines = list(
        pch=c(".","."), 
        cex=4, 
        type=c("p","p"), 
        col=col[1:2]
      ), 
      columns = 2
    )
  , xlab = "Span"
  , ylab = "Mean squared error"
  , pch = 16
  , scale = list(y=list(relation="free"))
  , layout = c(4,3)
)
dev.off()

trellis.device(
  device = postscript, 
  file = file.path(local.output, "crossValid.dist.ps"),
  color = TRUE, 
  paper = "legal"
)
qq(drop ~ resid | factor(as.numeric(as.character(span)))
  , data = rst
  , subset = resid < 5
  , cex = 0.5
  , pch = 16
  , aspect = 1
  , layout = c(4,2)
  , scale = list(relation="free")
  , xlab = "Elevation degree = 2"
  , ylab = "Elevation degree = 1"
  , main = "Quantiles-quantiles of Mean Squared Error"
  , panel = function(x,y, ...){
      panel.abline(a=c(0,1), col="red", lwd=1)
      panel.qq(x,y,...)
  }
)
dev.off()

result <- ddply(
  .data = rst,
  .variable = c("drop", "span"),
  .fun = summarise,
  resid = mean(resid)
)

trellis.device(
  device = postscript, 
  file = file.path(local.output, "crossValidMean.ps"),
  color = TRUE, 
  paper = "letter"
)
xyplot(resid ~ as.numeric(as.character(span))
  , data = result
  , group = drop
  , cex = 1
  , key = list(
      text = list(label=c(
        "qudratic elevation",
        "linear elevation"
      )),
      lines = list(
        pch=c(1,1), 
        cex=1, 
        type=c("p","p"), 
        col=col[1:2]
      ), 
      columns = 2
    )
  , xlab = "Span"
  , ylab = "MSE"
  , main = "Mean Squared Error vs. Span"
  , pch = 1
  , type = "b"
  , scale = list(x=list(at=seq(0.01,0.05, by=0.004)))
  , panel = function(x,y, ...) {
      panel.xyplot(x,y,...)
      panel.abline(h=min(y), v=x[which.min(y)], lty=1, lwd=0.5, col="black")
  }
)
dev.off()


