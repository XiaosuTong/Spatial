source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$loess <- "loess02"
#par$family <- "gaussian"
par$span <- 0.025
par$family <- "symmetric"
#par$span <- 0.05
par$degree <- 2
par$multiple <- 0
par$file <- "bystation.all"

source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")

## calculate the stl2 fitting for each station
par$parameters <- list(
	sw = "periodic",
	tw = 109,
	sd = 1,
	td = 2,
	inner = 10,
	outer = 0
) 
job <- list()
if (par$multiple == 0) {
  job$map <- expression({
	  lapply(seq_along(map.values), function(r) {
		  map.values[[r]] <- map.values[[r]][order(map.values[[r]]$time), ]
		  Index <- which(!is.na(map.values[[r]]$tmax))
		  map.values[[r]]$fitted[Index] <- map.values[[r]]$tmax[Index]
		  map.values[[r]]$flag <- as.numeric(!is.na(map.values[[r]]$tmax))
		  v.stl <- stl2(
        x = map.values[[r]]$fitted, 
        t = map.values[[r]]$time, 
        n.p = 12, 
        s.window = parameters$sw, 
        s.degree = parameters$sd, 
        t.window = parameters$tw, 
        t.degree = parameters$td, 
        inner = parameters$inner, 
        outer = parameters$outer
      )$data
		  if (sum(!is.na(map.values[[r]]$tmax)) == 576) {
		  	rhcounter("stations", "_576_", 1)
        v.stl$type <- 1
		  } else {
		  	rhcounter("stations", "no_576_", 1)	 
        v.stl$type <- 0 	
		  }
			value <- cbind(map.values[[r]], subset(v.stl, select = -c(weights, sub.labels, raw)))
			rhcollect(map.keys[[r]], value)
	  })
  })
} else {
	job$map <- expression({
	  lapply(seq_along(map.values), function(r) {
		  map.values[[r]] <- map.values[[r]][order(map.values[[r]]$time), ]
		  Index <- which(!is.na(map.values[[r]]$tmax))
		  map.values[[r]]$fitted[Index] <- map.values[[r]]$tmax[Index]
		  map.values[[r]]$flag <- as.numeric(!is.na(map.values[[r]]$tmax))
		  v.stl <- do.call("cbind", 
      	stl2(
          x = map.values[[r]]$fitted, 
          t = map.values[[r]]$time, 
          n.p = 12, 
          s.window = parameters$sw, 
          s.degree = parameters$sd, 
          t.window = parameters$tw, 
          t.degree = parameters$td, 
          fc.window = c(parameters$fcw, parameters$ssw), 
          fc.degree = c(parameters$fcd, parameters$ssd), 
          inner = parameters$inner, 
          outer = parameters$outer
      	)[c("data","fc")]
      )
		  if (sum(!is.na(map.values[[r]]$tmax)) == 576) {
		  	rhcounter("stations", "_576_", 1)
		  	v.stl$type <- 1
			} else {
		  	rhcounter("stations", "no_576_", 1)
		  	v.stl$type <- 1
			}
			names(v.stl)[grep("fc.fc", names(v.stl))] <- c("fc.low", "fc.middle")
			value <- cbind(
				map.values[[r]], 
				subset(v.stl, select = -c(data.raw, data.weights, data.trend, data.remainder, data.sub.labels))
			)
			rhcollect(map.keys[[r]], value)
	  })
  })
}
job$parameters <- list(
  parameters = par$parameters,
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
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), 
		paste(par$loess, "bystation.all", sep=".")
	), 
	type = "sequence"
)
job$output <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", par$span, sep=""),
		paste(par$loess, par$file, "stl", par$multiple, sep=".")
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72
)
job$mon.sec <- 5
job$jobname <- file.path(
	rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", par$span, sep=""),
	paste(par$loess, par$file, "stl", par$multiple, sep=".")
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


load("~/Projects/Spatial/NCAR/RData/stations.100.RData")
job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r){
		value <- map.values[[r]]
		value$station.id <- map.keys[[r]]
		value$time1 <- c(rep(0:119, 4), 0:95)
    value$factor <- factor(
      rep(1:5, c(rep(120, 4), 96)),
      levels = c(5:1)
    )
		if (sum(!is.na(map.values[[r]]$tmax)) <= 17) { ## stations have less or equal 17 obs will be group 1, 100 stations
			rhcounter("stations", "_1_", 1)
			rhcollect(0, value)
		} else if (map.keys[[r]] %in% stations) { ## use the old 100 stations with most obs, will be group 2
			rhcounter("stations", "_2_", 1)
      rhcollect(1, value)
		}
	})
})
job$parameters <- list(
  stations = stations.tmax
)
job$input <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", par$span, sep=""),
		paste(par$loess, par$file, "stl", par$multiple, sep=".")
	), 
	type = "sequence"
)
job$output <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", par$span, sep=""),
		paste(par$loess, par$file, "stl", par$multiple, "compare", sep=".") 
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 72
)
job$mon.sec <- 5
job$jobname <- file.path(
	rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", par$span, sep=""),
	paste(par$loess, par$file, "stl", par$multiple, "compare", sep=".")
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


################################################
##             Ploting                        ##
################################################
library(lattice)
library(magrittr)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

rst <- rhread(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", par$span, sep=""),
		paste(par$loess, par$file, "stl", par$multiple, "compare", sep=".") 
	)
)

result <- do.call(rbind, lapply(rst, function(r){r[[2]][,c("remainder","type")]}))

trellis.device(
  device = postscript, 
  file = file.path(
    local.output, paste(par$dataset, "loess.stl.compare", "ps", sep = ".")
  ), 
  color=TRUE, 
  paper="legal"
)
	b <- qqmath( ~ remainder| factor(type)
		, data = result
		, distribution = qunif
		, aspect = 1
		, pch = 16
		, cex = 0.3
		, strip=strip.custom(factor.levels=c("w/", "w/o"))
		, xlab = "f-value"
		, ylab = "Remainder"
		, main = "Quantiles of Remainder from STL+"
		, panel = function(x, ...){
			panel.abline(h = 0, col = "black", lwd = 0.5)
			panel.abline(
				h = seq(-10, 10, by = 5), 
				v = seq(0, 1, by = 0.2), 
				col = "lightgray", lwd = 0.5
			)
			panel.qqmath(x,...)
		}
	)
	print(b)
dev.off()

ylab <- "Maximum Temperature (degrees centigrade)"

##  trend + seasonal component with original obs vs. month
trellis.device(
  device = postscript, 
  file = file.path(
    local.output, paste(par$dataset, "trend+seasonal.vs.month", "ps", sep = ".")
  ), 
  color=TRUE, 
  paper="legal"
)
lapply(rst, function(r) {
  xyplot((trend+seasonal) ~ time1 | factor
    , data = r[[2]]
    , layout = c(1,5)
    , strip = FALSE
    , aspect = "xy"
    , scale = list(
        y = list(
          relation = 'same', 
          tick.number=4, 
          alternating=TRUE
        ), 
        x = list(
          at = seq(0, 143, by=12), 
          relation = 'same'
        )
      )
    , key = list(
        text = list(label=c(
        	"imputed value",
          "raw obs",
          "trend+seasonal"
        )),
        lines = list(
          pch = 1, 
          cex = 0.5, 
          lwd = 1.5, 
          type = c("p","p","l"), 
          col = col[3:1]
        ), 
        columns = 3
      )
    , main = list(label = "Trend and Seasonal component vs. Month")
    , sub = paste(
        "Station", unique(r[[2]]$station.id)
      )
    , xlim = c(0,119)
    , ylab = list(label = ylab)
    , xlab = list(label = "Month")
    , panel = function(x, y, subscripts ,...) {
    	  data1 <- subset(r[[2]][subscripts, ], flag == 1)
    	  data2 <- subset(r[[2]][subscripts, ], flag == 0)
        panel.xyplot(
          x = sort(x), 
          y = y[order(x)], 
          type="b", col=col[1], pch=16,cex =0.5,lwd=1, ...
        )
        panel.xyplot(
          x = data1$time1, 
          y = data1$fitted,
          type="p", col=col[2], cex=0.5, ...
        )
        panel.xyplot(
        	x = data2$time1,
        	y = data2$fitted,
        	type="p", col=col[3], cex=0.5, ...
        )
        panel.abline(
          v=seq(0,119, by=12), color="lightgrey", lty=3, lwd=0.5
        )
      }
  )
})
dev.off()

for (k in 1:2) {

  if (k == 1) {
	  mainlab = "Remainder from Non-Imputed Value vs. Month"
  } else {
	  mainlab = "Remainder from Imputed Value vs. Month"
  }

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "remainder.vs.month.g", k, "ps", sep = ".")
    ), 
    color = TRUE,
    paper = "legal"
  )
  lapply(rst[(1:100)+(k-1)*100], function(r) {
    b <- xyplot(remainder ~ time,
      , data = r[[2]]
      , layout = c(1,1)
      , strip = FALSE
      , pch  = 16
      , cex  = 0.4
      , key = list(
          text = list(label=c(
            "remainder",
            "loess soomthing: span=0.15, degree=2",
            "loess soomthing: span=0.35, degree=1"
          )),
          lines = list(
            pch = 16, 
            cex = 0.5, 
            lwd = 1.5, 
            type = c("p","l","l"), 
            col = col[1:3]
          ), 
          columns = 3
        )
      , xlim = c(0, 576)
      , main = mainlab
      , ylab = list(label = ylab)
      , sub = paste(
          "Station", unique(r[[2]]$station.id)
        )    
      , xlab = list(label = "Month")
      , panel = function(x, y,...) {
          panel.abline(h=0, col="black", lty=1, lwd=0.5)
          panel.xyplot(x, y,...)
          panel.loess(x,y,degree=2,span=0.15, col=col[2], evaluation=200,...)
          panel.loess(x,y,degree=1,span=0.35, col=col[3], evaluation=200,...)
        }
    )
    print(b)
  })
  dev.off()

}

result <- do.call(rbind, lapply(rst, function(r){r[[2]]}))

for (k in 1:2) {

  if (k == 1) {
    mainlab = "Normal Quantiles of Remainder from Non-Imputed Value"
  } else {
    mainlab = "Normal Quantiles of Remainder from Imputed Value"
  }

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "QQ.remainder.g", k, "ps", sep = ".")
    ), 
    color = TRUE,
    paper = "legal"
  )
    a <- qqmath(~ remainder | factor(station.id)
      , data = result[(1:57600+(k-1)*57600), ]
      , distribution = qnorm
      , layout = c(5,3)
      , pch  = 16
      , aspect = 1
      , cex  = 0.3
      , main = mainlab
      , scale = list(y=list(relation="free"))
      , xlab = list(label = "Unit normal quantile")
      , ylab = list(label = "Remainder")
#     , prepanel = prepanel.qqmathline
      , panel = function(x, ...) {
          panel.abline(v=seq(-4,4,2), h=seq(-9,9,3), col="lightgrey", lty=1, lwd=0.5)
          panel.qqmathline(x, y=x)
          panel.qqmath(x, ...)
        }
    )
    print(a)
  dev.off()

}


## trend component and yearly mean vs.month
mav <- function(x,n=5){

  data.frame(mv=matrix(filter(x$mean,rep(1/n,n), sides=2), ncol=1), stringsAsFactors=FALSE)

}

dr <- result %>% ddply(
  .variable = c("station.id", "year"),
  .fun = summarise,
  mean = mean(fitted, na.rm=TRUE)
) %>% ddply(
  .variable = "station.id",
  .fun= mav
)

result<- result[order(result$station.id),]
result <- dr[rep(row.names(dr), each=12), "mv", drop=FALSE] %>% cbind(result)

for (k in 1:2) {

  if (k == 1) {
    mainlab = "Trend from Non-Imputed Value vs. Month"
  } else {
    mainlab = "Trend from Imputed Value vs. Month"
  }

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "trend.yealymeanmv.vs.month.g", k, "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  b <- xyplot( mv ~ time | factor(station.id) 
    , data = result[(1:57600+(k-1)*57600), ]
    , xlab = list(label = "Month", cex = 1.2)
    , ylab = list(label = ylab, cex = 1.2)
    , main = mainlab
    , xlim = c(0, 576)
    , pch = 16
    , layout = c(5,2)
#   , aspect="xy"
    , key=list(
        text = list(label=c("trend component","moving average of yearly mean")), 
        lines = list(pch=16, cex=0.7, lwd=1.5, type=c("l","p"), col=col[1:2]),
        columns=2
      )
    , scales = list(
        y = list(relation = 'free'), 
        x=list(at=seq(0, 576, by=120), relation = 'same')
      )
    , prepanel = function(x,y,subscripts,...){
        v <- result[(1:57600+(k-1)*57600), ][subscripts,]
        ylim <- range(v$mv, na.rm = TRUE)
        ans <- prepanel.default.xyplot(v$time, v$trend, ...)
        ans$ylim <- range(ans$ylim, ylim)
        ans
      }
    , panel = function(x, y, subscripts, ...) {
        panel.xyplot(
          x = x[seq(1, 576, by=12)],
          y = y[seq(1, 576, by=12)],
          type="p", col=col[2], cex = 0.5, ...
        )
        panel.xyplot(x, result[(1:57600+(k-1)*57600), ]$trend[subscripts], type="l", col=col[1], ...)
      }
  )
  print(b)
  dev.off()

}


## remainder vs year condtional on month in year
m.od <- c(
  "Jan","Feb","Mar","Apr","May","June",
  "July","Aug", "Sep", "Oct", "Nov", "Dec"
)

for (k in 1:2) {

  if (k == 1) {
    mainlab = "Remainder from Non-Imputed Value vs. Year"
  } else {
    mainlab = "Remainder from Imputed Value vs. Year"
  }

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "remainder.vs.year.g", k, "ps", sep = ".")
    ),    
    color = TRUE, 
    paper = "legal"
  )
  lapply(rst[(1:100)+(k-1)*100], function(r) {

    b <- xyplot( remainder ~ (as.numeric(year)-1949) | factor(month, levels= m.od)
      , data = r[[2]]
      , xlab = list(label = "Year", cex = 1.2)
      , ylab = list(label = ylab, cex = 1.2)
      , main = mainlab
      , sub = paste(
          "Station", unique(r[[2]]$station.id)
        ) 
      , pch = 16
      , cex = 0.5
      , key=list(
          text = list(label=c("remainder","loess smoothing: span=0.85, degree=1")), 
          lines = list(pch=16, cex=0.7, lwd=1.5, type=c("l","p"), col=col[1:2]),
          columns=2
        )
      , layout = c(12,1)
      , strip = TRUE
      , scales = list(
          y = list(relation = 'same', alternating=TRUE), 
          x = list(tick.number=10, relation='same')
        )
      , panel = function(x,y,...){
          panel.abline(h=0, color="black", lty=1, lwd=0.5)
          panel.xyplot(x,y,...)
          panel.loess(x,y,span=0.85, degree=1, col=col[2],...)
        }
    )
    print(b)

  })
  dev.off()

}
for (k in 1:2) {

  if (k == 1) {
    mainlab = "Remainder from Non-Imputed Value vs. Year"
  } else {
    mainlab = "Remainder from Imputed Value vs. Year"
  }

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "remainder.vs.year2.g", k, "ps", sep = ".")
    ),    
    color = TRUE, 
    paper = "legal"
  )
  lapply(rst[(1:100)+(k-1)*100], function(r) {

    b <- xyplot( remainder ~ (as.numeric(year)-1949) | factor(month, levels=m.od)
      , data = r[[2]]
      , xlab = list(label = "Year", cex = 1.2)
      , ylab = list(label = ylab, cex = 1.2)
      , main = mainlab
      , sub = paste(
          "Station", unique(r[[2]]$station.id)
        ) 
      , pch = 16
      , cex = 0.5
      , type = "b"
      , layout = c(2,6)
      , strip = TRUE
      , scales = list(
          y = list(relation = 'same', alternating=TRUE), 
          x = list(tick.number=10, relation='same')
        )
      , panel = function(x,y,...){
          panel.abline(h=0, color="black", lty=1, lwd=0.5)
          panel.xyplot(x,y,...)
        }
    )
    print(b)

  })
  dev.off()

}


## Normal quantiles of remainder conditional on month in year
for (k in 1:2) {

  if (k == 1) {
    mainlab = "Normal quantiles Remainder from Non-Imputed Value"
  } else {
    mainlab = "Normal quantiles Remainder from Imputed Value"
  }

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "QQ.remainder.month.g", k, "ps", sep = ".")
    ),     
    color = TRUE, 
    paper = "legal"
  )
  lapply(rst[(1:100)+(k-1)*100], function(r) {
  
    a <- qqmath(~ remainder | month,
      , data = r[[2]]      
      , distribution = qnorm
      , aspect = 1
      , pch = 16
      , cex = 0.4
      , layout = c(4,3)
      , main = mainlab
      , sub = paste(
          "Station", unique(r[[2]]$station.id)
        ) 
      , xlab = list(label="Unit normal quantile")
      , ylab = list(label = "Remainder")
      , prepanel = prepanel.qqmathline,
      , panel = function(x, y,...) {
         panel.grid(lty=3, lwd=0.5, col="black",...)
         panel.qqmathline(x, y=x)
         panel.qqmath(x, y,...)
       }
    )
    print(a)

  })
  dev.off()

}

for (k in 1:2) {

  if (k == 1) {
    mainlab = "Remainder from Non-Imputed Value"
  } else {
    mainlab = "Remainder from Imputed Value"
  }

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "Acf.remainder.g", k, "ps", sep = ".")
    ),     
    color = TRUE, 
    paper = "legal"
  )
  lapply(rst[(1:100)+(k-1)*100], function(r) {

    ACF <- ddply(
      .data = r[[2]],
      .variables = "station.id",
      .fun = summarise,
      correlation = c(acf(remainder, plot=FALSE)$acf),
      lag = c(acf(remainder, plot=FALSE)$lag) 
    )

    b <- xyplot( correlation ~ lag
      , data = subset(ACF, lag!=0)
      , xlab = "Lag"
      , ylab = "ACF"      
      , main = mainlab
      , sub = paste(
          "Station", unique(ACF$station.id)
        ) 
      , type = "h"
      , panel = function(x,y,...) {
          panel.abline(h=0)
          panel.xyplot(x,y,...)
        }
    )
    print(b)

  })
  dev.off()

}