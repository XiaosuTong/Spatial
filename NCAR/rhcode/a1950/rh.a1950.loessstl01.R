source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$loess <- "loess04"
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
		  map.values[[r]]$flag <- as.numeric(!is.na(map.values[[r]]$tmax)) # flag is 1: fitted is obs; flag is 0: fitted is imputed
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
      .data = r[[2]][!is.na(r[[2]]$remainder), ],
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


#######################################################
##         Remainder to feed spatial loess           ##
#######################################################

job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
    lapply(1:nrow(map.values[[r]]), function(k) {
      key <- c(map.values[[r]]$year[k], map.values[[r]]$month[k])
      value <- map.values[[r]][k, c("station.id", "elev", "lon", "lat", "fitted", "flag", "remainder")]
      rhcollect(key, value)
    })
  })
})
job$reduce <- expression(
  pre = {
    combine <- data.frame()
  },
  reduce = {
    combine <- rbind(combine, do.call(rbind, reduce.values))
  },
  post = {
    rhcollect(reduce.key, combine)
  }
)
job$combiner <- TRUE
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
    paste(par$loess, par$file, "stl", par$multiple, "bymonth", sep=".")
  ), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 72
)
job$mon.sec <- 5
job$jobname <- file.path(
  rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", par$span, sep=""),
  paste(par$loess, par$file, "stl", par$multiple, "bymonth", sep=".")
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


rst <- rhread(
  file.path(
    rh.datadir, par$dataset, "spatial", "a1950", par$family, paste("sp", par$span, sep=""),
    paste(par$loess, par$file, "stl", par$multiple, "bymonth", sep=".")
  )
)

library(magrittr)

Dist.month <- function(rst = rst) {

  tmp <- do.call(rbind, lapply(rst,"[[",1)) %>%
    data.frame(stringsAsFactors=FALSE) 
  names(tmp) <- c("year", "month")
  data <- tmp[rep(1:nrow(tmp), each = 7738),] %>% 
    cbind(do.call("rbind", lapply(rst, "[[", 2)))

  Qrst <- ddply(
    .data = data,
    .variable = c("year", "month", "flag"),
    .fun = function(r) {
      r <- r[!is.na(r$remainder),]
      a <- sort(r$remainder)
      idx <- c(1:19, round(seq(20, nrow(r)-20, length.out = 160), (nrow(r)-20):nrow(r)))
      f.value <- (idx - 0.5) / length(a)
      value <- data.frame(
        remainder = a[idx], 
        fv = f.value
      )
    }
  )

  Qrst$month <- factor(
  Qrst$month, 
    levels = c(
      "Jan","Feb","Mar","Apr","May","June",
      "July","Aug", "Sep", "Oct", "Nov", "Dec"
    )
  )

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "Dist.remainder.month", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )

  for (k in 1950:1997) {
    a <- xyplot(remainder ~ fv| factor(flag, label=c("w/", "w/o"))*month*factor(year)
      , data = Qrst
      , subset = year == paste(k)
      , layout = c(6,2)
      , pch  = 16
      , cex  = 0.3
      , xlab = list(label = "f-value")
      , ylab = list(label = "Remainder")
      , main = "Quantiles of Remainder"
      , panel = function(x, y, subscripts, ...) {
          panel.abline(h=seq(-15,10,by=5), v=seq(0,1,by=0.2), col="lightgray", lwd=0.5)
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
        }
    )
    print(a)

  }

  dev.off()

  Qrange <- Qrst %>% ddply(
    .variable = c("year", "month", "flag"),
    .fun = summarise,
    time = (as.numeric(unique(year))-1950)*12 + match(substr(unique(month),1,3), month.abb),
    maxi = max(remainder),
    mini = min(remainder)
  )

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "Range.remainder.month", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  
  for (k in c("maxi", "mini")) {

    if (k == "maxi") {
      mainlab <- "Maximum of Remainder over Each Month"
      ord <- subset(Qrange, flag==0)$time[order(subset(Qrange, flag == 0)$maxi)]
    } else {
      mainlab <- "Minimum of Remainder over Each Month"
      ord <- subset(Qrange, flag==0)$time[order(subset(Qrange, flag == 0)$mini)]
    }

    a <- dotplot( get(k) ~ factor(time, levels=ord)
      , data = Qrange
      , pch = 1
      , group = flag
      , cex = 0.5
      , key = list(
          text = list(label=c(
            "w/ imputed value",
            "w/o imputed value"        
          )),
          lines = list(
            pch = 1, 
            cex = 0.5, 
            type = c("p","p"), 
            col = col[1:2]
          ), 
          columns = 2
        )
      , scale = list(x= list(draw = FALSE))
      , xlab = "Month"
      , ylab = "Remainder"
      , main = mainlab
      , levels.fos = FALSE
    )
    print(a)

  }
  dev.off()

}

loess.comp.plot <- function(loess="loess04" comp = "remainder", f = "symmetric", s=0.025) {

  file <- paste(loess, par$file, "stl", par$multiple, "bymonth", sep = ".")

  rst <- rhread(
    file.path(
      rh.datadir, par$dataset, "spatial", "a1950", f, paste("sp", s, sep=""), file
    )
  )

  tmp <- do.call(rbind, lapply(rst,"[[",1)) %>%
    data.frame(stringsAsFactors=FALSE) 
  tmp$X2 <- tmp$X2 %>% substr(1,3) %>% match(month.abb)
  mod <- tmp %>% with(tmp[order(X1, X2),]) %>% row.names() %>% as.numeric()

  ## a simple function to Capitalize the first character of a string vector
  simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
  }

  if (comp == "remainder") {
    ylab <- simpleCap(comp)
  } else {
    ylab <- paste(simpleCap(comp), "component")
  } 

  mainlab <- paste(simpleCap(comp), "vs. Latitude")

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.lat.lon", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ lat | equal.count(lon, 20, overlap=0) 
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Longitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
           relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Latitude"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.lat.elev", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ lat | equal.count(elev, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Elevation", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Latitude"
      , ylab = ylab 
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()

  mainlab <- paste(simpleCap(comp), "vs. Longitude")

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.lon.lat", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ lon | equal.count(lat, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Latitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Longitude"
      , ylab = ylab    
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.lon.elev", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ lon | equal.count(elev, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Elevation", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Longitude"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()
 
  mainlab <- paste(simpleCap(comp), "vs. Elevation")

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.elev.lat", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ elev2 | equal.count(lat, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Latitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Log (Elevation + 128) (log base 2 meter)"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.elev.lon", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ elev2 | equal.count(lon, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Longitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Log (Elevation + 128) (log base 2 meter)"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()

}
