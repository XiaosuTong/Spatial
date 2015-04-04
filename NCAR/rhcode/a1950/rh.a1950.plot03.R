#############################################
##spatial loess fit comparasion after 1950 ##
#############################################

## compare the loess fit without elev with fit without elev
## so compare loess02 with loess04
library(maps)
library(lattice)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$loess1 <- "loess02.bystation.all"
par$loess2 <- "loess04.bystation.all"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

Qsample <- function(num) {

  x <- num[!is.na(num)]
  x <- sort(x)
  len <- length(x)
  idx <- round(seq(1, len, length.out = 800))
  f.value <- (idx - 0.5) / len
  value <- data.frame(
  	resid = x[idx],
  	fv = f.value
  )

  value

}

job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		file <- Sys.getenv("mapred.input.file")
		k <- strsplit(tail(strsplit(file, "/")[[1]],3)[1:2], "[.]")
		key <- c(k[[1]][2], k[[2]][1])
		if ("05" %in% key[1]) {
			key[1] <- "0.05"
		} else {
			key[1] <- "0.025"
		}
		if ("loess04" %in% key[2]) {
			key[2] <- "w/ elev"
		} else {
			key[2] <- "w/o elev"
		}
		resid <- with(map.values[[r]], tmax-fitted) 
		value <- resid[!is.na(resid)]
		rhcollect(key, value)
	})
})
job$reduce <- expression(
	pre = {
		combine <- vector()
	},
	reduce = {
		combine <- c(combine, unlist(reduce.values))
	},
	post = {
		value <- Qsample(combine)
		value$elev <- reduce.key[1]
		value$span <- reduce.key[2]
		rhcollect(reduce.key, value)
	}
)
job$input <- rhfmt(
	c(
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
	    "symmetric", "sp0.05", par$loess1 
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
		  "symmetric", "sp0.05", par$loess2
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
	    "symmetric", "sp0.025", par$loess1 
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
		  "symmetric", "sp0.025", par$loess2
	  )
	),	 
	type = "sequence"
)
job$parameters <- list(
	Qsample = Qsample
)
job$output <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "a1950", "compare.elev"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 4
)
job$mon.sec <- 5
job$jobname <- 	file.path(
	rh.datadir, par$dataset, "spatial", "a1950", "compare.elev"
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		file <- Sys.getenv("mapred.input.file")
		k <- strsplit(tail(strsplit(file, "/")[[1]],3)[1:2], "[.]")
		key <- c(k[[1]][2], k[[2]][1])
		if ("05" %in% key[1]) {
			key[1] <- "0.05"
		} else {
			key[1] <- "0.025"
		}
		if ("loess04" %in% key[2]) {
			key[2] <- "w/ elev"
		} else {
			key[2] <- "w/o elev"
		}
		resid <- with(map.values[[r]], tmax-fitted) 
		value <- data.frame(matrix(summary(resid[!is.na(resid)]), nrow = 1))
		names(value) <- c("min", "25Q", "median", "mean", "75Q", "max")
		rhcollect(key, value)
	})
})
job$reduce <- expression(
	pre = {
		combine <- vector()
	},
	reduce = {
		combine <- rbind(combine, do.call(rbind, reduce.values))
	},
	post = {
		combine$elev <- reduce.key[1]
		combine$span <- reduce.key[2]
		rhcollect(reduce.key, combine)
	}
)
job$input <- rhfmt(
	c(
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
	    "symmetric", "sp0.05", "loess04" 
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
		  "symmetric", "sp0.05", "loess02"
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
	    "symmetric", "sp0.025", "loess04" 
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
		  "symmetric", "sp0.025", "loess02"
	  )
	),	 
	type = "sequence"
)
job$output <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", "compare.elev.summary"
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 4
)
job$mon.sec <- 5
job$jobname <- 	file.path(
	rh.datadir, par$dataset, "spatial", "a1950", "compare.elev.summary"
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)

compare <- function(file = "compare") {

	rst <- rhread(	
	  file.path(rh.datadir, par$dataset, "spatial", "a1950", "compare.elev")
	)
	result <- do.call(rbind, lapply(rst, "[[", 2))

  trellis.device(
	  device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "residual.compare.elev", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )

  a <- xyplot( resid ~ fv | factor(elev)*factor(span)
    , data = result
    , pch = 16
    , cex = 0.3
    , layout = c(4,1)
    , xlab = "f-value"
    , ylab = "Residuals"
    , main = "Quantiles of Residuals"
    , panel = function(x,y,...) {
    	  panel.abline(h=0, col="black", lwd = 0.5)
    	  panel.abline(h=seq(-20,20,by=5), v=seq(0,1,by=0.25), col="lightgrey", lwd=0.5)
    	  panel.xyplot(x,y,...)
    }
  )
  print(a)	
  
  dev.off()

  rst <- rhread(	
	  file.path(rh.datadir, par$dataset, "spatial", "a1950", "compare.elev.summary")
	)
	result <- do.call(rbind, lapply(rst, "[[", 2))

  simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
  }

	trellis.device(
	  device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "residual.compare.summary.elev", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )

	for (k in names(result)[1:6]) {
    a <- qqmath( ~ get(k) | factor(elev)*factor(span)
      , data = result
      , distribution = qunif
      , cex = 0.3
      , pch = 16
      , layout = c(4,1)
      , xlab = "f-value"
      , ylab = "Residuals"
      , main = paste("Quantiles of", simpleCap(k), "of Residuals")
      , panel = function(x,...) {
    	    panel.abline(h=0, col="black", lwd = 0.5)
    	    panel.abline(h=seq(-20,20,by=5), v=seq(0,1,by=0.25), col="lightgrey", lwd=0.5)
    	    panel.qqmath(x,...)
      }
    )
    print(a)
  }	
  
  dev.off()
}	


#######################################################
##         compare the Remainders                    ##
#######################################################

file <- "bystation.all.stl.0.bymonth"

Qsample <- function(data) {

  x <- data$remainder[!is.na(data$remainder)]
  x <- sort(x)
  len <- length(x)
  idx <- round(seq(1, len, length.out = 800))
  f.value <- (idx - 0.5) / len
  value <- data.frame(
  	resid = x[idx],
  	fv = f.value
  )

  value

}

job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		file <- Sys.getenv("mapred.input.file")
		k <- strsplit(tail(strsplit(file, "/")[[1]],3)[1:2], "[.]")
		key <- c(k[[1]][2], k[[2]][1])
		if ("05" %in% key[1]) {
			key[1] <- "0.05"
		} else {
			key[1] <- "0.025"
		}
		if ("loess04" %in% key[2]) {
			key[2] <- "w/ elev"
		} else {
			key[2] <- "w/o elev"
		}
		value <- map.values[[r]][, c("flag","remainder")]
		rhcollect(key, value)
	})
})
job$reduce <- expression(
	pre = {
		combine <- vector()
	},
	reduce = {
		combine <- rbind(combine, do.call(rbind, reduce.values))
	},
	post = {
		value <- ddply(
			.data = combine,
			.variable = "flag",
			.fun = Qsample
		)
		value$elev <- reduce.key[1]
		value$span <- reduce.key[2]
		rhcollect(reduce.key, value)
	}
)
job$setup <- expression(
	reduce = {
		library(plyr, lib.loc = lib.loc)
	}
)
job$input <- rhfmt(
	c(
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
	    "symmetric", "sp0.05", paste("loess02", file, sep=".")
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
		  "symmetric", "sp0.05", paste("loess04", file, sep=".")
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
	    "symmetric", "sp0.025", paste("loess02", file, sep=".") 
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
		  "symmetric", "sp0.025", paste("loess04", file, sep=".")
	  )
	),	 
	type = "sequence"
)
job$parameters <- list(
	Qsample = Qsample
)
job$output <- rhfmt(
	file.path(rh.datadir, par$dataset, "spatial", "a1950", "compare.imputed"), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 4
)
job$mon.sec <- 5
job$jobname <- 	file.path(
	rh.datadir, par$dataset, "spatial", "a1950", "compare.imputed"
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		file <- Sys.getenv("mapred.input.file")
		k <- strsplit(tail(strsplit(file, "/")[[1]],3)[1:2], "[.]")
		key <- c(k[[1]][2], k[[2]][1])
		if ("05" %in% key[1]) {
			key[1] <- "0.05"
		} else {
			key[1] <- "0.025"
		}
		if ("loess04" %in% key[2]) {
			key[2] <- "w/ elev"
		} else {
			key[2] <- "w/o elev"
		}
		value <- map.values[[r]][, c("flag","remainder")]
		value <- ddply(
			.data = value,
			.variable = "flag",
			.fun = summarise,
			min = summary(remainder)[1],
			Q25 = summary(remainder)[2],
			median = summary(remainder)[3],
			mean = summary(remainder)[4],
			Q75 = summary(remainder)[5],
			max = summary(remainder)[6]
		)
		rhcollect(key, value)
	})
})
job$reduce <- expression(
	pre = {
		combine <- vector()
	},
	reduce = {
		combine <- rbind(combine, do.call(rbind, reduce.values))
	},
	post = {
		combine$elev <- reduce.key[1]
		combine$span <- reduce.key[2]
		rhcollect(reduce.key, combine)
	}
)
job$input <- rhfmt(
	c(
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
	    "symmetric", "sp0.05", paste("loess02", file, sep=".")
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
		  "symmetric", "sp0.05", paste("loess04", file, sep=".")
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
	    "symmetric", "sp0.025", paste("loess02", file, sep=".") 
	  ),
	  file.path(
		  rh.datadir, par$dataset, "spatial", "a1950", 
		  "symmetric", "sp0.025", paste("loess04", file, sep=".")
	  )
	),	 
	type = "sequence"
)
job$setup <- expression(
	map = {
		library(plyr, lib.loc = lib.loc)
	},
	reduce = {
		library(plyr, lib.loc = lib.loc)
	}
)
job$output <- rhfmt(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", "compare.impute.summary"
	), 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 4
)
job$mon.sec <- 5
job$jobname <- 	file.path(
	rh.datadir, par$dataset, "spatial", "a1950", "compare.impute.summary"
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)

compare.impute <- function(file = "compare") {

	rst <- rhread(	
	  file.path(rh.datadir, par$dataset, "spatial", "a1950", "compare.imputed")
	)
	result <- do.call(rbind, lapply(rst, "[[", 2))

  trellis.device(
	  device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "remainder.compare.imputed", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )

  a <- xyplot( resid ~ fv | factor(elev)*factor(span)
    , data = result
    , pch = 16
    , cex = 0.4
    , group = flag
    , key=list(
          text = list(label=c("w/ imputed value","w/o imputed value")), 
          lines = list(pch=16, cex=0.7, lwd=1.5, type=c("p","p"), col=col[1:2]),
          columns=2
      ) 
    , layout = c(4,1)
    , xlab = "f-value"
    , ylab = "Remainder"
    , main = "Quantiles of Remainder"
    , panel = function(x,y,...) {
    	  panel.abline(h=0, col="black", lwd = 0.5)
    	  panel.abline(h=seq(-20,20,by=5), v=seq(0,1,by=0.25), col="lightgrey", lwd=0.5)
    	  panel.xyplot(x,y,...)
    }
  )
  print(a)	
  
  dev.off()

  rst <- rhread(	
	  file.path(rh.datadir, par$dataset, "spatial", "a1950", "compare.impute.summary")
	)
	result <- do.call(rbind, lapply(rst, "[[", 2))

  simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
  }

	trellis.device(
	  device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, "remainder.compare.summary.imputed", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )

	for (k in names(result)[2:7]) {
    a <- qqmath( ~ get(k) | factor(elev)*factor(span)
      , data = result
      , distribution = qunif
      , cex = 0.3
      , group = flag
      , pch = 16
      , key=list(
          text = list(label=c("w/ imputed value","w/o imputed value")), 
          lines = list(pch=16, cex=0.7, lwd=1.5, type=c("p","p"), col=col[1:2]),
          columns=2
        ) 
      , layout = c(4,1)
      , xlab = "f-value"
      , ylab = "Remainder"
      , main = paste("Quantiles of", simpleCap(k), "of Remainder")
      , panel = function(x,...) {
    	    panel.abline(h=0, col="black", lwd = 0.5)
    	    panel.abline(h=seq(-20,20,by=5), v=seq(0,1,by=0.25), col="lightgrey", lwd=0.5)
    	    panel.qqmath(x,...)
      }
    )
    print(a)
  }	
  
  dev.off()
}	
