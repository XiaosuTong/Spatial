#############################################################
##Dot plot for each station to see the where are missing
#############################################################
source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
if(par$dataset == "precip") {
	ylab <- "Precipitation (millimeters)"
}else if(par$dataset == "tmax") {
	ylab <- "Maximum Temperature (degree centigrade)"
}else {
	ylab <- "Minimum Temperature (degree centigrade)"
}
job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
	v <- map.values[[r]]
	k <- map.values[[r]]$partion[1]
	v$obs <- !is.na(v[, par$dataset])
	month <- v$month
	levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
	month <- as.numeric(factor(month, levels=c(1:12)))
	date <- paste(v$year, month, "01", sep="-")
	v$date <- as.POSIXct(
		strptime(date, format = "%Y-%m-%d"), 
		format = '%Y%m%d', 
		tz = ""
	)
	v <- v[order(v$date),]
	v$time <- 0:1235
	rhcollect(k, v)
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
	trellis.device(
		device = postscript, 
		file = paste(
			"./tmp/"
			par$dataset 
			"_measurement_status_against_month_part_",
			reduce.key,
			".ps", 
			sep = ""
		), 
		color = TRUE, 
		paper = "legal"
	)
      b <- dotplot( as.numeric(obs) ~ time | station.id,
		data = combined,
        xlab = list(label = "Month", cex = 1.5),
        ylab = list(label = ylab, cex = 1.5),
        type = "p",
        col = c("black","red"),
        distribute.type= TRUE,
        groups = obs,
        pch=16,
        cex=0.3,
        layout = c(1,10),
        strip.left = TRUE,
        strip = FALSE,
        grib = TRUE,
        scales = list(
        	y = list(
             	relation = 'same', 
             	cex = 1, 
             	labels = c("0","1"), 
             	alternating = TRUE
            ), 
            x = list(
             	tick.number = 10, 
             	relation = 'same', 
             	cex = 1.2
            )
		),
        panel = function(x,y,...) {
            panel.abline( 
                v = seq(0,1235, by=60), 
                color = "lightgrey", 
                lty = 3, 
                lwd = 0.5
            )
            panel.dotplot(x,y,...)
        }
      )
      print(b)
	dev.off()
  }
)
job$setup <- expression(
    map = {
		library(lattice)
    },
)
job$parameters <- list(
	par = par,
	ylab = ylab
)
job$input <- rhfmt(
	paste(rh.datadir, par$dataset, sep = ""), 
	type = "sequence"
),
job$output <- rhfmt(
	paste(rh.output, par$dataset, sep = ""), 
	type = "sequence"
),
job$mapred <- list(
	mapred.map.tasks = 50, 
	mapred.reduce.tasks = 50
),
job$jobname <- paste(rh.output, par$dataset, sep=""),
job$readback <- FALSE
job$copyFiles <- TRUE
job.mr <- do.call("rhwatch", job)