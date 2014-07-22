#############################################################
##dot plot for each station to see the where are missing
#############################################################

#Set up the directory and load the data.
#library(lattice)
#library(plyr)
#local.datadir <- "~/Projects/Spatial/NCAR/RData/"
#local.output <- "/home/shaula/u16/tongx/Projects/Spatial/NCAR/output/"
#rh.datadir <- "/ln/tongx/Spatial/tmp/"
#The output directory for saving plots should have all permission since the plots are written to the output directory by a user related to hadoop.
#In the reduce step, the permission for the plots should be changed to be all +wrx.
#rh.output <- "/ln/tongx/Spatial/output/"

map <- expression({
  lapply(seq_along(map.values), function(r) {
      v <- map.values[[r]]
      k <- map.values[[r]]$partion[1]
      v$obs <- !is.na(v$tmax)
      month <- v$month
      levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
      month <- as.numeric(factor(month, levels=c(1:12)))
      date <- paste(v$year, month, "01", sep="-")
      v$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
      v <- v[order(v$date),]
      v$time <- 0:1235
      rhcollect(k, v)
  })
})


reduce <- expression(
   pre={
	combined <- data.frame()
   },
   reduce={
	combined <- rbind(combined, do.call(rbind, reduce.values))
   },
   post={
	library(lattice)
      	trellis.device(postscript, file = paste(local.output, "tmax_measurement_status_against_month_part_",reduce.key,".ps", sep = ""), color=TRUE, paper="legal")
        b <- dotplot( as.numeric(obs) ~ time | station.id,
             data = combined,
             xlab = list(label = "Month", cex = 1.5),
             ylab = list(label = "Maximum Temperature Measurement Status", cex = 1.5),
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
#             scales = list(y = list(relation = 'same', cex=1.2, alternating=TRUE), x=list(relation= 'same',format = "%b %Y", tick.number=10), cex=1.2),
             scales = list(y = list(relation = 'same', cex=1, labels=c("0","1"), alternating=TRUE), x=list(tick.number = 10, relation='same', cex=1.2)),
             panel = function(x,y,...) {
                  panel.abline( v=seq(0,1235, by=60), color="lightgrey", lty=3, lwd=0.5)
                  panel.dotplot(x,y,...)
             }
        )
        print(b)
	dev.off()
	system(paste("chmod 777 ", local.output, "tmax_measurement_status_against_month_part_",reduce.key,".ps", sep = ""))
   }
)

ex <- rhwatch(
        map = map,
        reduce = reduce,
        input = rhfmt(paste(rh.datadir,"tmax",sep=""), type="sequence"),
        output = rhfmt(paste(rh.output,"tmax",sep=""), type="sequence"),
        mapred = list(mapred.map.tasks = 50, mapred.reduce.tasks = 50),
        jobname = "tmax",
        readback = FALSE)

