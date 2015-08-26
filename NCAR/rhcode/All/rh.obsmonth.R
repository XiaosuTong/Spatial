source("~/Rhipe/rhinitial.R")

Machine <- "rossmann"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

par <- list()
par$dataset <- "tmax"

if (Machine == "adhara") {
  root <- "/ln/tongx/Spatial/tmp"
} else if (Machine == "rossmann") {
  root <- "/wsc/tongx/Spatial/tmp"
}
if(par$dataset == "precip") {
  Nstations <- 11918
} else {
  Nstations <- 8125
}

job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
    count <- sum(!is.na(map.values[[r]]$resp))
    value <- data.frame(
      count = count,
      year = map.keys[[r]][1],
      month = map.keys[[r]][2],
      stringsAsFactors = FALSE
    )
    rhcollect(1, value)
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
  file.path(root, par$dataset, "All", "bymonth"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(root, par$dataset, "All", "monthcount"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 1,  #cdh3,4
  mapreduce.job.reduces = 1  #cdh5
)
job$readback <- FALSE
job$combiner <- TRUE
job$jobname <- file.path(root, par$dataset, "All", "monthcount")
job.mr <- do.call("rhwatch", job)

rst <- rhread(file.path(root, par$dataset, "All", "monthcount"))[[1]][[2]]

library(plyr)
library(lattice)
month <- c(
  "Jan","Feb","Mar","Apr","May","June",
  "July","Aug", "Sep", "Oct", "Nov", "Dec"
)
rst$month <- factor(rst$month, levels = month)
date <- paste(rst$year, as.numeric(rst$month), "01", sep="-")
rst$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
rst <- arrange(rst, year, month)
start <- head(rst$date, 1)
end <- tail(rst$date, 1)


trellis.device(
  device = postscript, 
  file = file.path(local.output, "obs.month.ps"), 
  color = TRUE, 
  paper = "letter"
)
  b <- xyplot(log2(count) ~ date
    , data = rst
    , type = "l"
    , xlab = list(label = "Date", cex = 1.5)
    , ylab = list(label = "Log Base 2 Number of Observation", cex = 1.5)
    , ylim = c(9, 13.5)
    , aspect = "xy"
    , scales = list(
        y = list(cex = 1.2),
        x = list(format = "%b %Y", at = c(seq(start, end, by="240 month")), cex = 1.2)
      )   
    , panel = function(x,y,...) {
    	  panel.abline(v = seq(start, end, by="120 month"), h = seq(9.5, 13.5, 0.5), lwd = 0.5, col = "lightgray")
        panel.abline(h = log2(Nstations), col = "red", lty=1, lwd = 1)
        panel.xyplot(x,y,...)
      }
  )
  print(b)
dev.off()