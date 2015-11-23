rst <- rhread("/ln/tongx/Spatial/tmp/tmax/All/Jan1982")
IDec1981 <- which(unlist(lapply(lapply(rst, "[[", 1), function(r){all(r==c("1981","Dec"))})))
IJan1982 <- which(unlist(lapply(lapply(rst, "[[", 1), function(r){all(r==c("1982","Jan"))})))
IFeb1982 <- which(unlist(lapply(lapply(rst, "[[", 1), function(r){all(r==c("1982","Feb"))})))
miss <- subset(rst[[IDec1981]][[2]], !(station.id %in% rst[[IJan1982]][[2]]$station.id))
me <- as.character(miss$station.id)


job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
    value <- subset(map.values[[r]], station.id %in% target)
    m <- mean(value$resp, na.rm = TRUE)
    value <- data.frame(
      mean = m,
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
  file.path(rh.root, par$dataset, "All", "bymonth"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.root, par$dataset, "All", "monthmean"), 
  type = "sequence"
)
job$parameters <- list(
  target = me
)
job$mapred <- list(
  mapred.reduce.tasks = 1,  #cdh3,4
  mapreduce.job.reduces = 1  #cdh5
)
job$readback <- FALSE
job$combiner <- TRUE
job$jobname <- file.path(rh.root, par$dataset, "All", "monthmean")
job.mr <- do.call("rhwatch", job)


rst <- rhread(file.path(rh.root, par$dataset, "All", "monthmean"))[[1]][[2]]

library(plyr)
library(lattice)
rst$month <- factor(rst$month, levels = month.abb)
date <- paste(rst$year, as.numeric(rst$month), "01", sep="-")
rst$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
rst <- arrange(rst, year, month)
start <- head(rst$date, 1)
end <- tail(rst$date, 1)

rst <- subset(rst, year >= 1931)
rst$resp <- ts(rst$mean, start=c(1931,1), frequency=12)

trellis.device(
  device = postscript, 
  file = file.path(local.root, "output", "monthmean.ps"), 
  color = TRUE, 
  paper = "letter"
)
  b <- xyplot(rst$resp
    , type = "l"
    , strip = FALSE
    , xlab = list(label = "Date", cex = 1.5)
    , ylab = list(label = ylab, cex = 1.5)
    , cut = list(number = 5, overlap = 0)
    , scales = list(
        y = list(cex = 1.2),
        x = list(at=seq(1931, 1998, 2), cex = 1.2)
      )
    , ylim = c(-20,20)
    , panel = function(x,y,...) {
        panel.abline(v = seq(1931, 1998, by=4), lwd = 0.5, col = "lightgray")
        panel.abline(h=min(rst$resp), lwd=0.5, col="red")
        panel.xyplot(x,y,...)
        panel.points(x=1982, y=min(rst$resp), col="red")
        panel.text(x=1982, y=min(rst$resp), labels=round(min(rst$resp),2), adj=c(0,0))
      }
  )
  print(b)
dev.off()
