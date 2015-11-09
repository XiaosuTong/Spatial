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
  file.path(rh.root, par$dataset, "All", "bymonth"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.root, par$dataset, "All", "monthcount"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 1,  #cdh3,4
  mapreduce.job.reduces = 1  #cdh5
)
job$readback <- FALSE
job$combiner <- TRUE
job$jobname <- file.path(rh.root, par$dataset, "All", "monthcount")
job.mr <- do.call("rhwatch", job)

rst <- rhread(file.path(rh.root, par$dataset, "All", "monthcount"))[[1]][[2]]

library(plyr)
library(lattice)
rst$month <- factor(rst$month, levels = month.abb)
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

trellis.device(
  device = postscript,
  file = file.path(local.root, "output", "obs.month.byseason.ps"),
  color = TRUE,
  paper = "letter"
)
  b <- xyplot(log2(count) ~ as.numeric(year)
    , data = rst
    , type = "l"
    , lwd = 2.5
    , xlab = list(label = "Year", cex = 1.5)
    , ylab = list(label = "Log Base 2 Number of Observation", cex = 1.5)
    , ylim = c(9, 13.5)
    , scales = list(
        y = list(cex = 1.2),
        x = list(at = c(seq(1895, 2000, by=20)), cex = 1.2)
      )  
    , key=list(
        text = list(label=c("Jan","Feb", "Dec")),
        lines = list(lwd=2.5, type="l", col=col[1:3], lty=c(1,2,4)), 
        columns = 3
      )
    , panel = function(x,y,...) {
        panel.abline(v = seq(1895,2000, by=10), h = seq(9.5, 13.5, 0.5), lwd = 0.5, col = "lightgray")
        panel.abline(h = log2(Nstations), col = "red", lty=1, lwd = 1)
        for(i in c("Jan","Dec","Feb")) {
          sub <- subset(rst, month == i)
          if (i == "Jan") {
            panel.xyplot(x=as.numeric(sub$year), y=log2(sub$count), lty = 1, col=col[1],...)
          } else if(i == "Feb") {
            panel.xyplot(x=as.numeric(sub$year), y=log2(sub$count), lty = 2, col=col[2],...)
          } else {
            panel.xyplot(x=as.numeric(sub$year), y=log2(sub$count), lty = 4, col=col[3],...)
          }
        }
      }
  )
  print(b)
dev.off()
