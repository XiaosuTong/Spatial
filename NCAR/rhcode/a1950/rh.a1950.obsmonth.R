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
  file.path(rh.root, par$dataset, "a1950", "bymonth"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.root, par$dataset, "a1950", "monthcount"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 1,  #cdh3,4
  mapreduce.job.reduces = 1  #cdh5
)
job$readback <- FALSE
job$combiner <- TRUE
job$jobname <- file.path(rh.root, par$dataset, "a1950", "monthcount")
job.mr <- do.call("rhwatch", job)

rst <- rhread(file.path(rh.root, par$dataset, "a1950", "monthcount"))[[1]][[2]]

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
  file = file.path(local.root, "output", "a1950.obs.month.ps"), 
  color = TRUE, 
  paper = "letter"
)
  b <- xyplot(log2(count) ~ date
    , data = rst
    , type = "l"
    , xlab = list(label = "Date", cex = 1.5)
    , ylab = list(label = "Log Base 2 Number of Observation", cex = 1.5)
    , scales = list(
        y = list(cex = 1.2),
        x = list(format = "%b %Y", at = c(seq(start, end, by="120 month")), cex = 1.2)
      )   
    , panel = function(x,y,...) {
        panel.abline(v = seq(start, end, by="120 month"), h = seq(9.5, 13.5, 0.1), lwd = 0.5, col = "lightgray")
        panel.xyplot(x,y,...)
      }
  )
  print(b)
dev.off()

rst <- rhread("/ln/tongx/spatem/tmax/a1950/bymonth")
trellis.device(
  device = postscript, 
  file = file.path(local.root, "output", paste(par$dataset, "a1950", "status", "ps", sep=".")), 
  color = TRUE, 
  paper = "letter"
)
for(i in 1:576) {
  b <- xyplot( lat ~ lon
    , data = rst[[i]][[2]]
    , groups = factor(as.numeric(!is.na(resp)))
    , xlab = list(label="Longitude", cex = 1.5)
    , ylab = list(label="Latitude", cex = 1.5)
    , sub = paste(rst[[i]][[1]][2], rst[[i]][[1]][1])
    , key = list(
        type = "p", 
        text = list(label=c("missing","valid")),  
        points = list(cex=1, pch=16, col=col[1:2]), 
        columns = 2
      )
    , pch = 16
    , xlim = c(-125.5, -66.5)
    , aspect = 0.66
    , cex = 0.4
    , scales = list(
        y = list(cex = 1.2),x = list(cex = 1.2)
      )
    , panel = function(x,y,...) {
        panel.abline(h=seq(25,50,5),v=seq(-120,-70,10), lty=1, lwd=0.5, col="lightgray")
        panel.polygon(us.map$x,us.map$y)   
        panel.xyplot(x,y,...)
        panel.text(x=-123, y=26, labels=paste("Missing:", sum(is.na(rst[[i]][[2]]$resp))), adj=c(0,0))
        panel.text(x=-123, y=25.5, labels=paste("Valid:    ", sum(!is.na(rst[[i]][[2]]$resp))), adj=c(0,1))
      }
  )
  print(b)
}
dev.off()