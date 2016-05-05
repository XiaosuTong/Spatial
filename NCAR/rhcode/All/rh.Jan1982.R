###################################################################
##  Check the drop of count on Jan 1982, and Dec 1981, Feb 1982  ##
###################################################################
job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
    if(all(map.keys[[r]] == c("1982","Jan")) || all(map.keys[[r]] == c("1982","Feb"))|| all(map.keys[[r]] == c("1981","Dec")) ) {
      value <- map.values[[r]]
      value <- subset(value, !is.na(resp))
      rhcollect(map.keys[[r]], value)
    }  
  })
})
job$reduce <- expression(
  pre = {
    combine <- data.frame()
  },
  reduce = {
    combine <- rbind(combine, do.call("rbind", reduce.values))
  },
  post = {
    rhcollect(reduce.key, combine)
  }
)
job$input <- rhfmt(
  file.path(rh.root, par$dataset, "All", "bymonth"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.root, par$dataset, "All", "Jan1982"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 1,  #cdh3,4
  mapreduce.job.reduces = 1  #cdh5
)
job$readback <- FALSE
job$mon.sec <- 10
job$combiner <- TRUE
job$jobname <- file.path(rh.root, par$dataset, "All", "Jan1982")
job.mr <- do.call("rhwatch", job)


rst <- rhread("/ln/tongx/spatem/tmax/All/Jan1982")

IDec1981 <- which(unlist(lapply(lapply(rst, "[[", 1), function(r){all(r==c("1981","Dec"))})))
IJan1982 <- which(unlist(lapply(lapply(rst, "[[", 1), function(r){all(r==c("1982","Jan"))})))
IFeb1982 <- which(unlist(lapply(lapply(rst, "[[", 1), function(r){all(r==c("1982","Feb"))})))

miss <- subset(rst[[IDec1981]][[2]], !(station.id %in% rst[[IJan1982]][[2]]$station.id))
miss$group <- "miss" #miss is the stations that active in Dec 1981 but not in Jan 1982

add <- subset(rst[[IJan1982]][[2]], !(station.id %in% rst[[IDec1981]][[2]]$station.id))
add$group <- "add" #add is the stations that active in Jan 1982 but not in Dec 1981

addback <- subset(rst[[IFeb1982]][[2]], !(station.id %in% rst[[IJan1982]][[2]]$station.id))
addback$group <- "addback"#addback is the stations that active in Feb 1982 but not in Jan1982

miss2 <- subset(rst[[IJan1982]][[2]], !(station.id %in% rst[[IFeb1982]][[2]]$station.id))
miss2$group <- "miss2" #miss2 is the stations that active in Jan 1982 but not in Feb 1982

data <- rbind(miss, addback)
#data <- rbind(addback, miss, add, miss2)
data$group <- factor(data$group)

us.map <- map('state', plot = FALSE, fill = TRUE)

trellis.device(
  device = postscript, 
  file = file.path(local.root, "output", "obs.Jan1982.ps"), 
  color = TRUE, 
  paper = "letter"
)
  a <- xyplot(lat ~ lon
    , data = data
    , groups = group
    , xlab = list(label="Longitude", cex = 1.5)
    , ylab = list(label="Latitude", cex = 1.5)
    , scale = list(cex=1.2)
    , xlim = c(-125.5, -66.5)
    , aspect = 0.66
    , key = list(
        columns = 2, 
        cex = 1.2,
        points = list(pch=c(16,1), col = col[c(2,1)]), 
        text = list(
          label=c("Feb 1982 Active", "Jan 1982 Missed")#miss add miss2 addback
        )
      )
    , col = col[c(2, 1)], ##add addback miss miss2
    , pch = c(16,1),
    , cex = 1.2,
    , panel = function(x,y,...) {
        panel.abline(h=seq(25,50,5),v=seq(-120,-70,10), lty=1, lwd=0.5, col="lightgray")
        panel.polygon(us.map$x,us.map$y)   
        panel.xyplot(x,y,...)
      }
  )
  print(a)
dev.off()

