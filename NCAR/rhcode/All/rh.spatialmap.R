  ###############################################################
##  Create division by station which includes all full obs   ##
##  Then create the 100stations.RData which contain the 100  ##
##  stations.id  saved on HDFS                               ##
###############################################################
job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
    if(sum(is.na(map.values[[r]]$resp)) == 0) {
      value <- data.frame(
        station.id = map.keys[[r]],
        elev = as.numeric(attributes(map.values[[r]])$location["elev"]),
        lon = as.numeric(attributes(map.values[[r]])$location["lon"]),
        lat = as.numeric(attributes(map.values[[r]])$location["lat"]),
        stringsAsFactors = FALSE
      )
      rhcollect(1, value)
    }
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
    row.names(combined) <- NULL
    rhcollect(reduce.key, combined)
  }
)
job$input <- rhfmt(
  file.path(rh.root, par$dataset, "All", "bystation"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.root, par$dataset, "All", "fullobs.station"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 1,  #cdh3,4
  mapreduce.job.reduces = 1  #cdh5
)
job$readback <- FALSE
job$combiner <- TRUE
job$jobname <- file.path(rh.root, par$dataset, "All", "fullobs.station")
job.mr <- do.call("rhwatch", job)

rst <- rhread(file.path(rh.root, par$dataset, "All", "fullobs.station"))[[1]][[2]]

set.seed(100)
stations100 <- sample(rst$station.id, 100)

rhsave(stations100, file = file.path(rh.root, par$dataset, "100stations", "Rdata", "100stations.RData"))

####################################################################
##  Create divition by station which includes all a1950 stations  ##
##  Then save create the a1950.RData which contain the a1950      ##
##  stations.id  saved on HDFS                                    ##
####################################################################
job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
    sub <- subset(map.values[[r]], year >= 1950)
    if(sum(is.na(sub$resp)) < 576) {
      rhcollect(1, map.keys[[r]])
    }
  })
})
job$reduce <- expression(
  pre = {
    combine <- character()
  },
  reduce = {
    combine <- c(combine, unlist(reduce.values))
  },
  post = {
    rhcollect(reduce.key, combine)
  }
)
job$input <- rhfmt(
  file.path(rh.root, par$dataset, "All", "bystation"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.root, par$dataset, "All", "a1950.station"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 1,  #cdh3,4
  mapreduce.job.reduces = 1  #cdh5
)
job$readback <- FALSE
job$combiner <- TRUE
job$jobname <- file.path(rh.root, par$dataset, "All", "a1950.station")
job.mr <- do.call("rhwatch", job)
## It is a character vector with length 7738
stations.a1950 <- rhread(file.path(rh.root, par$dataset, "All", "a1950.station"))[[1]][[2]]

rhsave(stations.a1950, file = file.path(rh.root, par$dataset, "a1950", "Rdata", "a1950stations.RData"))

############################################################
##  spatial map of station location for fullobs stations  ##
##  and 100 stations                                      ##
############################################################

library(maps)
us.map <- map('state', plot = FALSE, fill = TRUE)

rhload(file.path(rh.root, par$dataset, "100stations", "Rdata", "100stations.RData"))
rst <- rhread(file.path(rh.root, par$dataset, "All", "fullobs.station"))[[1]][[2]]

trellis.device(
  device = postscript, 
  file = file.path(local.root, "output", "fullobs.station.ps"), 
  color = TRUE, 
  paper = "letter"
)
  b <- xyplot( lat ~ lon
    , data = rst
    , xlab = list(label="Longitude", cex = 1.5)
    , ylab = list(label="Latitude", cex = 1.5)
    , layout = c(1,1)
    , pch = 16
    , cex = 0.5
    , col = "red"
    , scales = list(
        y = list(cex = 1.2),
        x = list(cex = 1.2)
      )
#   , strip = strip.custom(par.strip.text= list(cex = 1.5)),
#   , par.settings = list(layout.heights = list(strip = 1.5)),
    , panel = function(x,y,...) {
        panel.polygon(us.map$x,us.map$y)   
        panel.xyplot(x,y,...)
      }
  )
  print(b)
dev.off()

trellis.device(
  device = postscript, 
  file = file.path(local.root, "output", "100stations.ps"), 
  color = TRUE, 
  paper = "letter"
)
  b <- xyplot( lat ~ lon
    , data = rst
    , subset = station.id %in% stations100
    , xlab = list(label="Longitude", cex = 1.5)
    , ylab = list(label="Latitude", cex = 1.5)
    , layout = c(1,1)
    , pch = 16
    , cex = 0.5
    , col = "red"
    , scales = list(
        y = list(cex = 1.2),
        x = list(cex = 1.2)
      )
#   , strip = strip.custom(par.strip.text= list(cex = 1.5)),
#   , par.settings = list(layout.heights = list(strip = 1.5)),
    , panel = function(x,y,...) {
        panel.polygon(us.map$x,us.map$y)   
        panel.xyplot(x,y,...)
      }
  )
  print(b)
dev.off()

#############################################################
##  Spatial map of station location for station after1950  ##
#############################################################
rhload(file.path(rh.root, "stationinfo", "USinfo.RData"))
rhload(file.path(rh.root, par$dataset, "a1950", "Rdata", "a1950stations.RData"))

if(par$dataset == "precip") {
  info <- USpinfo
} else {
  info <- UStinfo
}
trellis.device(
  device = postscript, 
  file = file.path(local.root, "output", "a1950stations.ps"), 
  color = TRUE, 
  paper = "letter"
)
  b <- xyplot( lat ~ lon
    , data = info
	, subset = station.id %in% a1950.stations
    , xlab = list(label="Longitude", cex = 1.5)
    , ylab = list(label="Latitude", cex = 1.5)
    , pch = 16
    , cex = 0.4
    , col = "red"
    , scales = list(
        y = list(cex = 1.2),x = list(cex = 1.2)
      )
    , panel = function(x,y,...) {
        panel.polygon(us.map$x,us.map$y)   
        panel.xyplot(x,y,...)
      }
  )
  print(b)
dev.off()

##########################################################
##  spatial map of station location for all stations.   ##
##  conditional on elevation                            ##   
##########################################################
library(plyr)
rhload(file.path(rh.root, "stationinfo", "USinfo.RData"))
if(par$dataset == "precip") {
  info <- USpinfo
  group <- c(rep(1:7, each = 1490), rep(8,1488))
} else {
  info <- UStinfo
  group <- c(rep(1:7, each = 1015), rep(8, 1020))
}

info <- arrange(info, elev)
info$group <- group

label <- ddply(
  .data = info,
  .vari = "group",
  .fun  = summarise,
  label = paste(min(elev), "m", " ~ ", max(elev), "m", sep="")
)

trellis.device(
  device = postscript, 
  file = file.path(local.root, "output", "allstationsone.ps"), 
  color = TRUE, 
  paper = "letter"
)
  b <- xyplot( lat ~ lon #| factor(group, label = label$label)
    , data = info
    , xlab = list(label="Longitude", cex = 1.5)
    , ylab = list(label="Latitude", cex = 1.5)
    , layout = c(1,1)
    , pch = 16
    , cex = 0.4
    , col = "red"
    , scales = list(
        y = list(cex = 1.2),
        x = list(cex = 1.2)
      )
#   , strip = strip.custom(par.strip.text= list(cex = 1.5)),
#   , par.settings = list(layout.heights = list(strip = 1.5)),
    , panel = function(x,y,...) {
        panel.polygon(us.map$x,us.map$y)   
        panel.xyplot(x,y,...)
      }
  )
  print(b)
dev.off()


################################################
##  diagnostic for elevation of all stations  ## 
################################################
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

rhload(file.path(rh.root, "stationinfo", "USinfo.RData"))
if(par$dataset == "precip") {
  info <- USpinfo
  group <- c(rep(1:7, each = 1490), rep(8,1488))
} else {
  info <- UStinfo
  group <- c(rep(1:7, each = 1015), rep(8, 1020))
}
info <- arrange(info, elev)
info$group <- group

rg <- ddply(
  .data = info,
  .vari = "group",
  .fun  = summarise,
  type  = c("min","max"),
  elev  = range(elev),
  label = rep(paste(min(elev),max(elev), sep="~"),2)
)

trellis.device(
  device = postscript, 
  file = file.path(local.root, "output", "elevrange.ps"),, 
  color=TRUE, 
  paper="letter"
)
  b <- dotplot( factor(label, levels = unique(label)) ~ elev
    , data = rg
    , xlab = list(label="Elevation (meter)", cex=1.5)
    , ylab = list(label="Level", cex=1.5)
    , groups = type
    , pch = 16
    , grib = TRUE
    , key = list(
        columns=2, 
        points= list(pch=16, col = col[1:2]), 
        text= list(label=c("Max","Min"), cex=1.2)
      )
    , scales = list(
        x = list(cex=1.2), 
        y = list(cex=1.2)
      )
    , panel = function(x,y,...) {
        panel.dotplot(x,y,...)
      }
  )
  print(b)
dev.off()


trellis.device(
  device = postscript, 
  file = file.path(local.root, "output", "QQelevation.ps"), 
  color=TRUE, 
  paper="letter"
)
  a <- qqmath( ~ elev,
    , data = info,
    , distribution = qunif,
   # , aspect = 1,
    , pch = 16,
    , cex = 0.3,
    , xlab = list(label="f-value", cex=1.5),
    , ylab = list(label="Elevation (meter)", cex=1.5),
    , scales = list(
        x = list(cex=1.2), 
        y = list(cex=1.2)
      )
    , prepanel = prepanel.qqmathline,
    , panel = function(x, y,...) {
        panel.abline( h = rg$elev, color = "lightgray", lty = 2, lwd = 0.5)
        panel.qqmath(x, y,...)
      }
  )
  print(a)
dev.off()

trellis.device(
  device = postscript, 
  file = file.path(local.root, "output", "QQelevation.log.ps"), 
  color=TRUE, 
  paper="letter"
)
  a <- qqmath( ~ log2(elev),
    , data = subset(info, elev >= 1)
    , distri = qunif
    , pch = 16,
    , cex = 0.3,
    , xlab = list(label="f-value", cex=1.5),
    , ylab = list(label="Log Base 2 Elevation", cex=1.5),
    , scales = list(
        x = list(cex=1.2), 
        y = list(cex=1.2)
      )
    , prepanel = prepanel.qqmathline,
    , panel = function(x, y,...) {
        panel.abline(h = log2(rg$elev[-1]), color="lightgrey", lty=2, lwd=0.5)
        panel.qqmath(x, y,...)
      }
  )
  print(a)
dev.off()
