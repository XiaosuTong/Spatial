source("~/Rhipe/ross.initial.R")

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
  file.path(root, par$dataset, "All", "bystation"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(root, par$dataset, "All", "fullobs.station"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 1,  #cdh3,4
  mapreduce.job.reduces = 1  #cdh5
)
job$readback <- FALSE
job$combiner <- TRUE
job$jobname <- file.path(root, par$dataset, "All", "fullobs.station")
job.mr <- do.call("rhwatch", job)

rst <- rhread(file.path(root, par$dataset, "All", "fullobs.station"))[[1]][[2]]

set.seed(100)
stations100 <- sample(rst$station.id, 100)

rhsave(stations100, file = file.path(root, par$dataset, "100stations", "Rdata", "100stations.RData"))




#######################################################
##
##
#######################################################

library(maps)
us.map <- map('state', plot = FALSE, fill = TRUE)

rhload(file.path(root, par$dataset, "100stations", "Rdata", "100stations.RData"))
rst <- rhread(file.path(root, par$dataset, "All", "fullobs.station"))[[1]][[2]]

trellis.device(
  device = postscript, 
  file = file.path(local.output, "fullobs.station.ps"), 
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
  file = file.path(local.output, "100stations.ps"), 
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


#################################################################
##
#################################################################

library(plyr)

rhload(file.path(root, "stationinfo", "USinfo.RData"))
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
  file = file.path(local.output, "allstations.ps"), 
  color = TRUE, 
  paper = "letter"
)
  b <- xyplot( lat ~ lon | factor(group, label = label$label)
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


####################################################################
##
####################################################################
library(lattice)
library(plyr)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

rhload(file.path(root, "stationinfo", "USinfo.RData"))
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
  file = file.path(local.output, "elevrange.ps"),, 
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
  file = file.path(local.output, "QQelevation.ps"), 
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
  file = file.path(local.output, "QQelevation.log.ps"), 
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


#################################################################################
##
#################################################################################
data1 <- subset(USppt, year=="1981"&month=="Dec")
data1 <- data1[!is.na(data1$precip),c("station.id","lat","lon")]
data2 <- subset(USppt, year=="1982"&month=="Jan")
data2 <- data2[!is.na(data2$precip),c("station.id", "lat", "lon")]
data3 <- subset(USppt, year=="1982"&month=="Feb")
data3 <- data3[!is.na(data3$precip),c("station.id", "lat", "lon")]
miss <- data1[!(data1$station.id %in% data2$station.id),]
miss$group <- "miss" #miss is the stations that active in Dec 1981 but not in Jan 1982
add <- data2[!(data2$station.id %in% data1$station.id),]
add$group <- "add" #add is the stations that active in Jan 1982 but not in Dec 1981
addback <- data3[!(data3$station.id %in% data2$station.id),]
miss2 <- data2[!(data2$station.id %in% data3$station.id),]
miss2$group <- "miss2" #miss2 is the stations that active in Jan 1982 but not in Feb 1982
addback$group <- "addback"#addback is the stations that active in Feb 1982 but not in Jan1982
data <- rbind(addback, miss, add, miss2)
data$group <- factor(data$group)
trellis.device(
    postscript, 
    file = paste(
        outputdir, 
        "spatial.precip.around1982", ".ps", sep = ""), 
    color = TRUE, 
    paper = "legal"
)
a <- xyplot(
    lat ~ lon,
#    lat ~ lon | group,
    data = data,
    xlab = list(label="Longitude"),
    ylab = list(label="Latitude"),
    main = "Spatial Location of Precipitation Stations",
    layout = c(1,1),
    key = list(
        columns = 4, 
        points = list(pch=c(1,16,1,16), col = col[c(1,3,4,2)]), 
        text = list(label=c("Jan 1982 Missed","Jan 1982 Added", "Feb 1982 Missed", 
            "Feb 1982 Added")#miss add miss2 addback
        )
    ),
    groups = group,
    col = col[c(3,2,1,4)], ##add addback miss miss2
    pch = c(16,16,1,1),
    cex = 0.7,
    panel = function(...) {
        panel.polygon(us.map$x,us.map$y)   
        panel.xyplot(...)
    }
)
print(a)
dev.off()

