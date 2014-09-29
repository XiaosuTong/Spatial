##################################################
##Copy raw text file to HDFS, 
##and create division by year, division by station
##################################################
source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "adhara"
par$dataset <- "tmax"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

for(x in formatC(1:103, width = 3, flag = "0")) {
  rhput(
    paste(local.raw,"/NCAR_pinfill/ppt.complete.Y", x, sep = ""), 
    paste(rh.datadir,"/Raw/precip/ppt.complete.Y", x, sep = "")
  )
}
for(x in formatC(1:103, width = 3, flag = "0")) {
  rhput(
    paste(local.raw, "/NCAR_tinfill/tmax.complete.Y", x, sep = ""), 
    paste(rh.datadir,"/Raw/tmax/tmax.complete.Y", x, sep = "")
  )
  rhput(
    paste(local.raw,"/NCAR_tinfill/tmin.complete.Y", x, sep = ""), 
    paste(rh.datadir,"/Raw/tmin/tmin.complete.Y", x, sep="")
  )
}


job <- list()
job$map <- expression({
  y <- do.call("rbind", 
    lapply(map.values, function(r) {
      row <- strsplit(r, " +")[[1]]
      c(row[1], row[2:13], substring(row[14], 1:12, 1:12))
    })
  )
  file <- Sys.getenv("mapred.input.file") #get the file name that Hadoop is reading
  k <- as.numeric(
    substr(tail(strsplit(tail(strsplit(file, "/")[[1]],1), "[.]")[[1]], 1), 2, 4)
  )
  miss <- as.data.frame(
    matrix(as.numeric(y[, (1:12) + 13]), ncol = 12)
  )
  tmp <- as.data.frame(
    matrix(as.numeric(y[, (1:12) + 1]), ncol = 12)
  )
  name <- y[, 1]
  tmp <- tmp/10
  tmp[miss == 1] <- NA
  names(tmp) <- c(
    "Jan", "Feb", "Mar", "Apr", "May", "June", 
    "July", "Aug", "Sep", "Oct", "Nov", "Dec"
  )
  tmp <- cbind(
    station.id = name, 
    tmp, 
    year = rep((k + 1894))
  )
  if(par$dataset != "precip") {
    value <- data.frame(
      station.id = rep(tmp$station.id, 12),
      elev = rep(UStinfo$elev, 12),
      lon = rep(UStinfo$lon, 12),
      lat = rep(UStinfo$lat, 12),
      month = rep(names(tmp)[2:13], each = dim(tmp)[1]),
      resp = c(
        tmp[, 2], tmp[, 3], tmp[, 4], tmp[, 5], tmp[, 6], tmp[, 7], 
        tmp[, 8], tmp[, 9], tmp[, 10], tmp[, 11], tmp[, 12], tmp[, 13]
      )
    )
  }else {
    value <- data.frame(
      station.id = rep(tmp$station.id, 12),
      elev = rep(USpinfo$elev, 12),
      lon = rep(USpinfo$lon, 12),
      lat = rep(USpinfo$lat, 12),
      year = rep(tmp$year,12),
      month = rep(names(tmp)[2:13], each = dim(tmp)[1]),
      resp = c(
        tmp[, 2], tmp[, 3], tmp[, 4], tmp[, 5], tmp[, 6], tmp[, 7], 
        tmp[, 8], tmp[, 9], tmp[, 10], tmp[, 11], tmp[, 12], tmp[, 13]
      )
    )    
  }
  rhcollect(unique(tmp$year), value)
})
job$shared <- c("/ln/tongx/Spatial/tmp/stationinfo/USinfo.RData")
job$setup <- expression(
  map = {
    load("USinfo.RData")
  }
)
job$parameters <- list(
  par = par
)
job$input <- rhfmt(
  file.path(
    rh.datadir,
    "Raw", par$dataset 
  ),
  type = "text"
)
job$output <- rhfmt(
  file.path(
    rh.datadir,
    par$dataset,
    "data",
    "byyear"
  ), 
  type = "sequence"
)
job$mapred <- list( 
  mapred.reduce.tasks = 72, 
  rhipe_map_buff_size = 8125 
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)

job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
    lapply(1:nrow(map.values[[r]]), function(x) {
      key <- as.character(map.values[[r]][x, 1])
      value <- map.values[[r]][x, -1]
      value$year <- map.keys[[r]]
      rhcollect(key, value)
    })
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
  file.path(
    rh.datadir,
    par$dataset,
    "data",
    "byyear"
  ), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(
    rh.datadir,
    par$dataset,
    "data",
    "bystation"
  ), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 72
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)

