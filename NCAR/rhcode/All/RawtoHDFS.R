##################################################
##Copy raw text file to HDFS, 
##and create division by month, division by station
##################################################
for(x in formatC(1:103, width = 3, flag = "0")) {
  rhput(
    paste(local.raw,"/NCAR_pinfill/ppt.complete.Y", x, sep = ""), 
    paste(rh.root,"/Raw/precip/ppt.complete.Y", x, sep = "")
  )
}
for(x in formatC(1:103, width = 3, flag = "0")) {
  rhput(
    paste(local.raw, "/NCAR_tinfill/tmax.complete.Y", x, sep = ""), 
    paste(rh.root,"/Raw/tmax/tmax.complete.Y", x, sep = "")
  )
  rhput(
    paste(local.raw,"/NCAR_tinfill/tmin.complete.Y", x, sep = ""), 
    paste(rh.root,"/Raw/tmin/tmin.complete.Y", x, sep="")
  )
}


job <- list()
job$map <- expression({
  if(par$dataset == "precip") {
    info <- USpinfo
  } else {
    info <- UStinfo
  }
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
  names(tmp) <- month.abb
  tmp <- cbind(station.id = name, tmp, year = rep((k + 1894)))
  value <- data.frame(
    station.id = rep(tmp$station.id, 12),
    elev = rep(info$elev, 12),
    lon = rep(info$lon, 12),
    lat = rep(info$lat, 12),
    month = rep(names(tmp)[2:13], each = dim(tmp)[1]),
    resp = c(
      tmp[, 2], tmp[, 3], tmp[, 4], tmp[, 5], tmp[, 6], tmp[, 7], 
      tmp[, 8], tmp[, 9], tmp[, 10], tmp[, 11], tmp[, 12], tmp[, 13]
    ),
    stringsAsFactors = FALSE
  )
  d_ply(
    .data = value,
    .vari = "month",
    .fun = function(r) {
      rhcollect(c(k+1894, unique(r$month)), subset(r, select = -c(month)))
    }
  )
  
})
job$shared <- file.path(rh.root, "stationinfo", "USinfo.RData")
job$setup <- expression(
  map = {
    load("USinfo.RData")
    library(plyr)
  }
)
job$parameters <- list(par = par)
job$input <- rhfmt(file.path(root, "Raw", par$dataset), type = "text")
job$output <- rhfmt(file.path(root, par$dataset, "All", "bymonth"), type = "sequence")
job$mapred <- list( 
  mapred.reduce.tasks = 100, #cdh3,4
  mapreduce.job.reduces = 100,  #cdh5
  rhipe_map_buff_size = Nstations #total number of stations
)
job$jobname <- file.path(root, par$dataset, "All", "bymonth")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)



job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r) {
    lapply(1:nrow(map.values[[r]]), function(x) {
      key <- c(as.character(map.values[[r]][x, "station.id"]), map.values[[r]][x, 2:4])
      value <- data.frame(
        resp = map.values[[r]][x, "resp"],
        year = map.keys[[r]][1],
        month = map.keys[[r]][2],
        stringsAsFactors = FALSE
      )
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
    attr(combined, "location") <- unlist(reduce.key)[2:4]
    rhcollect(reduce.key[[1]][1], combined) 
    ##when reduce.key is length 1, then after rhcollect, reduce.key is single value; 
    ##if reduce.key is length>1, then after rhcollect, reduce.key is a list
    #rhcollect(reduce.key, combined)  
  }
)
job$input <- rhfmt(
  file.path(root, par$dataset, "All", "bymonth"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(root, par$dataset, "All", "bystation"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 100,  #cdh3,4
  mapreduce.job.reduces = 100  #cdh5
)
job$readback <- FALSE
job$jobname <- file.path(root, par$dataset, "All", "bystation")
job.mr <- do.call("rhwatch", job)

