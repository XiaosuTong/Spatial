install.packages("ncdf")
library(ncdf)

# Example of how the paste function works:
getparameternames <- function(myfolder, mydate) {
  mync <- open.ncdf( paste(myfolder, mydate, "nc", sep=".") )
  mycount <- mync$nvars
  mynames <- sapply( c(1:mycount), function(x) {mync$var[[x]]$name} )
  close.ncdf(mync)
  mynames
}

# General function for getting the times from SATURN 03,
# in a given folder (myfolder), for a given date (mydate)

gettimes <- function(myfolder, mydate) {
  mync <- open.ncdf( paste(myfolder, mydate, "nc", sep=".") )
  mytimes <- get.var.ncdf(mync,"time")
  close.ncdf(mync)
  mytimes
}

  
# General function for getting the data from SATURN 03,
# in a given folder (myfolder), for a given date (mydate), for a given parameter (myparameter)
  
getdata <- function(myfolder, mydate, myparameter) {
  mync <- open.ncdf( paste("./tmp/", myfolder, mydate, "nc", sep=".") )
  myvec <- get.var.ncdf(mync, mync$var[[myparameter]])
  close.ncdf(mync)
  myvec
}

# myfun returns a long vector for all observation at one pressure level

myfun.air <- function(index, time.len, air.all){
  tmp <- lapply(1:time.len, function(rr) {
    as.vector(air.all[1:349, 1:277, index, rr]) #by column
  })
  do.call("c", tmp)
}

# Get the air variable in the Pressure data
getair <- function(myfolder, mydate, myparameter) {
  mync <- open.ncdf( file.path("./tmp", paste(myfolder, mydate, "nc", sep=".")) )
  air.all <- get.var.ncdf(mync, mync$var[[myparameter]])
  lon <- get.var.ncdf(mync, mync$var[["lon"]])
  lat <- get.var.ncdf(mync, mync$var[["lat"]])
  time.len <- dim(air.all)[4]
  data <- data.frame(
    lon = rep(as.vector(lon), times = time.len),
    lat = rep(as.vector(lat), times = time.len),
    time = rep(1:time.len, each = 349*277)
  )
  all.levels<- mlply(
    .data = data.frame(
      index = 1:29
    ),
    .fun  = myfun.air,
    time.len = time.len,
    air.all = air.all
  )
  close.ncdf(mync)
  data <- cbind(data, do.call("cbind", all.levels))
  data
}


myfun.sfc.air <- function(time.len, air.all){
  tmp <- lapply(1:time.len, function(rr) {
    as.vector(air.all[1:349, 1:277, rr]) #by column
  })
  do.call("c", tmp)
}

# Get the surface air temperature in Monolevel data
getsfc.air <- function(myfolder, mydate, myparameter) {
  mync <- open.ncdf( file.path("./tmp", paste(myfolder, mydate, "nc", sep=".")) )
  air.all <- get.var.ncdf(mync, mync$var[[myparameter]])
  lon <- get.var.ncdf(mync, mync$var[["lon"]])
  lat <- get.var.ncdf(mync, mync$var[["lat"]])
  time.len <- dim(air.all)[3]
  data <- data.frame(
    lon = rep(as.vector(lon), times = time.len),
    lat = rep(as.vector(lat), times = time.len),
    time = rep(1:time.len, each = 349*277)
  )
  data$air.sfc <- myfun.sfc.air(time.len, air.all)
  close.ncdf(mync)
  data
}

msys <- function(on){
  system(sprintf("wget  %s --directory-prefix ./tmp 2> ./errors", on))
  if(length(grep("(failed)|(unable)", readLines("./errors"))) > 0){
    stop(paste(readLines("./errors"), collapse="\n"))
  }
}
########################################################################

par <- list()
par$machine <- "gacrux"
if(par$machine == "gacrux") {
  rh.datadir <- "/ln/tongx/Spatial/NARR"
}

#par$myfolder <- "air"
#par$myparameter <- "air"
par$myfolder <- "air.sfc"
par$myparameter <- "air"
if(par$myfolder == "air") {
  par$address <- "ftp://ftp.cdc.noaa.gov/Datasets/NARR/pressure/"
}else if(par$myfolder == "air.sfc") {
  par$address <- "ftp://ftp.cdc.noaa.gov/Datasets/NARR/monolevel/"
}

job.air <- list()
job.air$map <- expression({
  lapply(map.values, function(r) {
    year <- 1978 + ceiling(r/12)
    month <- r - (year-1979)*12
    key <- paste(year, sprintf("%02d", month), sep = "")
    on <- sprintf(paste(par$address, par$myfolder, ".%s.nc", sep = ""), key)
    rhstatus(sprintf("Downloading %s", key))
    msys(on)
    rhstatus(sprintf("Downloaded %s", key))
    rhcounter("FILES", key, 1)
    value.all <- getair(par$myfolder, key, par$myparameter)
    system("rm ./tmp/*.nc")
    rhcounter("getair", key, 1)
    d_ply(
      .data = value.all,
      .variable = "time",
      .fun = function(k){
        key <- paste(key, sprintf("%03d", unique(k$time)), sep = "")
        value <- k[, !(names(k) %in% "time")]
        attr(value, "file") <- par$myfolder
        rhcollect(key, value)
      }
    )
  })
})
job.air$setup <- expression(
  map = {
    suppressMessages(library(plyr, lib.loc = lib.loc))
    suppressMessages(library(ncdf, lib.loc = lib.loc))
  }
)
job.air$parameters <- list(
  par = par,
  msys = msys,
  myfun.air = myfun.air,
  getair = getair,
  lib.loc = file.path(path.expand("~"), "R_LIBS")
)
job.air$input <- c(24, 12) 
job.air$output <- rhfmt(
  file.path(rh.datadir, par$myfolder, "bytime"), 
  type = "sequence"
)
job.air$mapred <- list(
  mapred.reduce.tasks = 72,
  mapred.task.timeout = 0
#  rhipe_reduce_buff_size = 10000
)
job.air$mon.sec <- 5
job.air$copyFiles <- TRUE
job.air$jobname <- file.path(rh.datadir, par$myfolder, "bytime")
job.air$readback <- FALSE


####################################################################
##
####################################################################
bound <- data.frame(
  month = 0:12,
  day = c(
    0,31,59,90,120,151,181,
    212,243,273,304,334,365
  )
)
bound.leap <- bound
bound.leap$day <- (bound.leap$day + c(0,0,rep(1,11)))*8
bound$day <- bound$day*8
leap <- c(1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2012)

job.sfc.air <- list()
job.sfc.air$map <- expression({
  lapply(map.values, function(r) {
    year <- 1978 + r
    on <- sprintf(paste(par$address, par$myfolder, ".%s.nc", sep = ""), year)
    rhstatus(sprintf("Downloading %s", year))
    msys(on)
    rhstatus(sprintf("Downloaded %s", year))
    rhcounter("FILES", year, 1)
    value.all <- getsfc.air(par$myfolder, year, par$myparameter)
    system("rm ./tmp/*.nc")
    rhcounter("getair", year, 1)
    d_ply(
      .data = value.all,
      .variable = "time",
      .fun = function(k){
        t <- unique(k$time)
        if(year %in% leap) {
          month <- max(bound.leap[bound.leap$day >= t, 1])
          ##find the month for a given time
          time <- t - bound.leap[month, "day"] 
          ##find the time in month for a given time
        }else {
          month <- min(bound[bound$day >= t, 1])  
          ##find the month for a given time
          time <- t - bound[month, "day"] 
          ##find the time in month for a given time
        }
        key <- paste(year, sprintf("%02d", month), sprintf("%03d", time), sep = "") 
        value <- k[, !(names(k) %in% "time")]
        attr(value, "file") <- par$myfolder         
        rhcollect(key, value)
      }
    )
  rm(list = ls())
  gc()  
  })
})
job.sfc.air$setup <- expression(
  map = {
    suppressMessages(library(plyr,lib.loc = lib.loc))
    suppressMessages(library(ncdf,lib.loc = lib.loc))
  }
)
job.sfc.air$parameters <- list(
  par        = par,
  msys       = msys,
  myfun.air  = myfun.sfc.air,
  getair     = getsfc.air,
  bound      = bound,
  bound.leap = bound.leap,
  leap       = leap,
  lib.loc    = file.path(path.expand("~"), "R_LIBS")
)
job.sfc.air$input <- c(24, 12) 
job.sfc.air$output <- rhfmt(
  file.path(rh.datadir, par$myfolder, "bytime"), 
  type = "sequence"
)
job.sfc.air$mapred <- list(
  mapred.reduce.tasks = 100,
  mapred.task.timeout = 0,
  rhipe_map_buff_size = 1
)
job.sfc.air$mon.sec <- 10
job.sfc.air$copyFiles <- TRUE
job.sfc.air$jobname <- file.path(rh.datadir, par$myfolder, "bytime")
job.sfc.air$readback <- FALSE


job.mr <- do.call("rhwatch", job.sfc.air)