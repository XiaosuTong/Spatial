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

#myfun returns a long vector for all observation at one pressure level

myfun <- function(index, air.all){
  tmp <- lapply(1:248, function(rr) {
    as.vector(air.all[1:349, 1:277, index, rr]) #by column
  })
  do.call("c", tmp)
}

getair <- function(myfolder, mydate, myparameter) {
  mync <- open.ncdf( file.path("./tmp", paste(myfolder, mydate, "nc", sep=".")) )
  air.all <- get.var.ncdf(mync, mync$var[[myparameter]])
  lon <- get.var.ncdf(mync, mync$var[["lon"]])
  lat <- get.var.ncdf(mync, mync$var[["lat"]])
  data <- data.frame(
    lon = rep(as.vector(lon), times = 248),
    lat = rep(as.vector(lat), times = 248),
    time = rep(1:248, each = 349*277)
  )
  all.levels<- mlply(
    .data = data.frame(
      index = 1:29
    ),
    .fun  = myfun,
    air.all = air.all
  )
  close.ncdf(mync)
  data <- cbind(data, do.call("cbind", all.levels))
  data
}

########################################################################

par <- list()
par$myfolder <- "air"
par$myparameter <- "air"
par$machine <- "gacrux"
if(par$myfolder == "air") {
  par$address <- "ftp://ftp.cdc.noaa.gov/Datasets/NARR/pressure/"
}
if(par$machine == "gacrux") {
  rh.datadir <- "/ln/tongx/Spatial/NARR"
}
job <- list()
job$map <- expression({
  msys <- function(on){
    system(sprintf("wget  %s --directory-prefix ./tmp 2> ./errors", on))
    if(length(grep("(failed)|(unable)", readLines("./errors"))) > 0){
      stop(paste(readLines("./errors"), collapse="\n"))
    }
  }
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
        rhcollect(key, k[, !(names(k) %in% "time")])
      }
    )
  })
})
job$setup <- expression(
  map = {
    suppressMessages(library(plyr,lib.loc = lib.loc))
    suppressMessages(library(ncdf,lib.loc = lib.loc))
  }
)
job$parameters <- list(
  par = par,
  myfun = myfun,
  getair = getair,
  lib.loc = file.path(path.expand("~"), "R_LIBS")
)
job$input <- c(2, 1) 
job$output <- rhfmt(
  file.path(rh.datadir, par$myfolder, "bytime"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 72
#  rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 5
job$copyFiles <- TRUE
job$jobname <- file.path(rh.datadir, par$myfolder, "bytime")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)