#######################################################
## Predict 36 observation ahead using                ##
## previous 600 observations for the 100 stations    ##
#######################################################
stationSplit <- function(reduce=100){
  
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      v <- map.values[[r]]
      v$time <- 0:1235
      lapply(1:601, function(i) {
        v <- v[i:(600 + i - 1 + 36),]
        rhcollect(c(map.keys[[r]], i), v)
      })
    })
  })
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "bystation"),
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "bystationSplit"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = reduce,  #cdh3,4
    mapreduce.job.reduces = reduce #cdh5
  )
  job$jobname <- file.path(rh.root, par$dataset, "100stations", "bystationSplit")
  job$mon.sec <- 10
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)

}

predict36 <- function(parameter, k, index) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      key <- c(map.keys[[r]][2], group)
      v <- map.values[[r]]
      v.raw <- tail(v, 36)
      v.model <- v
      v.model[601:636, "resp"] <- NA
      v.predict <- tail(stl2(
        x = v.model$resp, 
        t = v.model$time, 
        n.p = 12, 
        s.window = sw, 
        s.degree = sd, 
        t.window = tw, 
        t.degree = td, 
        fc.window = fcw,
        fc.degree = fcd,
        inner = inner, 
        outer = outer)$data, 36)
      value <- cbind(v.raw, v.predict)
      value$lag <- c(1:36)
      value <- subset(value, select = -c(weights, remainder, raw, sub.labels))
      rhcollect(key, value)
    })
  })
  job$setup <- expression(
    map = {
      library(lattice)
      library(yaImpute)
      library(stl2)
      library(plyr)
    })
  if (parameter[k,"tw"] != "periodic"){
    job$parameters <- list(
      sw = as.numeric(parameter[k,"sw"]), 
      tw = parameter[k,"tw"], 
      sd = parameter[k,"sd"], 
      td = parameter[k,"td"], 
      fcw = parameter[k, "fcw"],
      fcd = parameter[k, "fcd"],
      inner = 10, 
      outer = 0, 
      group = k
    )
  }else{
    job$parameters <- list(
      sw = parameter[k,"sw"], 
      tw = parameter[k,"tw"], 
      sd = parameter[k,"sd"], 
      td = parameter[k,"td"],
      fcw = parameter[k, "fcw"],
      fcd = parameter[k, "fcd"], 
      inner = 10, 
      outer = 0, 
      group = k
    )
  }
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "bystationSplit"),
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "STLtuning", index, paste("run", k, sep="")), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 50,  #cdh3,4
    mapreduce.job.reduces = 50,  #cdh5
    rhipe_reduce_buff_size = 10000
  )
  job$jobname <- paste(par$dataset, index, "predict", k)
  job$mon.sec <- 10
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)

}