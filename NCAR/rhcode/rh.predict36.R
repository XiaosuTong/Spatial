###################################################################
##  Input is the best choice spatial fitting from bymonth.fit,   ##
##  and the key is station.id. Output key is c(station.id, i),   ##
##  where i is the starting time index of time series for one    ##
##  STL prediction of 36 following obs.                          ##
##                                                               ##
##  for a1950 valid=270, tn=576,                                 ##
##  for 100stations, valid=600, tn=1236                          ##
###################################################################
stationSplit <- function(reduce=100, input, type="100stations", tn=1236, valid=600){
  
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      v <- map.values[[r]]
      v <- arrange(v, year, match(month, month.abb))
      v$time <- 0:(tn-1)
      lapply(1:(valid+1), function(i) {
        v <- v[i:(valid + i - 1 + 36),]
        rhcollect(c(map.keys[[r]], i), v)
      })
    })
  })
  job$parameters <- list(
    tn = tn,
    valid = valid
  )
  job$setup <- expression(
    map = {
      library(plyr)
    }
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, type, "bystationSplit"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = reduce,  #cdh3,4
    mapreduce.job.reduces = reduce #cdh5
  )
  job$jobname <- file.path(rh.root, par$dataset, type, "bystationSplit")
  job$mon.sec <- 10
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)

}

###########################################################
## Predict 36 observation ahead using                    ##
## previous 600/270 observations for the 100 stations    ##
## predict36 did not include fc component                ##
###########################################################
predict36 <- function(type, parameter, k, index, valid) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      key <- c(map.keys[[r]][1], group)
      v <- map.values[[r]]
      if(type == "a1950") {
        Index <- which(is.na(v$resp))
        v$resp[Index] <- v$fitted[Index]
      }
      v.raw <- tail(v, 36)
      v.model <- v
      v.model[valid+(1:36), "resp"] <- NA
      if(fc.flag) {
        v.predict <- tail(
          do.call("cbind", 
            stl2(x = v.model$resp, t = v.model$time, n.p = 12, 
              s.window = sw, s.degree = sd, t.window = tw, t.degree = td, 
              fc.window = c(fcw, scw), fc.degree = c(fcd, scd), 
              inner = inner, outer = outer)[c("data","fc")]), 36
        )
        names(v.predict)[grep("fc.fc", names(v.predict))] <- c("fc.first", "fc.second")
      } else {
        v.predict <- tail(stl2(
          x = v.model$resp, 
          t = v.model$time, 
          n.p = 12, s.window = sw, s.degree = sd, t.window = tw, t.degree = td, 
          inner = inner, outer = outer)$data, 36
        )
      }
      value <- cbind(v.raw, v.predict)
      value$lag <- c(1:36)
      if(fc.flag) {
        value <- subset(value, select = -c(data.weights, data.remainder, data.raw, fc.remainder, data.trend, data.sub.labels))
      } else{
        value <- subset(value, select = -c(weights, remainder, raw, sub.labels))
      }
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
  if (parameter[k,"sw"] != "periodic"){
    job$parameters <- list(
      valid = valid, type = type,
      sw = as.numeric(parameter[k,"sw"]), tw = parameter[k,"tw"], sd = parameter[k,"sd"], 
      td = parameter[k,"td"], fcw = parameter[k, "fcw"], fcd = parameter[k, "fcd"],
      scw = parameter[k, "scw"], scd = parameter[k, "scd"], inner = 10, outer = 0, group = k, 
      fc.flag = parameter[k, "fc.flag"]
    )
  }else{
    job$parameters <- list(
      valid = valid, type = type,
      sw = parameter[k,"sw"], tw = parameter[k,"tw"], sd = parameter[k,"sd"], 
      td = parameter[k,"td"], fcw = parameter[k, "fcw"], fcd = parameter[k, "fcd"], 
      scw = parameter[k, "scw"], scd = parameter[k, "scd"], inner = 10, outer = 0, group = k,
      fc.flag = parameter[k, "fc.flag"]
    )
  }
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, type, "bystationSplit"),
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, paste("run", k, sep="")), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 50,  #cdh3,4
    mapreduce.job.reduces = 50,  #cdh5
    rhipe_reduce_buff_size = 10000
  )
  job$jobname <- file.path(rh.root, par$dataset, type, "STLtuning", index, paste("run", k, sep=""))
  job$mon.sec <- 10
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)

}

##################################################################################
## The input files are the all runk files, the output file is by.stgrouplag     ##
## The key in the output is c(station.id, group, lap), in map output, the       ##
## value is one row dataframe which has residual and rest of columns. The       ##
## value from reduce function is the dataframe including all 601/271 replicates ##
## for one c(station.id, group, lap).                                           ##
################################################################################## 
lagResidual <- function(n, index, type){

  job<- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      value <- map.values[[r]]
      if(any(grepl("fc", names(value)))){
        value$residual <- with(value, resp-data.seasonal-fc.first-fc.second)
      }else{
        value$residual <- with(value, resp-seasonal-trend)
      }
      lapply(1:nrow(value), function(i){
        key <- c(map.keys[[r]], value$lag[i])
        rhcollect(key, subset(value, select=-c(lag))[i,])
      })
    })
  })
  job$reduce <- expression(
    pre={
      combined <- data.frame()
    },
    reduce={
      combined <- rbind(combined, do.call(rbind, reduce.values))
    },
    post={
      rhcollect(reduce.key, combined)
    }
  )

  files <- unlist(
    lapply(1:n, function(r) {
      file.path(rh.root, par$dataset, type, "STLtuning", index, paste("run", r, sep=""))
    })
  )
  job$input <- rhfmt(files, type = "sequence")
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "by.stagrouplag"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapreduce.job.reduces = 100,  #cdh5
    mapred.reduce.tasks = 100
  )
  job$jobname <- file.path(rh.root, par$dataset, type, "STLtuning", index, "by.stagrouplag")
  job$combiner <- TRUE
  job$readback <- FALSE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)

}

################################################################################
##  input key is c(station.id, group, lag), input value is data.frame         ##
##  of 601/271 replicates. Output key is 1 which is meaningless, output value ##
##  is the data.frame include 5 quantiles of each c(station.id, group, lag)   ##
################################################################################
lagResidQuan <- function(index, type, reduce=1) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r){
      v <- map.values[[r]]
      value <- data.frame(
        quantiles = as.numeric(quantile(
          x = v$residual, 
          probs = c(0.05, 0.25, 0.5, 0.75, 0.95)
        )),
        type = c(0.05, 0.25, 0.5, 0.75, 0.95)
      )
      value$station.id <- map.keys[[r]][1]
      value$group <- map.keys[[r]][2]
      value$lag <- map.keys[[r]][3]
      if(reduce == 1) {
        rhcollect(1, value)
      } else {
        rhcollect(map.keys[[r]][1], value)
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
      rhcollect(reduce.key, combined)
    }
  )  
  job$setup <- expression(
    map = {
      library(plyr)
  })
  job$parameters <- list(
    parameter = parameter,
    reduce = reduce
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "by.stagrouplag"),
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "residquant"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = reduce
  )
  job$jobname <- file.path(rh.root, par$dataset, type, "STLtuning", index, "residquant")
  job$readback <- FALSE
  job$combiner <- TRUE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)


}


############################################################################
##  input key is c(station.id, group, lag), input value is data.frame     ##
##  with 601/271 replicates. In the map, the mean and std of each         ##
##  c(station.id, group) is calculated over 601/271 replicates. Final     ##
##  output key is c(station.id, group), output value is a vector. First   ##
##  value is the mean over 36 of the mean of residual over 601 for each   ##
##  c(station. id, group). The second value is the mean over 36 of        ##
##  the standard deviation of residual over 601 for each                  ##
##  c(station.id, group).                                                 ##
##  Then the next job rbind all stations all group to be one data.frame   ##
##  by reading c(station.id, group) key-value pairs in "group.absmeanstd" ##
##  to create "overall.absmeanstd".                                       ##
##  The results of this job will be used for dotplot of error vs.         ##
##  station.id superpose on different group.                              ##
############################################################################
StdMean.group <- function(index, type, num) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      v <- map.values[[r]]
      value <- data.frame(
        means = mean(abs(v$residual)),
        std = 1.96*sd(v$residual)
      )
      key <- c(map.keys[[r]][1], map.keys[[r]][2])
      rhcollect(key, value)    
    })
  })
  job$reduce <- expression(
    pre = {
      combined <- data.frame()
    },
    reduce = {
      combined <- rbind(combined, do.call("rbind", reduce.values))
    },
    post = {
      overmean <- sum(combined$means)/num
      overstd <- sum(combined$std)/num
      rhcollect(reduce.key, c(overmean, overstd))
    }
  )
  job$parameters <- list(
    parameter = parameter,
    num = num
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "by.stagrouplag"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "group.absmeanstd"),
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 10
  )
  job$jobname <- file.path(rh.root, par$dataset, type, "STLtuning", index, "group.absmeanstd")
  job$readback <- FALSE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      value <- data.frame(
        station.id = map.keys[[r]][1],
        group = map.keys[[r]][2],
        mean.absmeans = map.values[[r]][1],
        mean.std = map.values[[r]][2],
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
      combined <- rbind(combined, do.call("rbind", reduce.values))
    },
    post = {
      rhcollect(reduce.key, combined)
    }
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "group.absmeanstd"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "overall.absmeanstd"),
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 10
  )
  job$jobname <- file.path(rh.root, par$dataset, type, "STLtuning", index, "overall.absmeanstd")
  job$readback <- FALSE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)

}



#################################################################################
##  The input file is by.stgrouplag, key is c(station.id, group, lap),         ##
##  value is the data.frame of 601/271 replicates.                             ##
##  The output file is the grouplag.absmeanstd, key is 1 if it is              ##
##  100stations, or is station.id for a1950.                                   ##
##  The value is dataframe that each row has mean of abs of error,             ##
##  mean of error, and SD of error for given station, given lap, given group.  ##
##  The results of this job will be used for lineplot of error vs. lag         ##
##  superpose on different groups.                                             ##
#################################################################################
StdMean.grouplag <- function(index, parameter, type) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      v <- map.values[[r]]
      value <- data.frame(
        means = mean(v$residual),
        std = 1.96*sd(v$residual),
        absmeans = mean(abs(v$residual))
      )
      value$station.id <- map.keys[[r]][1]
      value$group <- as.numeric(map.keys[[r]][2])
      value$lag <- as.numeric(map.keys[[r]][3])
      if(type=="100stations") {
        rhcollect(1, value)
      } else {
        rhcollect(map.keys[[r]][1], value)
      }
    })
  })
  job$reduce <- expression(
    pre={
      combined <- data.frame()
    },
    reduce={
      combined <- rbind(combined, do.call(rbind, reduce.values))
    },
    post={
      rhcollect(reduce.key, combined)
    }
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "by.stagrouplag"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "grouplag.absmeanstd"),
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 10
  )
  job$parameters <- list(
    parameter = parameter,
    type = type
  )
  job$jobname <- file.path(rh.root, par$dataset, type, "STLtuning", index, "grouplag.absmeanstd")
  job$readback <- FALSE
  job$combiner <- TRUE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)

}