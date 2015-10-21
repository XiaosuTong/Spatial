##for a1950 valid=270
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

#######################################################
## Predict 36 observation ahead using                ##
## previous 600 observations for the 100 stations    ##
## predict36 did not include fc component            ##
#######################################################
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

#################################################################################
## The input files are the all runk files, the output file is by.stgrouplag
## The key in the output is c(station.id, group, lap), in map output, the value is one row
## dataframe which has residual and rest of columns. The value from reduce function
## is the dataframe including all 601 replicates for one c(station.id, group, lap).
################################################################################# 
lagResidual <- function( n, index){

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
      file.path(rh.root, par$dataset, "100stations", "STLtuning", index, paste("run", r, sep=""))
    })
  )
  job$input <- rhfmt(files, type = "sequence")
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "STLtuning", index, "by.stgrouplag"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 100
  )
  job$jobname <- file.path(rh.root, par$dataset, "100stations", "STLtuning", index, "by.stgrouplag")
  job$combiner <- TRUE
  job$readback <- FALSE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)

}

############################################################################
##  input key is c(station.id, group, lag), input value is data.frame
##  of 601 replicates. Output key is 1 which is meaningless, output value
##  is the data.frame include 5 quantiles of each c(station.id, group, lag)
############################################################################
lagResidQuan <- function(index) {

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
  job$setup <- expression(
    map = {
      library(plyr)
  })
  job$parameters <- list(
    parameter = parameter
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "STLtuning", index, "by.stgrouplag"),
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "STLtuning", index, "residquant"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 1
  )
  job$jobname <- file.path(rh.root, par$dataset, "100stations", "STLtuning", index, "residquant")
  job$readback <- FALSE
  job$combiner <- TRUE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)


}


############################################################################
##
##  input key is c(station.id, group, lag), input value is data.frame
##  of 601 replicates. Output key is c(station.id, group), output value
##  is a vector. First value is the mean over 36 of the mean of residual 
##  over 601 for each station. The second value is the mean over 36 of 
##  the standard deviation of residual over 601 for each station.
##
############################################################################
StdMean.group <- function(index) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      v <- map.values[[r]]
      value <- data.frame(
        means = mean(abs(v$residual)),
        std = 1.96*sd(v$residual)
      )
      key <- map.keys[[r]][1:2]
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
      overmean <- sum(combined$means)/36
      overstd <- sum(combined$std)/36
      rhcollect(reduce.key, c(overmean, overstd))
    }
  )
  job$parameters <- list(
    parameter = parameter
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "STLtuning", index, "by.stgrouplag"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "STLtuning", index, "group.absmeanstd"),
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 10
  )
  job$jobname <- file.path(rh.root, par$dataset, "100stations", "STLtuning", index, "group.absmeanstd")
  job$readback <- FALSE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)

}


##################################################################################
##
##  The input file is by.stgrouplag, key is c(station.id, group, lap), 
##  value is the data.frame of 601 replicates.
##  The output file is the absmeanstd, key is 1 which is meaningless.
##  The value is dataframe that each row has mean of abs of error, mean of error, 
##  and SD of error for given station, given lap, given group.
##
##################################################################################
StdMean.grouplag <- function(index, parameter) {

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
      rhcollect(1, value)
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
    file.path(rh.root, par$dataset, "100stations", "STLtuning", index, "by.stgrouplag"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "100stations", "STLtuning", index, "grouplag.absmeanstd"),
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 10
  )
  job$parameters <- list(
    parameter = parameter
  )
  job$jobname <- file.path(rh.root, par$dataset, "100stations", "STLtuning", index, "grouplag.absmeanstd")
  job$readback <- FALSE
  job$combiner <- TRUE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)

}