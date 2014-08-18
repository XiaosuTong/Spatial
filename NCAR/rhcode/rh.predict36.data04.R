####################################################################################
##
## Predict 36 observation ahead using previous 600 observations for the 100 stations
## Different parameter set up are used, predict1 to predict9 are the results for
## different setting.
## tmax.100stations.RData which has tmax.stations dataframe is already on HDFS, for
## each map.keys/map.values, this dataframe will be called by using job$setup and
## job$shared.
## stations information is in USinfo.RData which is also on HDFS.
##
####################################################################################
#Emultiple1 is the experiment5
#final else, which is Emultiple is the experiment6.
#In data04 this script, fc.degree and fc.window are varied.
dataset <- "tmax"
index <- "Emultiple"

#Load the station.id for the 100 stations for Tmax
load(file.path(local.datadir, paste(dataset,"div.stations.RData", sep="")))
load(file.path(local.datadir, "stations.RData"))
if(index == "Emultiple1"){
	parameter <- list(
    list(
      sw = "periodic", 
      tw = 1855, 
      sd = 1, 
      td = 1, 
      fcw = 1855, 
      fcd = 1, 
      ssw = 121, 
      ssd = 2, 
      multiple = TRUE
    )
  )
}else if(index == "Emultiple") {
	parameter <- list(
    run1 = list(
      sw = "periodic", 
      tw = 1855, 
      sd = 1, 
      td = 1, 
      fcw = 1855, 
      fcd = 1, 
      ssw = 241, 
      ssd = 1, 
      multiple = TRUE
    ),
		run2 = list(
      sw = "periodic", 
      tw = 241, 
      sd = 1, 
      td = 1, 
      multiple = FALSE
    )
  )
}
######################################################################################
## The output runk is from the kth parameter setting. The key is c(station.id, group), 
## which group is from 1 to 9. For each parameter setting, each station, there will
## be 601 runs for 36 prediction, so the total key/value pairs are 60100. The value is 
## the corresponding dataframe which has stl fitting for coming 36 months based on 
## previous 600 months.
######################################################################################
for (k in 1:length(parameter)) {
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      i <- ceiling(map.keys[[r]]/601)
      j <- map.keys[[r]] - (i - 1)*601
      v <- subset(
        get(paste(dataset, ".100stations", sep = "")), 
        station.id == unique(get(paste(dataset, ".100stations", sep=""))$station.id)[i]
      )
      key <- c(
        levels(get(paste(dataset, ".100stations", sep = ""))$station.id)[i],
        k
      )
      v$month <- factor(
        v$month, 
        levels = c(
          "Jan", "Feb", "Mar", "Apr", "May", "June",
          "July","Aug", "Sep", "Oct", "Nov", "Dec"
        )
      )
      v <- v[order(v$year, v$month), ]
      v$time <- 0:1235
      v <- v[j:(600 + j - 1 + 36),]
      v.raw <- tail(v, 36)
      v.model <- v
      v.model[601:636, c(dataset)] <- NA
      if(flag) { 
        v.predict <- tail(
          do.call("cbind", 
            stl2(
              x = v.model[,c(dataset)], 
              t = v.model$time, 
              n.p = 12, 
              s.window = sw, 
              s.degree = sd, 
              t.window = tw, 
              t.degree = td, 
              fc.window = c(fcw, ssw), 
              fc.degree = c(fcd, ssd), 
              inner = inner, 
              outer = outer
            )[c("data","fc")]
          ), 
          36
        )
      }else{
        v.predict <- tail(
          stl2(
            x = v.model[,c(dataset)], 
            t = v.model$time, 
            n.p = 12, 
            s.window = sw, 
            s.degree = sd, 
            t.window = tw, 
            t.degree = td, 
            inner = inner, 
            outer = outer
          )$data, 
          36
        )
      }
      value <- cbind(v.raw, v.predict)
      names(value)[grep("fc.fc", names(value))] <- c("fc.trend", "fc.seasonal")
      value$lap <- c(1:36)
      elev <- unique(value$elev)
      lon <- unique(value$lon)
      lat <- unique(value$lat)
      name <- as.character(unique(value$station.name))
      if(flag) {
        value <- subset(
          value, 
          select = -c(
            data.weights, data.remainder, data.raw, 
            data.sub.labels, lon, lat, station.name, elev
          )
        )
      }else{
        value <- subset(
          value, 
          select = -c(
            weights, remainder, raw, sub.labels, 
            lon, lat, station.name, elev
          )
        )
      } 
      attributes(value)$elev  <- elev
      attributes(value)$lon <- lon
      attributes(value)$lat <- lat
      attributes(value)$station.name <- name
      rhcollect(key, value)
    })
  })
  job$setup <- expression(
    map = {
      library(lattice)
      library(yaImpute, lib.loc = lib.loc)
      library(stl2, lib.loc = lib.loc)
#     suppressMessages(library(bit,lib.loc=lib.loc))
      load(paste(dataset, ".100stations.RData", sep=""))
    }
  )
  job$shared <- c(
    file.path(rh.datadir,dataset,"100stations","Rdata", paste(dataset, ".100stations.RData", sep=""))
  )
  if (parameter[[k]]$multiple){
    job$parameters <- list(
      sw = parameter[[k]]$sw, 
      tw = parameter[[k]]$tw, 
      sd = parameter[[k]]$sd, 
      td = parameter[[k]]$td, 
      fcw = parameter[[k]]$fcw, 
      fcd = parameter[[k]]$fcd, 
      ssw = parameter[[k]]$ssw, 
      ssd = parameter[[k]]$ssd, 
      inner = 10, 
      outer = 0, 
      dataset = dataset, 
      flag = parameter[[k]]$multiple
    )
  }else{
    job$parameters <- list(
      sw = parameter[[k]]$sw, 
      tw = parameter[[k]]$tw, 
      sd = parameter[[k]]$sd, 
      td = parameter[[k]]$td, 
      inner = 10, 
      outer = 0, 
      dataset = dataset, 
      flag = parameter[[k]]$multiple
    )
  }
  job$input <- c(60100, 220)
  job$output <- rhfmt(
    file.path(rh.datadir, dataset, "100stations","sharepredict", index, paste("run", k, sep="")), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 72, 
    rhipe_reduce_full_size = 10000
  )
  job$jobname <- paste(dataset, "predict", k)
  job$mon.sec <- 10
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)
}


#################################################################################
## The input files are the all runk files, the output file is 36.lap.station.
## The key in the output is c(group, lap), in map output, the value is one row
## dataframe which has residual and rest of columns. The value from reduce function
## is the dataframe including all the stations for one c(group, lap).
################################################################################# 
job<- list()
#After 36.lap.station file, all information about stations, like lon, lat, elev
#is gone. For station information, load the USinfo.RData on HDFS.
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		lapply(1:dim(map.values[[r]])[1], function(i){
		  value <- map.values[[r]][i, ]
		  if(any(grepl("fc", names(value)))){
        value$residual <- 
          value[,c(dataset)] - value$data.seasonal - value$fc.trend - value$fc.seasonal
		  }else{
        value$residual <- value[,c(dataset)] - value$seasonal - value$trend
		  }
		  value$station.id <- map.keys[[r]][1]
		  key <- c(map.keys[[r]][2], value$lap)
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

files <- unlist(lapply(1:length(parameter), function(r) {
	file.path(rh.datadir, dataset, "100stations", "sharepredict", index, paste("run", r, sep="")) 
}))
job$parameters <- list(
  dataset = dataset
)
job$input <- rhfmt(files, type="sequence")
job$output <- rhfmt(
  file.path(rh.datadir, dataset,"100stations","sharepredict",index,"36.lap.station"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 72
)
job$jobname <- paste(dataset, "36 each lap")
job$readback <- FALSE
job$mon.sec <- 10
job.mr <- do.call("rhwatch", job)


#################################################################################
## The input file is the 36.lap.station, key is c(group, lap), value is the 
## dataframe for all 100 stations information and fitting for that c(group, lap).
## The output file is the 36.lap.quantile, key is the quantile from c(0.05, 0.25, 
## 0.5, 0.75, 0.95), the value is the dataframe includes residual, and lap, 
## station.id, group. 
#################################################################################
job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r){
    v <- map.values[[r]]
    tmp <- ddply(
      .data = v, 
      .variables = "station.id",
      .fun = summarise,
      lap = lap[1:5],
      quantiles = quantile(
        x = residual, 
        probs=c(0.05, 0.25, 0.5, 0.75, 0.95)
      ),
      type = c(0.05, 0.25, 0.5, 0.75, 0.95)
    )
    tmp$group <- rep(map.keys[[r]][1], dim(tmp)[1])
    lapply(1:dim(tmp)[1], function(i){
      key <- tmp[i, 4]
      value <- tmp[i, -4]
#	    for(j in 1:dim(parameter)[1]){
#	      if(value$group == j){ 
#         tmp$fcw <- parameter[j,5]
#         tmp$ssw <- parameter[j,7]
#	      }
# 		}
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
job$setup <- expression(
  map = {
    library(plyr)
  }
)
job$parameters <- list(
  parameter = parameter
)
job$input <- rhfmt(
  file.path(rh.datadir, dataset, "100stations", "sharepredict", index, "36.lap.station"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.datadir, dataset,"100stations","sharepredict",index,"36.lap.quantile"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 72
)
job$jobname <- paste(dataset, "36 each lap")
job$readback <- FALSE
job$mon.sec <- 10
job.mr <- do.call("rhwatch", job)


####################################################################################
## The input file is the 36.lap.station, key is c(group, lap), value is the 
## dataframe for all 100 stations information and fitting for that c(group, lap).
## The output file is the over.absmeanstd.lap.station, key is c(station.id, sw, tw), 
## value is the overall mean of mean(abs(residual)) over 36 lags.
## Get the dotplot of overmean against sw or tw conditional on station, as well as
## the scatter plot matrix of overmean.
####################################################################################
job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r){
    v <- map.values[[r]]
    tmp <- ddply(
      .data = v,
      .variables="station.id",
      .fun = summarise,
      lap = lap[1],
      means = mean(abs(residual)),
			std = 1.96*sd(residual)
    )
    tmp$group <- map.keys[[r]][1]
#   for(j in 1:dim(parameter)[1]){
#     if(unique(tmp$group) == j){
#       tmp$fcw <- parameter[j,5]
#       tmp$ssw <- parameter[j,7]
#     }
#   }
    myfun <- function(data){
#     rhcollect(c(data$station.id,data$fcw,data$ssw), c(data$means, data$std))
			rhcollect(c(data$station.id,data$group), c(data$means, data$std))
    }
    d_ply(
      .data = tmp,
      .variables = "station.id",
      .fun = myfun
    )
  })
})
job$reduce <- expression(
  pre = {
    total.mean <- 0
		total.std <- 0
  },
  reduce = {
    total.mean <- total.mean + sum(unlist(lapply(reduce.values, function(r){r[1]})))
		total.std <- total.std + sum(unlist(lapply(reduce.values, function(r){r[2]})))
  },
  post = {
    overmean <- total.mean/36
		overstd <- total.std/36
    rhcollect(reduce.key, c(overmean, overstd))
  }
)
job$setup <- expression(
  map = {
    library(plyr)
  }
)
job$parameters <- list(
  parameter = parameter
)
job$input <- rhfmt(
  file.path(rh.datadir, dataset,"100stations","sharepredict",index,"36.lap.station"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.datadir, dataset,"100stations","sharepredict",index,"over.absmeanstd.lap.station"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 72
)
job$jobname <- paste(dataset, "mean over lag")
job$readback <- FALSE
job$mon.sec <- 10
job.mr <- do.call("rhwatch", job)


##################################################################################
## The input file is the 36.lap.station, key is c(group, lap), value is the 
## dataframe for all 100 stations information and fitting for that c(group, lap).
## The output file is the absmeanstd.lap.station, key is 1 which is meaningless.
## The value is dataframe that each row has mean of abs of error, mean of error, 
## and SD of error for given station, given lap, given group.
##################################################################################
job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r){
    v <- map.values[[r]]
    tmp <- ddply(
      .data = v,
      .variables = "station.id",
      .fun = summarise,
      lap = lap[1],
      absmeans = mean(abs(residual)),
		 	means = mean(residual),
		 	std = 1.96*sd(residual)
    )
    tmp$group <- map.keys[[r]][1]
#   for(j in 1:dim(parameter)[1]){
#     if(unique(tmp$group) == j){
#       tmp$fcw <- parameter[j,5]
#       tmp$ssw <- parameter[j,7]
#     }
#   }
    rhcollect(1, tmp)
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
  }
)
job$input <- rhfmt(
  file.path(rh.datadir,dataset,"100stations","sharepredict",index,"36.lap.station"), 
  type = "sequence"
)
job$output <- rhfmt(
  file.path(rh.datadir, dataset,"100stations","sharepredict",index,"absmeanstd.lap.station"), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 72
)
job$jobname <- "mean of error"
job$readback <- FALSE
job$mon.sec <- 10
job.mr <- do.call("rhwatch", job)

