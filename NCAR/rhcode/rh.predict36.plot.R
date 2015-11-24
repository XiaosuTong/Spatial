####################################################################
##  Order the stations based on the deviation from the normal     ##
##  distribution. For each station, given lap, given group,       ##
##  the distribution of 601/271 prediction error is compared      ##
##  with normal distribution. The sum of abs deviation over all   ##
##  lag for one station under given group is calculated.          ##
##  The order of stations under given group is saved as           ##
##  group.divorderstations                                        ##
####################################################################
QQDivFromNormal <- function(index="E1", type="a1950") {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r){
      if(map.keys[[r]][1] %in% sample.a1950$station.id) {
        v <- map.values[[r]]
        key <- c(map.keys[[r]][1], map.keys[[r]][2])
        v$group <- map.keys[[r]][1]
        yy <- quantile(v$residual, c(0.25, 0.75))
        xx <- qnorm(c(0.25, 0.75))
        r <- diff(yy)/diff(xx)
        x <- qnorm(ppoints(length(v$residual)))
        y <- r*x + yy[1] - xx[1]*r
        div <- sum(abs(sort(v$residual) - y))
        rhcollect(key, div)
        rhcounter("Station","sample", 1)
      }
    })
  })
  job$reduce <- expression(
    pre = {
      total <- 0
    },
    reduce = {
      total <- total + sum(unlist(reduce.values))
    },
    post = {
      rhcollect(reduce.key, total)
      rhcounter("Station", "reduce", 1)
    }
  )
  job$shared <- c(file.path(rh.root, par$dataset, type, "Rdata", "sample.a1950.RData"))
  job$setup <- expression(
    map = {
      library(plyr)
      load("sample.a1950.RData")
    }
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "by.stagrouplag"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "divorderstations"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 72
  )
  job$jobname <- file.path(rh.root, par$dataset, type, "STLtuning", index, "divorderstations")
  job$readback <- FALSE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r){
      value <- data.frame(
        station.id = map.keys[[r]][1], 
        leaf = with(sample.a1950, leaf[which(station.id == map.keys[[r]][1])]),
        div = as.numeric(map.values[[r]]), 
        stringsAsFactors=FALSE
      )
      rhcollect(map.keys[[r]][2], value)
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
      combined <- arrange(combined, div)
      rhcollect(reduce.key, combined)
    }
  )
  job$shared <- c(file.path(rh.root, par$dataset, type, "Rdata", "sample.a1950.RData"))
  job$setup <- expression(
    map = {load("sample.a1950.RData")},
    reduce = {library(plyr)}
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "divorderstations"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "group.divorderstations"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 10
  )
  job$jobname <- file.path(rh.root, par$dataset, type, "STLtuning", index, "group.divorderstations")
  job$readback <- FALSE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)

}


#################################################################################
##  For each group, get a normal quantile plot for each station conditional    ##
##  on lag. Under given group, stations are ordered by diviation from normal   ##
##  distribution.                                                              ##
##  There is an issue that the combined in reduce function is over 128MB       ##
##  which will give an error for out of memory. So we only plot in the         ##
##  reduce step, no data is saved in this step. The station is ordered by      ##
##  the amount of divation from normal, and the order is called from           ##
##  group.divorderstations.                                                    ##
#################################################################################
## plotEngine is the plotting function called in reduce function for each group.
## orderdf is the data.frame contains station order and station leaf.
plotEngine <- function(data, dataset, key, orderdf) { 

  data$lag <- as.numeric(data$lag)
  trellis.device(
    device = postscript, 
    file = file.path(".", "tmp", paste("QQ.error", dataset, "group", key, "ps", sep= ".")), 
    color = TRUE, 
    paper = "legal"
  )
  for(i in orderdf$station.id){
    b <- qqmath( ~ residual | lag
      , data = subset(data, station.id==i)
      , xlab = list(label = "Unit normal quantile", cex = 1.5)
      , ylab = list(label = ylab, cex = 1.5)
      , sub = paste("Station", i, "from cell", with(orderdf, leaf[which(station.id==i)]))
      , type = "p"
      , aspect = 1
      , col = "red"
      , layout = c(9,4)
      , scale = list(
          x = list(cex=1.2, at = seq(-2,2,2)),
          y = list(cex=1.2, at = seq(-10,10,5))
        )
      , pch=16
      , cex=0.3
      , prepanel = prepanel.qqmathline
      , panel = function(x,...) {
          panel.qqmathline(x, distribution = qnorm)
          panel.qqmath(x,...)
        }
    )
    print(b)
  }
  dev.off()

}
QQstationlag <- function(param=parameter, index="E1", type="a1950") {

  stationOrder <- rhread(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "group.divorderstations")
  )
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r){
      if(map.keys[[r]][1] %in% sample.a1950$station.id) {
        value <- subset(map.values[[r]], select = c(residual))
        value$lag <- map.keys[[r]][3]
        value$station.id <- map.keys[[r]][1]
        group <- as.numeric(map.keys[[r]][2])
        value$sw <- param[group, "sw"]
        value$tw <- param[group, "tw"]
        rhcollect(group, value)
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
      I <- which(unlist(lapply(div.stations, "[[", 1))==reduce.key)
      plotEngine(data=combined, dataset=dataset, key=reduce.key, orderdf=div.stations[[I]][[2]])
    }
  )
  job$shared <- c(file.path(rh.root, par$dataset, type, "Rdata", "sample.a1950.RData"))
  job$setup <- expression(
    map = {load("sample.a1950.RData")},
    reduce = {
      load("sample.a1950.RData")
      library(lattice)
    }
  )
  job$parameters <- list(
    param = param, 
    ylab = ylab,
    dataset = par$dataset,
    plotEngine = plotEngine,
    div.stations = stationOrder
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "by.stagrouplag"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "meanerrorqqplot"), 
    type="sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 10
  )
  job$jobname <- file.path(rh.root, par$dataset, type, "STLtuning", index, "meanerrorqqplot")
  job$readback <- FALSE
  job$mon.sec <- 10
  job$copyFiles <- TRUE
  job.mr <- do.call("rhwatch", job)  

}

#####################################################################################
## The input file is the grouplag.absmeanstd, key is 1 which is meaningless,       ##
## value is the dataframe which has absmean, mean, and SD of prediction error for  ##
## a given station, give lap, and given group.                                     ##
## Get the scatter plot of mean of ABS(ERROR) cross 601 replicates vs lag          ##
## conditional on group for each station.                                          ##
#####################################################################################
subsetStations <- function(type, index) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      if(map.keys[[r]][1] %in% sample.a1950$station.id) {
        rhcollect(1, map.values[[r]])
        rhcounter("Station","sample", 1)
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
  job$shared <- c(file.path(rh.root, par$dataset, type, "Rdata", "sample.a1950.RData"))
  job$setup <- expression(
    map = {load("sample.a1950.RData")},
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "grouplag.absmeanstd"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "grouplagSample.absmeanstd"), 
    type="sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 1
  )
  job$jobname <- file.path(rh.root, par$dataset, type, "STLtuning", index, "grouplagSample.absmeanstd")
  job$readback <- FALSE
  job$mon.sec <- 15
  job.mr <- do.call("rhwatch", job)  

}

errorVsLag <- function(type, index, var, target) {

  ## grouplagSample.absmeanstd is one key-value pair which includes all 128 stations mean of given group given lag
  rst <- rhread(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "grouplagSample.absmeanstd")
  )[[1]][[2]]  
  
  rst <- cbind(rst, parameter[c(rst$group), var]) 

  if(("periodic" %in% parameter$sw) & ("sw" %in% var)){
    rst$sw <- factor(rst$sw, 
      levels=c(sort(as.numeric(unique(rst$sw)[which(unique(rst$sw)!="periodic")])), "periodic")
    )
    rst[,var[which(var != "sw")]] <- as.factor(rst[,var[which(var != "sw")]])
  }else{
    rst[,var[1]] <- as.factor(rst[,var[1]])
    rst[,var[2]] <- as.factor(rst[,var[2]])
  }
  
  rst <- arrange(rst, station.id, group, lag)

  if(target == "means"){
    ylab <- "Mean of Error"
  } else if(target == "absmeans") {
    ylab <- "Mean of Absolute Value of Prediction Error"
  } else {
    ylab <- "Standard Deviation of Prediction Error"
  }

  rhload(file.path(rh.root, par$dataset, type, "Rdata", "sample.a1950.RData"))
  num <- with(rst, length(levels(get(var[2]))))
  trellis.device(
    device = postscript, 
    file   = file.path(local.root, "output", paste(par$dataset, target,"vs.lag", var[1],"ps", sep=".")), 
    color  = TRUE, 
    paper  = "letter"
  )
  for(i in 1:128) {
    b <- xyplot( get(target) ~ lag | get(var[1])
      , data = subset(rst, station.id == with(sample.a1950, station.id[which(leaf==i)]))
      , sub = paste("Station", with(sample.a1950, station.id[which(leaf==i)]), "from cell", i)
      , xlab = list(label = "Lag", cex = 1.5)
      , ylab = list(label = ylab, cex = 1.5)
      , groups = get(var[2])
      , key = list(
          type = "l", 
          text = list(label=paste(var[2],"=",with(rst, levels(get(var[2]))), sep="")),  
          lines = list(lwd=2, col=col[1:num]), 
          columns = num
        )
      , type = "b"
      , lwd = 1.5
      , cex = 0.7
      , scales = list(
          x = list(at=seq(from=0, to=36, by=6), cex=1.2), 
          y = list(relation="sliced", cex=1.2)
        )
      , layout = c(num,1)
      , panel = function(x,y,...) {
          panel.xyplot(x,y,...)
          panel.abline(h=0, v=seq(0,36, by=12), color="black", lty=1, lwd=0.5)
        }
    )
    print(b)
  }
  dev.off()
  num <- with(rst, length(levels(get(var[1]))))
  trellis.device(
    device = postscript, 
    file   = file.path(local.root, "output", paste(par$dataset, target, "vs.lag", var[2],"ps", sep=".")), 
    color = TRUE, 
    paper="letter"
  )
  for(i in 1:128) {
    b <- xyplot( get(target) ~ lag | get(var[2])
      , data = subset(rst, station.id == with(sample.a1950, station.id[which(leaf==i)]))
      , sub = paste("Station", with(sample.a1950, station.id[which(leaf==i)]), "from cell", i)
      , xlab = list(label = "Lag", cex = 1.5),
      , ylab = list(label = ylab, cex = 1.5),
      , groups = get(var[1])
      , key = list(
          type = "l", 
          text = list(label=paste(var[1],"=",with(rst, levels(get(var[1]))), sep="")),  
          lines = list(lwd=2, col=col[1:num]), 
          columns = num
        )
      , type = "b"
      , lwd = 1.5
      , cex = 0.7
      , scales = list(
          x = list(at=seq(from=0, to=36, by=6), cex=1.2), 
          y = list(relation="sliced", cex=1.2)
        )
      , layout = c(num,1)
      , panel = function(x,y,...) {
          panel.xyplot(x,y,...)
          panel.abline(h=0, v=seq(0,36, by=12), color="black", lty=1, lwd=0.5)
        }
    )
    print(b)
  }
  dev.off()

}

####################################################################################
## The input file is the mean.lap.station, key is 1 which is meaningless, value is the 
## dataframe which has mean of 36 mean(abs(error)) and mean of 36 1.96*SD(error) for
## each station. So there are 900 rows in the dataframe.
## Get the dotplot of overmean against sw or tw conditional on station, as well as
## the scatter plot matrix of overmean.
####################################################################################
overallErrorVsStation <- function(type, index, var, target, sub) {

  rst <- rhread(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "overall.absmeanstd")
  )[[1]][[2]]  

  rst <- cbind(rst, parameter[c(rst$group), var]) 

  if(sub) {
    rhload(file.path(rh.root, par$dataset, type, "Rdata", "sample.a1950.RData"))
    rst <- subset(rst, station.id %in% sample.a1950$station.id)
  }
  if(("periodic" %in% parameter$sw) & ("sw" %in% var)){
    rst$sw <- factor(rst$sw, 
      levels=c(sort(as.numeric(unique(rst$sw)[which(unique(rst$sw)!="periodic")])), "periodic")
    )
    rst[,var[which(var != "sw")]] <- as.factor(rst[,var[which(var != "sw")]])
  }else{
    rst[,var[1]] <- as.factor(rst[,var[1]])
    rst[,var[2]] <- as.factor(rst[,var[2]])
  }
  
  if(target == "mean.absmeans") {
    ylab <- "Mean of Absolute Value of Prediction Error"
  } else {
    ylab <- "Mean of Standard Deviation of Error"
  }

  myTable <- tapply(rst$mean.absmeans, rst$group, mean)
  Min <- as.numeric(names(myTable)[which.min(myTable)])
  od.mean <- as.character(arrange(subset(rst, group == Min), get(target))$station.id)
  rst$station.id <- factor(rst$station.id, levels=od.mean)

  num <- with(rst, length(levels(get(var[2]))))
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste(par$dataset, target, "error", var[1], "ps", sep=".")), 
    color = TRUE, 
    paper = "letter"
  )
    b <- dotplot( get(target) ~ station.id | get(var[1])
      , data = rst
      , ylab = list(label = ylab, cex = 1.5)
      , xlab = list(label = "Stations", cex = 1.5)
      , groups = get(var[2])
      , key = list(
          type = "p", pch = 1,
          text = list(label=paste(var[2],"=",with(rst, levels(get(var[2]))), sep="")),  
          points = list(col=col[1:num]), 
          columns = num
        )
      , layout = c(with(rst, length(levels(get(var[1])))), 1)
      , scales = list(
          x = list(draw = FALSE), y = list(relation = "same", cex=1.2)
        )
      , panel = function(x,y,...) {
          #panel.abline(h=seq(0,4,0.5), lwd=0.5, col="lightgray")
          panel.dotplot(x,y,levels.fos=NULL, pch = 1, cex = 0.5,...) 
        }
    )
    print(b)
  dev.off()

  num <- with(rst, length(levels(get(var[1]))))
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste(par$dataset, target, "error", var[2],"ps", sep=".")), 
    color = TRUE, 
    paper = "letter"
  )
    b <- dotplot( get(target) ~ station.id | get(var[2])
      , data = rst
      , ylab = list(label = ylab, cex = 1.5)
      , xlab = list(label = "Stations", cex = 1.5)
      , groups = get(var[1])
      , key = list(
          type = "p", pch = 1,
          text = list(label=paste(var[1],"=",with(rst, levels(get(var[1]))), sep="")),  
          points = list(col=col[1:num]), 
          columns = num
        )
      , layout = c(with(rst, length(levels(get(var[2])))), 1)
      , scales = list(
          x = list(draw = FALSE), y = list(relation = "same", cex=1.2)
        )
      , panel = function(x,y,...) {
          #panel.abline(h=seq(0,4,0.5), lwd=0.5, col="lightgray")
          panel.dotplot(x,y,levels.fos=NULL, pch = 1, cex = 0.5,...) 
        }
    )
    print(b)
  dev.off()

}




