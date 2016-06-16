
plotEng.raw <- function(data, station, leaf) {

  data$factor <- factor(
    x = rep(paste("Period", 1:6), c(rep(108,5), 36)),
    levels = paste("Period", c(6:1))
  )
  data$time <- c(rep(0:107, times = 5), 0:35) 

  b <- xyplot( resp ~ time | factor
    , data = data
    , xlab = list(label = "Month", cex = 1.5)
    , ylab = list(label = ylab, cex = 1.5)
    , sub = list(label=paste("Station ", station, "from cell", leaf), cex=1.2)
    , layout = c(1,6)
    , strip = FALSE,
    , xlim = c(0, 107)
    , key=list(
        cex = 1.2,
        text = list(label=c("spatial smoothed value", "temporal fitted value")), 
        lines = list(pch=16, cex=0.7, lwd=1.5, type=c("p","l"), col=col[c(1:2)]),
        columns=2
      )
    , scales = list(
        y = list(tick.number=4, cex=1.2), 
        x = list(at=seq(0, 107, by=12), relation='same', cex=1.2)
      )
    , panel = function(x,y,subscripts,...) {
        panel.abline(v=seq(0,108, by=12), color="lightgrey", lty=3, lwd=0.5)
        panel.xyplot(x, y, type="p", col=col[1], pch=16, cex=0.5, ...)
        panel.xyplot(data[subscripts,]$time, (data[subscripts,]$trend+data[subscripts,]$seasonal+data[subscripts,]$Rspa), type="l", col=col[2], lwd=1, ...)   
        panel.xyplot(data[subscripts,]$time, (data[subscripts,]$trend+data[subscripts,]$seasonal), type="l", col=col[3], lwd=1, lty=3, ...)   
      }
  )
  return(b)

}

tmp <- arrange(rstbystat[[100]][[2]], date)
plotEng.raw(tmp, 1,1)


plotEng.contour <- function(data) {

  b <- levelplot( (seasonal+trend+Rspa) ~ lon * lat
    , data = data
    , region = TRUE
    , col.regions = colorRampPalette(c("blue", "yellow","red"))
    , xlab = list(label="Longitude", cex=1.5)
    , ylab = list(label="Latitude", cex=1.5)
    , xlim = c(-125, -66.5)
    , scale = list(cex=1.2)
    , panel = function(x, y, z, ...) {
       panel.levelplot(x,y,z,...) 
       panel.polygon(us.map$x,us.map$y, border = "black")
     }
  )

}
us.map <- map('state', plot = FALSE, fill = TRUE)

plotEng.contour(rstbymth[[11]][[2]])



#################################
##  mapreduce.task.io.sort.mb  ## 
#################################
load("inputsize/io_sort.RData")
rst <- arrange(rst, rep, io_sort, insize)
rst$x <- rep(c(1,2,3,4), each=24)
trellis.device(
  device = postscript, 
  file = file.path("./io_sort.ps"), 
  color=TRUE, 
  paper="letter"
)
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  )
  sub <- subset(rst, job == i)
  b <- xyplot(log2(elapsed) ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="task.io.sort.mb", cex=1.5)
    , scale = list(x=list(at=c(1,2,3,4), labels=c("128","256","512", "768")), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )
  print(b)
}
dev.off()
trellis.device(
  device = postscript, 
  file = file.path("./iosort_spillratio.ps"), 
  color=TRUE, 
  paper="letter"
)
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  )
  sub <- subset(rst, job == i)
  b <- xyplot(spillRatio ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , aspect = 2
    , main = mainlab
    , ylab = list(label="Spill factor", cex=1.5)
    , xlab = list(label="task.io.sort.mb", cex=1.5)
    , scale = list(x=list(at=c(1,2), labels=c("128", "512")), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )
  print(b)
}
dev.off()

########################################
##  mapreduce.map.sort.spill.percent  ##
########################################
load("inputsize/iospillpercent.RData")
rst$x <- rep(c(1,2,3), each=24)
trellis.device(
  device = postscript, 
  file = file.path("./iospillpercent.ps"), 
  color=TRUE, 
  paper="letter"
)
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  )  
  sub <- subset(rst, job == i)
  b <- xyplot(log2(elapsed) ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="map.sort.spill.percent ", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c("0.2","0.8","1")), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )
  print(b)
}
dev.off()
trellis.device(
  device = postscript, 
  file = file.path("./iospillpercent_spillratio.ps"), 
  color=TRUE, 
  paper="letter"
)
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  )
  sub <- subset(rst, job == i)
  b <- xyplot(spillRatio ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , aspect = 2
    , main = mainlab
    , ylab = list(label="Spill factor", cex=1.5)
    , xlab = list(label="map.sort.spill.percent ", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c("0.2","0.8","1")), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )
  print(b)
}
dev.off()


#####################################################
##   mapreduce.job.reduce.slowstart.completedmaps  ##
#####################################################
load("inputsize/slowstart.RData")
rst$x <- rep(c(1,2,3), each=24)
trellis.device(
  device = postscript, 
  file = file.path("./slowstart.ps"), 
  color=TRUE, 
  paper="letter"
)
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  )
  sub <- subset(rst, job == i)
  b <- xyplot(log2(elapsed) ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="job.reduce.slowstart", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c(0.2, 0.5, 0.8)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )
  print(b)
}
dev.off()


############################
##  mapreduce.job.reduce  ##
############################
load("inputsize/redtask.RData")
trellis.device(
  device = postscript, 
  file = file.path("./redtask.ps"), 
  color=TRUE, 
  paper="letter"
)
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  )
  sub <- subset(rst, job == i)
  b <- xyplot(log2(elapsed) ~ redtask | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="job.reduce", cex=1.5)
    , scale = list(x=list(at=c(90, 179, 269, 358)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )
  print(b)
}
dev.off()


######################################
##   mapreduce.task.io.sort.factor  ##
######################################
load("inputsize/iosortfactor.RData")
rst$x <- rep(c(1,2,3,4), each=24)
trellis.device(
  device = postscript, 
  file = paste0("./", "iosortfactor", ".ps"), 
  color=TRUE, 
  paper="letter"
)  
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  )    
  sub <- subset(rst, job == i)
  b <- xyplot(log2(elapsed) ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="task.io.sort.factor", cex=1.5)
    , scale = list(x=list(at=c(1,2,3,4), labels=c(2,10,20,50)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)
}
dev.off()


################################################
##   mapreduce.reduce.shuffle.parallelcopies  ##
################################################
load("inputsize/parallcopy.RData")
rst$x <- rep(c(1,2,3), each=24)
trellis.device(
  device = postscript, 
  file = paste0("./", "parallcopy", ".ps"), 
  color=TRUE, 
  paper="letter"
)  
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  ) 
  sub <- subset(rst, job == i)
  b <- xyplot(log2(elapsed) ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="reduce.parallelcopies", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c(5,20,50)),cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)
}
dev.off()


##############################################
##   mapreduce.reduce.input.buffer.percent  ##
##############################################
load("inputsize/redinpercent.RData")
trellis.device(
  device = postscript, 
  file = paste0("./", "redinpercent", ".ps"), 
  color=TRUE, 
  paper="letter"
)  
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  ) 
  sub <- subset(rst, job == i)
  b <- xyplot(log2(elapsed) ~ redinpercent | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="reduce.input.buffer.percent", cex=1.5)
    , scale = list(x=list(at=c(0,0.5,1)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)
}
dev.off()


#####################################################
##  mapreduce.reduce.shuffle.input.buffer.percent  ##
#####################################################
load("inputsize/shuffleinput.RData")
rst$x <- rep(c(1,2,3,4), each=24)
trellis.device(
  device = postscript, 
  file = paste0("./", "shuffleinput", ".ps"), 
  color=TRUE, 
  paper="letter"
)  
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  ) 
  sub <- subset(rst, job == i)
  b <- xyplot(log2(elapsed) ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="reduce.shuffle.input.buffer.percent", cex=1.5)
    , scale = list(x=list(at=c(1,2,3,4), labels=c(0.4,0.7,0.9,1)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)
}
dev.off()


##############################################
##  mapreduce.reduce.shuffle.merge.percent  ##
##############################################
load("inputsize/shufflepercent.RData")
rst$x <- rep(c(1,2,3), each=24)
trellis.device(
  device = postscript, 
  file = paste0("./", "shufflepercent", ".ps"), 
  color=TRUE, 
  paper="letter"
)  
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  )
  sub <- subset(rst, job == i)
  b <- xyplot(log2(elapsed) ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="reduce.shuffle.merge.percent ", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c(0.3,0.6,1)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)
}
dev.off()


####################################
##  two factors for map spilling  ##
####################################
load("inputsize/mapmemory.RData")
rst$x <- rep(c(3,1,2),each=24)
trellis.device(
  device = postscript, 
  file = paste0("./", "mapmemory", ".ps"), 
  color=TRUE, 
  paper="letter"
)  
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  )
  sub <- subset(rst, job == i)
  b <- xyplot(log2(elapsed) ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="map memory setting", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c("160/1","200/0.8","800/0.2")), cex=1)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)
}
dev.off()
trellis.device(
  device = postscript, 
  file = file.path("./mapmeory_spillratio.ps"), 
  color=TRUE, 
  paper="letter"
)
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  )
  sub <- subset(rst, job == i)
  b <- xyplot(spillRatio ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , aspect = 2
    , ylim = c(0.9,2.1)
    , main = mainlab
    , ylab = list(label="Spill factor", cex=1.5)
    , xlab = list(label="map memory setting", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c("160/1","200/0.8","800/0.2")), cex=1)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )
  print(b)
}
dev.off()


###################################
##  two factors for reduce spill ##
###################################
load("inputsize/redmemory.RData")
rst$x <- rep(c(2,1,3),each=24)
trellis.device(
  device = postscript, 
  file = paste0("./", "redmemory", ".ps"), 
  color=TRUE, 
  paper="letter"
)  
for(i in c("readin","swaptoloc","swaptotime","swaptoloc2")) {
  mainlab <- paste(
    "Step", 
    which(c("readin","spofit", "swaptoloc", "stlfit", "swaptotime", "sprfit", "swaptoloc2") == i)
  )
  sub <- subset(rst, job == i)
  b <- xyplot(log2(elapsed) ~ x | as.factor(insize)
    , data = sub
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="reduce memory setting", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c("0.45/1","0.9/0.5","1/0.45")), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)
}
dev.off()


load("inputsize/jump.RData")
rst$x <- rep(c(1,2,3,4), each=6)
trellis.device(
  device = postscript, 
  file = paste0("./", "jump", ".ps"), 
  color=TRUE, 
  paper="letter"
)  
  b <- xyplot(log2(elapsed) ~ x | as.factor(insize)
    , data = rst
    , layout = c(6,1)
    , main = mainlab
    , ylab = list(label="log base 2 elapsed time", cex=1.5)
    , xlab = list(label="jump", cex=1.5)
    , scale = list(x=list(at=c(1,2,3,4), labels=c(5,10,50,100)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.xyplot(x,y,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)
dev.off()


####################################################################
##                                                                ##
##  pilot experiment results, only input size=1365 is considered  ##
##                                                                ##
####################################################################
rst <- subset(rst, insize==1365)
rst$job <- c("Step I","Step III","Step V","Step VII")

trellis.device(
  device = postscript, 
  file = file.path("~/io_sort.ps"), 
  color=TRUE, 
  paper="letter"
)

  b <- xyplot(log2(elapsed) ~ x | as.factor(job)
    , data = rst
    , layout = c(4,1)
    , ylab = list(label="Log base 2 elapsed time", cex=1.5)
    , xlab = list(label="ISM", cex=1.5)
    , scale = list(x=list(at=c(1,2,3,4), labels=c("128","256","512", "768")), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.abline(h=c(8.5,9,9.5), v=c(1,2,3,4), col="lightgray", lwd=1)
      panel.xyplot(x,y,cex=1.1,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )
  print(b)

dev.off()

rst$x <- rep(c(1,2,3), each=4)
trellis.device(
  device = postscript, 
  file = file.path("~/iospillpercent.ps"), 
  color=TRUE, 
  paper="letter"
)

  b <- xyplot(log2(elapsed) ~ x | as.factor(job)
    , data = rst
    , layout = c(4,1)
    , ylab = list(label="Log base 2 elapsed time", cex=1.5)
    , xlab = list(label="SSP", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c("0.2","0.8","1")), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.abline(h=c(8.5,9,9.5), v=c(1,2,3), col="lightgray", lwd=1)
      panel.xyplot(x,y, cex=1.1,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )
  print(b)

dev.off()

trellis.device(
  device = postscript, 
  file = paste0("~/", "parallcopy", ".ps"), 
  color=TRUE, 
  paper="letter"
)  

  b <- xyplot(log2(elapsed) ~ x | as.factor(job)
    , data = rst
    , layout = c(4,1)
    , ylab = list(label="Log base 2 elapsed time", cex=1.5)
    , xlab = list(label="SPC", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c(5,20,50)),cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.abline(h=c(8.5,9,9.5), v=c(1,2,3), col="lightgray", lwd=1)
      panel.xyplot(x,y,cex=1.1,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)

dev.off()


trellis.device(
  device = postscript, 
  file = file.path("~/slowstart.ps"), 
  color=TRUE, 
  paper="letter"
)

  b <- xyplot(log2(elapsed) ~ x | as.factor(job)
    , data = rst
    , layout = c(4,1)
    , ylab = list(label="Log base 2 elapsed time", cex=1.5)
    , xlab = list(label="RSC", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c(0.2, 0.5, 0.8)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.abline(h=c(8.5,9,9.5), v=c(1,2,3), col="lightgray", lwd=1)
      panel.xyplot(x,y,cex=1.1,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )
  print(b)

dev.off()

trellis.device(
  device = postscript, 
  file = file.path("~/redtask.ps"), 
  color=TRUE, 
  paper="letter"
)

  b <- xyplot(log2(elapsed) ~ x | as.factor(job)
    , data = rst
    , layout = c(4,1)
    , ylab = list(label="Log base 2 elapsed time", cex=1.5)
    , xlab = list(label="RTSK", cex=1.5)
    , scale = list(x=list(at=c(1,2,3,4), labels=c(90, 179, 269, 358)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.abline(h=c(8.5,9,9.5), v=c(1,2,3,4), col="lightgray", lwd=1)
      panel.xyplot(x,y,cex=1.1,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )
  print(b)

dev.off()

trellis.device(
  device = postscript, 
  file = paste0("~/", "redinpercent", ".ps"), 
  color=TRUE, 
  paper="letter"
)  

  b <- xyplot(log2(elapsed) ~ redinpercent | as.factor(job)
    , data = rst
    , layout = c(4,1)
    , ylab = list(label="Log base 2 elapsed time", cex=1.5)
    , xlab = list(label="IBP", cex=1.5)
    , scale = list(x=list(at=c(0,0.5,1)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.abline(h=c(8.5,9,9.5), v=c(0,0.5,1), col="lightgray", lwd=1)
      panel.xyplot(x,y,cex=1.1,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)

dev.off()

trellis.device(
  device = postscript, 
  file = paste0("~/", "shuffleinput", ".ps"), 
  color=TRUE, 
  paper="letter"
)  

  b <- xyplot(log2(elapsed) ~ x | as.factor(job)
    , data = rst
    , layout = c(4,1)
    , ylab = list(label="Log base 2 elapsed time", cex=1.5)
    , xlab = list(label="SIBP", cex=1.5)
    , scale = list(x=list(at=c(1,2,3,4), labels=c(0.4,0.7,0.9,1)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.abline(h=c(8.5,9,9.5), v=c(1,2,3,4), col="lightgray", lwd=1)
      panel.xyplot(x,y,cex=1.1,...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)

dev.off()

trellis.device(
  device = postscript, 
  file = paste0("~/", "shufflepercent", ".ps"), 
  color=TRUE, 
  paper="letter"
)  

  b <- xyplot(log2(elapsed) ~ x | as.factor(job)
    , data = rst
    , layout = c(4,1)
    , ylab = list(label="Log base 2 elapsed time", cex=1.5)
    , xlab = list(label="SMP", cex=1.5)
    , scale = list(x=list(at=c(1,2,3), labels=c(0.3,0.6,0.9)), cex=1.2)
    , panel = function(x,y, subscripts, ...) {
      panel.abline(h=c(8.5,9,9.5), v=c(1,2,3), col="lightgray", lwd=1)
      panel.xyplot(x,y,cex=1.1, ...)
      panel.average(x,y,fun = mean, horizontal = FALSE, col="red",...)
    }
  )  
  print(b)

dev.off()