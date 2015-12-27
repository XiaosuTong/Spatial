##################################################################
##  Compare two specific stlplus model based on the remainders  ##
##  mod1 and mod2 are two vectors, each of which contains first ##
##  the model file name on HDFS, second the model index on the  ##
##  webpage. If summary is TRUE, specific for type="a1950", a   ##
##  summary statistic will be calculated for each of 7738       ##
##  stations, and finally one data.frame will be combined to    ##
##################################################################
fitcompare <- function(type = "a1950", mod1, mod2, summary=FALSE) {
  if (type == "a1950") {
    ## If summary is TRUE, we will calculate a summary statistic for each station
    ## and combine all station to one data.frame
    if (summary) {  

      job <- list()
      job$map <- expression({
        lapply(seq_along(map.keys), function(r) {
          value <- data.frame(
            station.id = map.keys[[r]], 
            meanerror=mean(map.values[[r]]$remainder, na.rm=TRUE),
            stringsAsFactors = FALSE
          )    
          rhcollect(1, value)
        })
      })
      job$reduce <- expression(
        pre = {    
          combine <- data.frame()
        },
        reduce = {    
          combine <- rbind(combine, do.call("rbind", reduce.values))
        },
        post = {    
          rhcollect(reduce.key, combine)
        }
      )
      files <- c(
        file.path(rh.root, par$dataset, type, "STL", mod1[1]),
        file.path(rh.root, par$dataset, type, "STL", mod2[1])
      )
      job$input <- rhfmt(files, type = "sequence")
      job$output <- rhfmt(
        file.path(rh.root, par$dataset, type, "STLcompare", paste(mod1[1], "VS", mod2[1], ".summary", sep="")),
        type = "sequence"
      )
      job$mapred <- list(
        mapreduce.job.reduces = 1,  #cdh5
        mapred.reduce.tasks = 1
      )
      job$jobname <- file.path(rh.root, par$dataset, type, "STLcompare", paste(mod1[1], mod2[1], sep="VS"))
      job$combiner <- TRUE
      job$readback <- FALSE
      job$mon.sec <- 10
      job.mr <- do.call("rhwatch", job)  

    } else { ## else if summary is FALSE, only keep 128 stations and all residual of each 128 stations to one data.frame  

      job <- list()
      job$map <- expression({
        lapply(seq_along(map.keys), function(r) {  
        	if(map.keys[[r]] %in% sample.a1950$station.id) {
            value <- map.values[[r]][, "remainder", drop=FALSE]
            value$station.id <- map.keys[[r]]
            file <- Sys.getenv("mapred.input.file")
            span <- substr(tail(strsplit(file, "/")[[1]],3)[2], 3, 7)
            rhcollect(1, span)
        	}
        })
      })
      job$reduce <- expression(
        pre = {  
          combine <- data.frame()
        },
        reduce = {  
          combine <- rbind(combine, do.call("rbind", reduce.values))
        },
        post = {  
          rhcollect(reduce.key, combine)
        }
      )
      job$shared <- file.path(rh.datadir, par$dataset, type, "Rdata", "sample.a1950.RData")
      job$setup <- expression(
        map = {
        	library(plyr)
        	load("sample.a1950.RData")
        }
      )
      job$parameters <- list(
        mod1 = mod1[1],
        mod2 = mod2[1]
      )
      files <- c(
        file.path(rh.root, par$dataset, type, "STL", mod1[1]),
        file.path(rh.root, par$dataset, type, "STL", mod2[1])
      )
      job$input <- rhfmt(files, type = "sequence")
      job$output <- rhfmt(
        file.path(rh.root, par$dataset, type, "STLcompare", paste(mod1[1], mod2[1], sep="VS")),
        type = "sequence"
      )
      job$mapred <- list(
        mapreduce.job.reduces = 1,  #cdh5
        mapred.reduce.tasks = 1
      )
      job$jobname <- file.path(rh.root, par$dataset, type, "STLcompare", paste(mod1[1], mod2[1], sep="VS"))
      job$combiner <- TRUE
      job$readback <- FALSE
      job$mon.sec <- 10
      job.mr <- do.call("rhwatch", job)  

    }
  } else {
    ## for 100stations, no need for MapReduce job, just read back the fit result from 
    ## HDFS
    rst1 <- rhread(file.path(rh.root, par$dataset, type, "STL", mod1[1]))[[1]][[2]]
    rst2 <- rhread(file.path(rh.root, par$dataset, type, "STL", mod2[1]))[[1]][[2]]
    
    rst <- merge(
      x = subset(rst1, select = c(station.id, remainder, date)), 
      y = subset(rst2, select = c(station.id, remainder, date)),
      by = c("station.id","date")
    )

  }
  
  fitcompare.plot(rst, type=type, mod1=mod1[2], mod2=mod2[2])

}

############################################################
##  Plotting function called in the fitcompare function.  ##
##  mod1 and mod2 here are the model index on the webpage ##
############################################################
## the input data is a data.frame containing two columns, station.id, remainder, and date
fitcompare.plot <- function(rst, type, mod1, mod2) {
  
  ## first plot of comparison is QQ plot of remainders
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("QQ.compare", "ps", sep=".")), 
    color = TRUE, 
    paper = "letter"
  )
  b <- xyplot( remainder.x ~ remainder.y | station.id
    , data = rst
    , col = "red"
    , pch = 16
    , cex = 0.3
    , aspect = 1
    , layout = c(1,1)
    , scale = list(cex = 1.2, relation="free")
    , xlab = list(label=paste(mod1, "remainder"), cex=1.5)
    , ylab = list(label=paste(mod2, "remainder"), cex=1.5)
    , panel = function(x,y,...) {
        panel.abline(a=c(0,1), col="black", lwd=0.5, lty=1)
        panel.xyplot(sort(x), sort(y),...)
      }
  )
  print(b)
  dev.off()

  ## Second plot of comparison is the scatter plot of remainders
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("scatter.compare", "ps", sep=".")), 
    color = TRUE, 
    paper = "letter"
  )
  b <- xyplot( remainder.x ~ remainder.y | station.id
    , data = rst
    , col = "red"
    , pch = 16
    , cex = 0.3
    , aspect = 1
    , layout = c(1,1)
    , scale = list(cex = 1.2, relation="free")
    , xlab = list(label=paste(mod1, "remainder"), cex=1.5)
    , ylab = list(label=paste(mod2, "remainder"), cex=1.5)
    , panel = function(x,y,...) {
        panel.abline(a=c(0,1), col="black", lwd=0.5, lty=1)
        panel.xyplot(x, y,...)
      }
  )
  print(b)
  dev.off()

  ## Third plot of comparison is the dotplot of mean of remainders
  remainderMean <- ddply(
    .data = rst,
    .vari = "station.id",
    .fun = summarise,
    mean.x = mean(abs(remainder.x), na.rm=TRUE),
    mean.y = mean(abs(remainder.y), na.rm=TRUE),
    std.x = sd(remainder.x, na.rm=TRUE),
    std.y = sd(remainder.y, na.rm=TRUE)
  )
  trellis.device(
    device = postscript, 
    file = file.path(local.root, "output", paste("dotplot.compare", "ps", sep=".")), 
    color = TRUE, 
    paper = "letter"
  )
  for(r in c("mean","std")) {
    target <- paste(r, "y", sep=".")
    if (r == "mean") {
      ylab <- "Mean of Absolute Value of Remainder"
    } else {
      ylab <- "Standard Deviation of Remainder"
    }
    remainderMean <- arrange(remainderMean, get(target))
    remainderMean$station.id <- factor(remainderMean$station.id, levels=remainderMean$station.id)
    b <- dotplot( get(target) ~ station.id
      , data = remainderMean
      , scale = list(cex = 1.2, x=list(draw=FALSE))
      , xlab = list(label = "Station", cex=1.5)
      , ylab = list(label = ylab, cex=1.5)
      , key=list(
          text = list(label=c(mod1, mod2)),
          lines = list(pch=16, type="p", col=col[2:1]), 
          columns = 2
        )
      , panel = function(x,y,...) {
          panel.dotplot(x, y, col = col[1], pch=16, ...)
          panel.dotplot(x, remainderMean[, paste(r, "x", sep=".")], col=col[2], ...)
        }
    )
    print(b)
  }
  dev.off()


}