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

#######################################################################
##Order the stations based on the deviation from the normal distribution
##For each station, given lap, given group, the distribution of 601 
##prediction error is compared with normal distribution. The sum of abs
##deviation over all lap for one station under given group is calculated.
#######################################################################
####################################################################################
## The input file is the absmeanstd.lap.station, key is 1 which is meaningless, 
## value is the dataframe which has absmean, mean, and SD of prediction error for 
## a given station, give lap, and given group. So there are 32400 rows in the dataframe.
## Get the scatter plot of mean of abs(error) cross 601 replicates vs lag 
## conditional on group for each station.
####################################################################################
QQDivFromNormal <- function(index, type) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r){
      if(map.keys[[r]][1] %in% sample.a1950$station.id) {
        v <- map.values[[r]]
        v$group <- map.keys[[r]][1]
        yy <- quantile(v$residual, c(0.25, 0.75))
        xx <- qnorm(c(0.25, 0.75))
        r <- diff(yy)/diff(xx)
        x <- qnorm(ppoints(length(v$residual)))
        y <- r*x + yy[1] - xx[1]*r
        div <- sum(abs(sort(v$residual) - y))
        rhcollect(map.keys[[r]][1:2], div)
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
  job$setup <- expression(
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
  job$jobname <- file.path(rh.root, par$dataset, type, "STLtuning", index, "divorderstations")
  job$readback <- FALSE
  job$mon.sec <- 10
  job.mr <- do.call("rhwatch", job)

}


#################################################################################
## Get the QQ plot of mean and 1.96*std of error cross 601 replicates.
#################################################################################
plotEngine <- function(data, dataset, key) {

  trellis.device(
    device = postscript, 
    file = file.path(".", "tmp", paste("QQ.error", dataset, "group", key, "ps", sep= ".")), 
    color = TRUE, 
    paper = "legal"
  )
  for(i in unique(data$station.id)){
    b <- qqmath( ~ residual | as.numeric(lag)
      , data = subset(data, station.id==i)
      , xlab = list(label = "Unit normal quantile", cex = 1.5)
      , ylab = list(label = ylab, cex = 1.5)
      , sub = paste("Station", i)
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
      plotEngine(data=combined, dataset=dataset, key=reduce.key)
    }
  )
  job$shared <- c(file.path(rh.root, par$dataset, type, "Rdata", "sample.a1950.RData"))
  job$setup <- expression(
    map = {load("sample.a1950.RData")},
    reduce = {library(lattice)}
  )
  job$parameters <- list(
    param = param, 
    ylab = ylab,
    dataset = par$dataset,
    plotEngine = plotEngine
    #div.stations= get(paste(index, "div.stations", sep="."))
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


#################################################################################
##For each lap, get a normal quantile plot for each station conditional on group.
##There is an issue that the combined in reduce function is over 128MB which will
##give an error for out of memory. So we only plot in the reduce step, no data is
##saved in this step. The station is ordered by the amount of divation from normal,
##and the order is called from div.station.
#################################################################################

rst <- rhread(
    file.path(
        rh.datadir,
        par$dataset,
        "100stations",
        "sharepredict",
        index,
        "absmeanstd.lap.station"
    )
)
combined <- rst[[1]][[2]]
library(lattice)
combined <- combined[with(combined, order(station.id, sw, tw, lap)), ]
if(index != "E1"){
	combined$sw <- factor(
        combined$sw, 
        levels=c(sort(as.numeric(unique(combined$sw)[which(unique(combined$sw)!="periodic")])), "periodic")
    )
	combined$tw <- as.factor(combined$tw)
}else{
	combined$sw <- as.factor(combined$sw)
	combined$tw <- as.factor(combined$tw)
}
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("QQ_plot_of_mean_abserror_", dataset, ".ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- qqmath( ~ absmeans | group,	
		data = combined,
		distribution = qunif,
		aspect = 1,
        layout = c(3 ,3),
		xlab = list(label = "f-value", cex = 1.2),
		ylab = list(label = ylab, cex = 1.2),
		pch = 1,
		cex = 0.3,
	)
	print(b)
dev.off()
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("scatter_plot_of_mean_abserror_", dataset, "_vs_lag_sw.ps", sep="")
    ), 
    color=TRUE, 
    paper="legal"
)
	b <- xyplot( absmeans ~ lap | sw*as.factor(station.id),
        data = combined,
        xlab = list(label = "Lag", cex = 1.2),
        ylab = list(label = ylab, cex = 1.2),
        groups = as.factor(tw),
        key=list(
            type = "l", 
            text = list(label=levels(as.factor(combined$tw))),  
            lines = list(lwd=1.5, col=col[1:3]), 
            columns=3
        ),
        type = "b",
        scales = list(
            x=list(at=seq(from=0, to=36, by=6)), 
            y=list(relation="sliced")
        ),
        layout = c(3,1),
        pch = 1,
        panel = function(x,y,...) {
            panel.xyplot(x,y,...)
            panel.abline(h=0, v=seq(0,36, by=12), color="black", lty=1, lwd=0.5)
        }
	)
	print(b)
dev.off()
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("scatter_plot_of_mean_abserror_", dataset,"_vs_lag_tw.ps", sep="")
    ), 
    color = TRUE, 
    paper="legal"
)
	b <- xyplot( absmeans ~ lap | tw*as.factor(station.id),
        data = combined,
        xlab = list(label = "Lag", cex = 1.2),
        ylab = list(label = ylab, cex = 1.2),
        groups = as.factor(sw),
        key=list(
            type = "l", 
            text = list(label=levels(as.factor(combined$sw))),  
            lines = list(lwd=1.5, col=col[1:3]), 
            columns = 3
        ),
        type = "b",
        scales = list(
            x = list(at = seq(from = 0, to = 36, by = 6)), 
            y = list(relation = "sliced")
        ),
        layout = c(3,1),
        pch = 1,
        panel = function(x,y,...) {
            panel.xyplot(x,y,...)
            panel.abline(h=0, v=seq(0,36, by=12), color="black", lty=1, lwd=0.5)
        }
    )
	print(b)
dev.off()


#################################################################################
## Get the scatter plot of mean of error cross 601 replicates vs lag 
## conditional on group for each station.
#################################################################################
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("scatter_plot_of_mean_error_", dataset, "_vs_lag_sw.ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- xyplot( means ~ lap | sw*as.factor(station.id),
        data = combined,
        xlab = list(label = "Lag", cex = 1.2),
        ylab = list(label = ylab, cex = 1.2),
        groups = as.factor(tw),
        key=list(
            type = "l", 
            text = list(label=levels(as.factor(combined$tw))),  
            lines = list(lwd=1.5, col=col[1:3]), 
            columns = 3
        ),
        type = "b",
        ylim = c(min(combined$means), max(combined$means)),
        scales = list(
            x = list(at = seq(from = 0, to = 36, by = 6)), 
            y = list(relation = "same")
        ),
        layout = c(3,1),
        pch = 1,
        panel = function(x,y,...) {
            panel.xyplot(x,y,...)
            panel.abline(h=0, v=seq(0,36, by=12), color="black", lty=1, lwd=0.5)
        }
    )
    print(b)
dev.off()
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("scatter_plot_of_mean_error_", dataset,"_vs_lag_tw.ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- xyplot( means ~ lap | tw*as.factor(station.id),
        data = combined,
        xlab = list(label = "Lag", cex = 1.2),
        ylab = list(label = ylab, cex = 1.2),
        groups = as.factor(sw),
        key=list(
            type = "l", 
            text = list(label=levels(as.factor(combined$sw))),  
            lines = list(lwd=1.5, col=col[1:3]), 
            columns = 3
        ),
        type = "b",
        ylim = c(min(combined$means), max(combined$means)),
        scales = list(
            x = list(at = seq(from = 0, to = 36, by = 6)),
            y = list(relation = "same")
        ),
        layout = c(3,1),
        pch = 1,
        panel = function(x,y,...) {
            panel.xyplot(x,y,...)
            panel.abline(h=0, v=seq(0,36, by=12), color="black", lty=1, lwd=0.5)
        }
    )
    print(b)
dev.off()


#################################################################################
## Get the scatter plot of 1.96*std of error cross 601 replicates vs lag 
## conditional on group for each station.
#################################################################################
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("scatter_plot_of_std_error_", dataset, "_vs_lag_sw.ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- xyplot( std ~ lap | sw*as.factor(station.id),
        data = combined,
        xlab = list(label = "Lag", cex = 1.2),
        ylab = list(label = ylab, cex = 1.2),
        groups = as.factor(tw),
        key=list(
            type = "l", 
            text = list(label=levels(as.factor(combined$tw))),  
            lines = list(lwd=1.5, col=col[1:3]), 
            columns = 3
        ),
        type = "b",
        scales = list(
            x = list(at = seq(from = 0, to = 36, by = 6)), 
            y = list(relation = "sliced")
        ),
        layout = c(3,1),
        pch=1,
        panel = function(x,y,...) {
            panel.xyplot(x,y,...)
            panel.abline(h=0, v=seq(0,36, by=12), color="black", lty=1, lwd=0.5)
        }
    )
    print(b)
dev.off()

trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("scatter_plot_of_std_error_", dataset,"_vs_lag_tw.ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- xyplot( std ~ lap | tw*as.factor(station.id),
        data = combined,
        xlab = list(label = "Lag", cex = 1.2),
        ylab = list(label = ylab, cex = 1.2),
        groups = as.factor(sw),
        key=list(
            type = "l", 
            text = list(label=levels(as.factor(combined$sw))),  
            lines = list(lwd=1.5, col=col[1:3]), 
            columns = 3
        ),
        type = "b",
        scales = list(
            x = list(at = seq(from = 0, to = 36, by = 6)), 
            y = list(relation = "sliced")
        ),
        layout = c(3,1),
        pch = 1,
        panel = function(x,y,...) {
            panel.xyplot(x,y,...)
            panel.abline(h=0, v=seq(0,36, by=12), color="black", lty=1, lwd=0.5)
        }
    )
    print(b)
dev.off()


####################################################################################
## The input file is the mean.lap.station, key is 1 which is meaningless, value is the 
## dataframe which has mean of 36 mean(abs(error)) and mean of 36 1.96*SD(error) for
## each station. So there are 900 rows in the dataframe.
## Get the dotplot of overmean against sw or tw conditional on station, as well as
## the scatter plot matrix of overmean.
####################################################################################
rst <- rhread(
    file.path(
        rh.datadir, 
        dataset, 
        "100stations",
        "sharepredict",
        index,
        "over.absmeanstd.lap.station"
    )
)
result <- as.data.frame(
    do.call(rbind, lapply(rst, "[[", 1)), 
    stringsAsFactors=FALSE
)
result <- cbind(
    result, 
    as.data.frame(do.call(rbind, lapply(rst, "[[", 2)))
)
names(result) <- c("station.id","sw","tw", "overmean", "overstd")
#result$sw <- as.numeric(result$sw)
result$tw <- as.numeric(result$tw)
if(index != "E1"){
    result$sw <- factor(
        result$sw, 
        levels = c(
            sort(as.numeric(unique(result$sw)[which(unique(result$sw)!="periodic")])), 
            "periodic"
        )
    )
    result$tw <- as.factor(result$tw)
}else{
    result$sw <- as.factor(result$sw)
    result$tw <- as.factor(result$tw)
}
od.mean <- names(sort(tapply(result$overmean, result$station.id, mean)))
result$station.id <- factor(result$station.id, levels=od.mean)
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("dotplot_of_overmean_error_", dataset,"_tw.ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- dotplot( overmean ~ station.id | tw,
        data = result,
        ylab = list(label = "Mean of absolute prediction error", cex = 1.2),
        xlab = list(label = "Stations", cex = 1.2),
        groups = sw,
        key = list(
            type = "p", 
            text = list(label=levels(as.factor(result$sw))),  
            points = list(col=col[1:3]), 
            columns = 3, 
            pch = 1
        ),
        type = "p",
        layout = c(1,1),
        scales = list(
            x = list(draw = FALSE), 
            y = list(relation = "free")
        )
    )
    print(b)
dev.off()

trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("dotplot_of_overmean_error_", dataset,"_sw.ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
    b <- dotplot( overmean ~ station.id | sw,
        data = result,
        ylab = list(label = "Mean of absolute prediction error", cex = 1.2),
        xlab = list(label = "Stations", cex = 1.2),
        groups = as.factor(tw),
        key = list(
            type = "l", 
            text = list(label=levels(as.factor(result$tw))),  
            points = list(pch=1, col=col[1:3]), 
            columns = 3
        ),
        type = "p",
        layout = c(1,1),
        scales = list(
            x = list(draw = FALSE), 
            y = list(relation = "free")
        )
    )
    print(b)
dev.off()

od.std <- names(sort(tapply(result$overstd, result$station.id, mean)))
result$station.id <- factor(result$station.id, levels=od.std)
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("dotplot_of_overstd_error_", dataset,"_tw.ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
    b <- dotplot( overstd ~ station.id | tw,
        data = result,
        ylab = list(label = "Mean of 1.96*SD of prediction error", cex = 1.2),
        xlab = list(label = "Stations", cex = 1.2),
        groups = as.factor(sw),
        key=list(
            type = "p", 
            text = list(label=levels(as.factor(result$sw))),  
            points = list(col=col[1:3]), 
            columns = 3, 
            pch = 1
        ),
        type = "p",
        layout = c(1,1),
        scales = list(
            x = list(draw = FALSE), 
            y = list(relation = "free")
        ),
    )
    print(b)
dev.off()

trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("dotplot_of_overstd_error_", dataset,"_sw.ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
    b <- dotplot( overstd ~ station.id | as.factor(sw),
        data = result,
        ylab = list(label = "Mean of 1.96*SD prediction error", cex = 1.2),
        xlab = list(label = "Stations", cex = 1.2),
        groups = as.factor(tw),
        key = list(
            type = "l", 
            text = list(label = levels(as.factor(result$tw))),  
            points = list(pch = 1, col = col[1:3]), 
            columns = 3 
        ),
        type = "p",
        layout = c(1,1),
        scales = list(
            x = list(draw = FALSE), 
            y = list(relation = "free")
        ),
    )
    print(b)
dev.off()



