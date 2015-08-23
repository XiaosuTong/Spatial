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

#The output directory for saving plots should have all permission since the plots are 
#written to the output directory by a user related to hadoop.
#In the reduce step, the permission for the plots should be changed to be all +wrx.

#The dataset for Tmax on HDFS is already Ordered the stations by the observation counts

dataset <- "tmax"
index <- "E3"

if(index == "E1"){
    parameter <- expand.grid(
        sw = c(51,73,93), 
        tw = c(617,865,1113), 
        td=2, 
        sd=1
    )
}else if(index == "E2"){
    parameter <- expand.grid(
        sw=c(25, 125, "periodic"), 
        tw=c(617, 865, 1113), 
        td=2, 
        sd=1
    )
    parameter$sw <- as.character(parameter$sw)
}else if(index == "E3"){
    parameter <- expand.grid(
        sw=c(25, 125, "periodic"), 
        tw=c(121, 361, 1141), 
        td=2, 
        sd=1
    )
    parameter$sw <- as.character(parameter$sw)
}else{
    parameter <- expand.grid(
        sw=c(25, 125, "periodic"), 
        tw=c(361, 617, 1141)
    )
    parameter$sw <- as.character(parameter$sw)
}


lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

if(dataset == "tmax") {
    ylab <-  "Maximum Temperature (degrees centigrade)"
}else if(dataset == "tmin") {
    ylab <- "Minimum Temperature (degrees centigrade)"
}else {
    ylab <- "Precipitation (millimeters)"
}


####################################################################################
## The input file is the absmeanstd.lap.station, key is 1 which is meaningless, 
## value is the dataframe which has absmean, mean, and SD of prediction error for 
## a given station, give lap, and given group. So there are 32400 rows in the dataframe.
## Get the scatter plot of mean of abs(error) cross 601 replicates vs lag 
## conditional on group for each station.
####################################################################################
rst <- rhread(
    file.path(
        rh.datadir,
        dataset,
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
## Get the QQ plot of mean and 1.96*std of error cross 601 replicates.
#################################################################################
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("QQ_plot_of_mean_error_", dataset, ".ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- qqmath( ~ means | group,
        data = combined,
        distribution = qunif,
        aspect = 1,
        layout = c(3, 3),
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
        paste("QQ_plot_of_std_error_", dataset, ".ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
    b <- qqmath( ~ std | group,
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



#################################################################################
##For each lap, get a normal quantile plot for each station conditional on group.
##There is an issue that the combined in reduce function is over 128MB which will
##give an error for out of memory. So we only plot in the reduce step, no data is
##saved in this step. The station is ordered by the amount of divation from normal,
##and the order is called from div.station.
#################################################################################
job <- list()
job$map <- expression({
  lapply(seq_along(map.keys), function(r){
	value <- map.values[[r]]
	names(value)[8] <- "lag"
	value$group <- map.keys[[r]][1]
        for(j in 1:dim(parameter)[1]){
            if(unique(value$group) == j){
                value$sw <- parameter[j,1]
                value$tw <- parameter[j,2]
            }
        }
	rhcollect(map.keys[[r]][1], value)
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
        trellis.device(
            device = postscript, 
            file = paste(
                "./tmp", "/QQ_plot_of_error_", 
                dataset, 
                "_group", 
                reduce.key, 
                ".ps", 
                sep = ""
            ), 
            color = TRUE, 
            paper = "legal"
        )
        for(i in as.character(div.stations[[as.numeric(reduce.key)]][[2]]$station.id)){
          b <- qqmath( ~ residual | lag,
            data = subset(combined, station.id==i),
            xlab = list(label = "Unit normal quantile", cex = 1.2),
            ylab = list(label = paste("Station",i, ylab), cex = 1.2),
            type = "p",
            aspect = 1,
            col = "red",
            layout = c(9,4),
            pch=16,
            cex=0.3,
            prepanel = prepanel.qqmathline,
            panel = function(x,...) {
                panel.qqmathline(x, distribution = qnorm)
                panel.qqmath(x,...)
            }
          )
          print(b)
        }
        dev.off()
    }
)
job$setup <- expression(
    reduce = {
        library(lattice)
    }
)
job$parameters <- list(
    parameter = parameter, 
    ylab = ylab,
    div.stations= get(paste(index, "div.stations", sep="."))
)
job$input <- rhfmt(
    file.path(
        rh.datadir, 
        dataset,
        "100stations",
        "sharepredict",
        index,
        "36.lap.station"
    ), 
    type = "sequence"
)
job$output <- rhfmt(
    file.path(rh.output, dataset), 
    type="sequence"
)
job$mapred <- list(
    mapred.reduce.tasks = 10
)
job$jobname <- paste(dataset, "abs error quantile")
job$readback <- FALSE
job$mon.sec <- 10
job$copyFiles <- TRUE
job.mr <- do.call("rhwatch", job)

