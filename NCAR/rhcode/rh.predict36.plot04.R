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
source("~/Rhipe/rhinitial.R")
dataset <- "tmax"
index <- "E5"
par <- list()
par$machine <- "gacrux"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

#Load the station.id for the 100 stations for Tmax
#load(file.path(local.datadir, paste(dataset,"div.stations.RData", sep="")))
if(index == "E5"){
    parameter <- list(
        run1=list(
            sw="periodic", 
            tw=1855, 
            sd=1, 
            td=1, 
            fcw=1855, 
            fcd=1, 
            ssw=241, 
            ssd=1, 
            multiple=TRUE
        ),
        run2=list(
            sw="periodic", 
            tw=241, 
            sd=1, 
            td=1, 
            multiple=FALSE
        ),
        run3=list(
            sw="periodic", 
            tw=1855, 
            sd=1,
            td=1, 
            fcw=1855, 
            fcd=1, 
            ssw=121, 
            ssd=2, 
            multiple=TRUE
        )
    )
}
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col
if(dataset == "tmax") {
ylab <- "Maximum Temperature"
}else if(dataset == "tmin") {
ylab <- "Minimum Temperature"
}else {ylab <- "Precipitation"}


####################################################################################
## The input file is the absmeanstd.lap.station, key is 1 which is meaningless, 
## value is the dataframe which has absmean, mean, and SD of prediction error for 
## a given station, give lap, and given group. So there are 32400 rows in the dataframe.
## Get the scatter plot of mean of abs(error) cross 601 replicates vs lag 
## conditional on group for each station.
####################################################################################
rst <- rhread(
    file.path(
        rh.datadir,dataset,
        "100stations",
        "sharepredict",
        index,
        "absmeanstd.lap.station"
    )
)
combined <- rst[[1]][[2]]
library(lattice)
combined$group <- as.numeric(combined$group)
combined <- combined[with(combined, order(station.id, group, lap)), ]
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("QQ_plot_of_mean_abserror_", dataset, ".ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- qqmath( ~ absmeans,	
		data = combined,
		distribution = qunif,
		aspect = 1,
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
        paste("scatter_plot_of_mean_abserror_", dataset, "_vs_lag.ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- xyplot( absmeans ~ lap | as.factor(station.id),
        data = combined,
        xlab = list(label = "Lag", cex = 1.2),
        ylab = list(label = ylab, cex = 1.2),
        groups = group,
        key = list(
            type = "l", 
            text = list(
                label = c("two frequency components","one frequency component")
            ),  
            lines = list(
                lwd = 1.5, 
                col = col[1:length(unique(combined$group))]
            ), 
            columns = length(unique(combined$group))
        ),
        type = "b",
        scales = list(
            x = list(at=seq(from=0, to=36, by=6)), 
            y = list(relation="sliced")
        ),
        layout = c(4,1),
        pch=1,
        cex =0.3,
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
	b <- qqmath( ~ means,
        data = combined,
        distribution = qunif,
        aspect = 1,
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
    b <- qqmath( ~ std,
        data = combined,
        distribution = qunif,
        aspect = 1,
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
        paste("scatter_plot_of_mean_error_", dataset, "_vs_lag.ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- xyplot( means ~ lap | as.factor(station.id),
        data = combined,
        xlab = list(label = "Lag", cex = 1.2),
        ylab = list(label = ylab, cex = 1.2),
        groups = group,
        key = list(
            type = "l", 
            text = list(
                label = c("two frequency components","one frequency component")
            ),
            lines = list(
                lwd = 1.5, 
                col = col[1:length(unique(combined$group))]
            ), 
            columns = length(unique(combined$group))
        ),
        type = "b",
        ylim = c(min(combined$means), max(combined$means)),
        scales = list(
            x = list(at=seq(from=0, to=36, by=6)), 
            y = list(relation="same")
        ),
        layout = c(4,1),
        pch=1,
        cex=0.3,
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
        paste("scatter_plot_of_std_error_", dataset, "_vs_lag.ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- xyplot( std ~ lap | as.factor(station.id),
        data = combined,
        xlab = list(label = "Lag", cex = 1.2),
        ylab = list(label = ylab, cex = 1.2),
        groups = group,
        key = list(
            type="l", 
            text= list(
                label=c("two frequency components","one frequency component")
            ),  
            lines = list(
                lwd = 1.5, 
                col = col[1:length(unique(combined$group))]
            ), 
            columns = length(unique(combined$group))
        ),
        type = "b",
        scales = list(
            x = list(at=seq(from=0, to=36, by=6)), 
            y = list(relation="sliced")
        ),
        layout = c(4,1),
        pch=1,
        cex =0.3,
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
result <- as.data.frame(do.call(rbind, lapply(rst, "[[", 1)), stringsAsFactors=FALSE)
result <- cbind(result, as.data.frame(do.call(rbind, lapply(rst, "[[", 2))))
names(result) <- c("station.id","group","overmean", "overstd")
result$group <- as.numeric(result$group)
result <- result[with(result, order(station.id, group)), ]
od.mean <- names(sort(tapply(result$overmean, result$station.id, mean)))
result$station.id <- factor(result$station.id, levels=od.mean)
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("dotplot_of_overmean_error_", dataset,".ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
	b <- dotplot( overmean ~ station.id,
        data = result,
        ylab = list(label = "Mean of absolute prediction error", cex = 1.2),
        xlab = list(label = "Stations", cex = 1.2),
        groups = group,
        key=list(
            type = "p", 
            text = list(
                label=c("two frequency components","one frequency component")
            ),  
            points = list(col=col[1:length(unique(result$group))]), 
            columns = length(unique(result$group)), 
            pch = 1
        ),
        type = "p",
        layout = c(1,1),
        scales = list(
            x = list(draw=FALSE), 
            y = list(relation="free")
        ), 
    )
	print(b)
dev.off()

od.std <- names(sort(tapply(result$overstd, result$station.id, mean)))
result$station.id <- factor(result$station.id, levels=od.std)
trellis.device(
    device = postscript, 
    file = file.path(
        local.output, 
        paste("dotplot_of_overstd_error", dataset,".ps", sep="")
    ), 
    color = TRUE, 
    paper = "legal"
)
    b <- dotplot( overstd ~ station.id,
        data = result,
        ylab = list(label = "Mean of 1.96*SD of prediction error", cex = 1.2),
        xlab = list(label = "Stations", cex = 1.2),
        groups = group,
		key = list(
            type = "p", 
            text = list(
                label=c("two frequency components","one frequency component")
            ),  
            points = list(
                col = col[1:length(unique(result$group))]
            ), 
            columns = length(unique(result$group)), 
            pch = 1
        ),
        type = "p",
        layout = c(1,1),
        scales = list(
            x=list(draw=FALSE), 
            y=list(relation="free")
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
div.stations <- rhread(
    file.path(
        rh.datadir, 
        dataset, 
        "100stations",
        "sharepredict",
        index, 
        "group.orderstations"
    )
)
job <- list()
job$map <- expression({
    lapply(seq_along(map.keys), function(r){
	   value <- map.values[[r]]
	   names(value)[grep("lap", names(value))] <- "lag"
	   value$group <- map.keys[[r]][1]
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
            pch = 16,
            cex = 0.3,
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
    div.stations = div.stations
)
job$input <- rhfmt(
    file.path(rh.datadir, dataset,"100stations","sharepredict",index,"36.lap.station"), 
    type = "sequence"
)
job$output <- rhfmt(
    file.path(rh.output, dataset), 
    type = "sequence"
)
job$mapred <- list(
    mapred.reduce.tasks = 10
)
job$jobname <- paste(dataset, "abs error quantile")
job$readback <- FALSE
job$mon.sec <- 10
job$copyFiles <- TRUE
job.mr <- do.call("rhwatch", job)
