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
index <- "Emultiple2"

#Load the station.id for the 100 stations for Tmax
load(file.path(local.datadir, paste(dataset,"div.stations.RData", sep="")))
load(file.path(local.datadir, "stations.RData"))
if(index == "Emultiple"){
    parameter <- expand.grid(sw=c(25, 125, "periodic"), tw=1141, td=1, sd=1, fcw=c(121, 361, 1141), fcd=1, ssw=825, ssd=2)
    parameter$sw <- as.character(parameter$sw)
}else{
    parameter <- expand.grid(sw="periodic", tw=1141, td=1, sd=1, fcw=c(495, 811, 1141), fcd=1, ssw=c(361,511,825), ssd=2)
    parameter$sw <- as.character(parameter$sw)
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
rst <- rhread(file.path(rh.datadir,dataset,"100stations","sharepredict",index,"absmeanstd.lap.station"))
combined <- rst[[1]][[2]]
library(lattice)
combined <- combined[with(combined, order(station.id, fcw, ssw, lap)), ]
#if(any(grepl("periodic", unique(combined$sw)))){
#combined$sw <- factor(combined$sw, levels=c(sort(as.numeric(unique(combined$sw)[which(unique(combined$sw)!="periodic")])), "periodic"))
#combined$ssw <- as.factor(combined$ssw)
#}else{
#combined$sw <- as.factor(combined$sw)
#combined$ssw <- as.factor(combined$ssw)
#}
combined$ssw <- as.factor(combined$ssw)
combined$fcw <- as.factor(combined$fcw)
trellis.device(postscript, file = file.path(local.output, paste("QQ_plot_of_mean_abserror_", dataset, ".ps", sep="")), color=TRUE, paper="legal")
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
trellis.device(postscript, file = file.path(local.output, paste("scatter_plot_of_mean_abserror_", dataset, "_vs_lag_fcw.ps", sep="")), color=TRUE, paper="legal")
	b <- xyplot( absmeans ~ lap | fcw*as.factor(station.id),
             data = combined,
             xlab = list(label = "Lag", cex = 1.2),
             ylab = list(label = ylab, cex = 1.2),
             groups = as.factor(ssw),
             key=list(type="l", text=list(label=levels(as.factor(combined$ssw))),  lines=list(lwd=1.5, col=col[1:length(levels(combined$ssw))]), columns=length(levels(combined$ssw))),
             type = "b",
             scales = list(x=list(at=seq(from=0, to=36, by=6)), y=list(relation="sliced")),
             layout = c(length(levels(combined$fcw)),1),
             pch=1,
	     	 cex =0.3,
             panel = function(x,y,...) {
                  panel.xyplot(x,y,...)
                  panel.abline(h=0, v=seq(0,36, by=12), color="black", lty=1, lwd=0.5)
             }
	)
	print(b)
dev.off()
trellis.device(postscript, file = file.path(local.output, paste("scatter_plot_of_mean_abserror_", dataset,"_vs_lag_ssw.ps", sep="")), color=TRUE, paper="legal")
	b <- xyplot( absmeans ~ lap | ssw*as.factor(station.id),
             data = combined,
             xlab = list(label = "Lag", cex = 1.2),
             ylab = list(label = ylab, cex = 1.2),
             groups = as.factor(fcw),
             key=list(type="l", text=list(label=levels(as.factor(combined$fcw))),  lines=list(lwd=1.5, col=col[1:length(levels(combined$fcw))]), columns=length(levels(combined$fcw))),
             type = "b",
             scales = list(x=list(at=seq(from=0, to=36, by=6)), y=list(relation="sliced")),
             layout = c(length(levels(combined$ssw)),1),
             pch=1,
	     cex = 0.3,
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
trellis.device(postscript, file = file.path(local.output, paste("QQ_plot_of_mean_error_", dataset, ".ps", sep="")), color=TRUE, paper="legal")
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
trellis.device(postscript, file = file.path(local.output, paste("QQ_plot_of_std_error_", dataset, ".ps", sep="")), color=TRUE, paper="legal")
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
trellis.device(postscript, file = file.path(local.output, paste("scatter_plot_of_mean_error_", dataset, "_vs_lag_fcw.ps", sep="")), color=TRUE, paper="legal")
	b <- xyplot( means ~ lap | fcw*as.factor(station.id),
             data = combined,
             xlab = list(label = "Lag", cex = 1.2),
             ylab = list(label = ylab, cex = 1.2),
             groups = ssw,
             key=list(type="l", text=list(label=levels(combined$ssw)),  lines=list(lwd=1.5, col=col[1:length(levels(combined$ssw))]), columns=length(levels(combined$ssw))),
             type = "b",
	     	 ylim = c(min(combined$means), max(combined$means)),
             scales = list(x=list(at=seq(from=0, to=36, by=6)), y=list(relation="same")),
             layout = c(length(levels(combined$fcw)),1),
             pch=1,
	     	 cex=0.3,
             panel = function(x,y,...) {
                  panel.xyplot(x,y,...)
                  panel.abline(h=0, v=seq(0,36, by=12), color="black", lty=1, lwd=0.5)
             }
        )
        print(b)
dev.off()
trellis.device(postscript, file = file.path(local.output, paste("scatter_plot_of_mean_error_", dataset,"_vs_lag_ssw.ps", sep="")), color=TRUE, paper="legal")
	b <- xyplot( means ~ lap | ssw*as.factor(station.id),
             data = combined,
             xlab = list(label = "Lag", cex = 1.2),
             ylab = list(label = ylab, cex = 1.2),
             groups = fcw,
             key=list(type="l", text=list(label=levels(combined$fcw)),  lines=list(lwd=1.5, col=col[1:length(levels(combined$fcw))]), columns=length(levels(combined$fcw))),
             type = "b",
             ylim = c(min(combined$means), max(combined$means)),
             scales = list(x=list(at=seq(from=0, to=36, by=6)), y=list(relation="same")),
             layout = c(length(levels(combined$ssw)),1),
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
trellis.device(postscript, file = file.path(local.output, paste("scatter_plot_of_std_error_", dataset, "_vs_lag_fcw.ps", sep="")), color=TRUE, paper="legal")
	b <- xyplot( std ~ lap | fcw*as.factor(station.id),
             data = combined,
             xlab = list(label = "Lag", cex = 1.2),
             ylab = list(label = ylab, cex = 1.2),
             groups = ssw,
             key=list(type="l", text=list(label=levels(combined$ssw)),  lines=list(lwd=1.5, col=col[1:length(levels(combined$ssw))]), columns=length(levels(combined$ssw))),
             type = "b",
             scales = list(x=list(at=seq(from=0, to=36, by=6)), y=list(relation="sliced")),
             layout = c(length(levels(combined$fcw)),1),
             pch=1,
	     cex =0.3,
             panel = function(x,y,...) {
                  panel.xyplot(x,y,...)
                  panel.abline(h=0, v=seq(0,36, by=12), color="black", lty=1, lwd=0.5)
             }
        )
        print(b)
dev.off()
trellis.device(postscript, file = file.path(local.output, paste("scatter_plot_of_std_error_", dataset,"_vs_lag_ssw.ps", sep="")), color=TRUE, paper="legal")
	b <- xyplot( std ~ lap | ssw*as.factor(station.id),
             data = combined,
             xlab = list(label = "Lag", cex = 1.2),
             ylab = list(label = ylab, cex = 1.2),
             groups = fcw,
             key=list(type="l", text=list(label=levels(combined$fcw)),  lines=list(lwd=1.5, col=col[1:length(levels(combined$fcw))]), columns=length(levels(combined$fcw))),
             type = "b",
             scales = list(x=list(at=seq(from=0, to=36, by=6)), y=list(relation="sliced")),
             layout = c(length(levels(combined$ssw)),1),
             pch=1,
	     cex=0.3,
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
rst <- rhread(file.path(rh.datadir, dataset, "100stations","sharepredict",index,"over.absmeanstd.lap.station"))
result <- as.data.frame(do.call(rbind, lapply(rst, "[[", 1)), stringsAsFactors=FALSE)
result <- cbind(result, as.data.frame(do.call(rbind, lapply(rst, "[[", 2))))
names(result) <- c("station.id","fcw","ssw", "overmean", "overstd")

result$ssw <- as.factor(as.numeric(result$ssw))
result$fcw <- as.factor(as.numeric(result$fcw))
#if(any(grepl("periodic", unique(result$sw)))){
#        result$sw <- factor(result$sw, levels=c(sort(as.numeric(unique(result$sw)[which(unique(result$sw)!="periodic")])), "periodic"))
#        result$ssw <- as.factor(result$ssw)
#}else{
#        result$sw <- as.factor(result$sw)
#        result$ssw <- as.factor(result$ssw)
#}
od.mean <- names(sort(tapply(result$overmean, result$station.id, mean)))
result$station.id <- factor(result$station.id, levels=od.mean)
trellis.device(postscript, file = file.path(local.output, paste("dotplot_of_overmean_error_", dataset,"_fcw.ps", sep="")), color=TRUE, paper="legal")
	b <- dotplot( overmean ~ station.id | fcw,
             data = result,
             ylab = list(label = "Mean of absolute prediction error", cex = 1.2),
             xlab = list(label = "Stations", cex = 1.2),
             groups = ssw,
             key=list(type="p", text=list(label=levels(result$ssw)),  points=list(col=col[1:length(levels(result$ssw))]), columns=length(levels(result$ssw)), pch=1),
             type = "p",
             layout = c(1,1),
             scales = list(x=list(draw=FALSE), y=list(relation="free")), 
#             panel = function(x,y,...) {
#                  panel.xyplot(x,y,...)
#                  panel.abline(h=0, v=seq(0,36, by=12), color="lightgrey", lty=3, lwd=0.5)
 #            }
        )
	print(b)
dev.off()

trellis.device(postscript, file = file.path(local.output, paste("dotplot_of_overmean_error_", dataset,"_ssw.ps", sep="")), color=TRUE, paper="legal")
        b <- dotplot( overmean ~ station.id | ssw,
             data = result,
             ylab = list(label = "Mean of absolute prediction error", cex = 1.2),
             xlab = list(label = "Stations", cex = 1.2),
             groups = fcw,
             key=list(type="l", text=list(label=levels(result$fcw)),  points=list(pch=1, col=col[1:length(levels(result$fcw))]), columns=length(levels(result$fcw))),
             type = "p",
             layout = c(1,1),
             scales = list(x=list(draw=FALSE), y=list(relation="free")),
#                  panel.xyplot(x,y,...)
#                  panel.abline(h=0, v=seq(0,36, by=12), color="lightgrey", lty=3, lwd=0.5)
 #            }
        )
        print(b)
dev.off()

od.std <- names(sort(tapply(result$overstd, result$station.id, mean)))
result$station.id <- factor(result$station.id, levels=od.std)
trellis.device(postscript, file = file.path(local.output, paste("dotplot_of_overstd_error_", dataset,"_fcw.ps", sep="")), color=TRUE, paper="legal")
        b <- dotplot( overstd ~ station.id | fcw,
             data = result,
             ylab = list(label = "Mean of 1.96*SD of prediction error", cex = 1.2),
             xlab = list(label = "Stations", cex = 1.2),
             groups = ssw,
             key=list(type="p", text=list(label=levels(result$ssw)),  points=list(col=col[1:length(levels(result$ssw))]), columns=length(levels(result$ssw)), pch=1),
             type = "p",
             layout = c(1,1),
             scales = list(x=list(draw=FALSE), y=list(relation="free")),
#             panel = function(x,y,...) {
#                  panel.xyplot(x,y,...)
#                  panel.abline(h=0, v=seq(0,36, by=12), color="lightgrey", lty=3, lwd=0.5)
 #            }
        )
        print(b)
dev.off()

trellis.device(postscript, file = file.path(local.output, paste("dotplot_of_overstd_error_", dataset,"_ssw.ps", sep="")), color=TRUE, paper="legal")
        b <- dotplot( overstd ~ station.id | ssw,
             data = result,
             ylab = list(label = "Mean of 1.96*SD prediction error", cex = 1.2),
             xlab = list(label = "Stations", cex = 1.2),
             groups = fcw,
             key=list(type="l", text=list(label=levels(result$fcw)),  points=list(pch=1, col=col[1:length(levels(result$fcw))]), columns=length(levels(result$fcw))),
             type = "p",
             layout = c(1,1),
             scales = list(x=list(draw=FALSE), y=list(relation="free")),
#             panel = function(x,y,...) {
#                  panel.xyplot(x,y,...)
#                  panel.abline(h=0, v=seq(0,36, by=12), color="lightgrey", lty=3, lwd=0.5)
 #            }
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
                    value$td <- parameter[j,3]
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
        library(lattice)
        trellis.device(postscript, file = paste(local.output, "/QQ_plot_of_error_", dataset, "_group", reduce.key, ".ps", sep = ""), color=TRUE, paper="legal")
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
        system(paste("chmod 777 ", local.output, "/QQ_plot_of_error_", dataset, "_group", reduce.key, ".ps", sep = ""))
   }
)


job$parameters <- list(parameter = parameter, div.stations= get(paste(index, "div.stations", sep=".")))
job$input <- rhfmt(file.path(rh.datadir, dataset,"100stations","sharepredict",index,"36.lap.station"), type="sequence")
job$output <- rhfmt(file.path(rh.output, dataset), type="sequence")
job$mapred <- list(mapred.reduce.tasks = 72)
job$jobname <- paste(dataset, "abs error quantile")
job$readback <- FALSE
job$mon.sec <- 10
job.mr <- do.call("rhwatch", job)


#####################################################################################
##The input file is 36.lap.quantile, the key is quantile(0.05, 0.25, 0.5, 0.75, 0.95)
##For each quantile, get a scatter plot of error vs lag conditional on group
#####################################################################################
job <- list()
job$map <- expression({
    lapply(seq_along(map.values), function(r){
        rhcollect(map.keys[[r]], map.values[[r]])
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
        library(lattice)
	combined <- combined[with(combined, order(station.id, lap)), ]
        trellis.device(postscript, file = paste(local.output, "/scatter_plot_of_error_", dataset, "_tw_quantile", reduce.key, ".ps", sep = ""), color=TRUE, paper="legal")
        for(i in unique(combined$station.id)){
          b <- xyplot( abs(quantiles) ~ lap | as.factor(tw),
             data = subset(combined, station.id==i),
             xlab = list(label = "Lag", cex = 1.2),
             ylab = list(label = paste("Station",i, dataset), cex = 1.2),
	     main = list(label = paste("Absolute errors of prediction:", reduce.key), cex=1.3),
	     groups = as.factor(sw),
             key=list(type="l", text=list(label=levels(as.factor(subset(combined, station.id==i)$sw))),  lines=list(lwd=1.5, col=col[1:3]), columns=3),
             type = "b",
             #aspect = 1,
	     scales = list(x=list(at=seq(from=0, to=36, by=6))), 
             layout = c(3,1),
             pch=1,
	     panel = function(x,y,...) {
                  panel.xyplot(x,y,...)
             }

          )
          print(b)
        }
        dev.off()
        trellis.device(postscript, file = paste(local.output, "/scatter_plot_of_error_", dataset, "_sw_quantile", reduce.key, ".ps", sep = ""), color=TRUE, paper="legal")
        for(i in unique(combined$station.id)){
          b <- xyplot( abs(quantiles) ~ lap | as.factor(sw),
             data = subset(combined, station.id==i),
             xlab = list(label = "Lag", cex = 1.2),
             ylab = list(label = paste("Station",i, dataset), cex = 1.2),
             main = list(label = paste("Absolute errors of prediction:", reduce.key), cex=1.3),
             groups = as.factor(tw),
             key=list(type="l", text=list(label=levels(as.factor(subset(combined, station.id==i)$tw))),  lines=list(lwd=1.5, col=col[1:3]), columns=3),
             type = "b",
             #aspect = 1,
             scales = list(x=list(at=seq(from=0, to=36, by=6))),
             layout = c(3,1),
             pch=1,
             panel = function(x,y,...) {
                  panel.xyplot(x,y,...)
             }

          )
          print(b)
        }
        dev.off()
        system(paste("chmod 777 ", local.output, "/scatter_plot_of_error_",dataset,"_tw_quantile", reduce.key, ".ps", sep = ""))
        system(paste("chmod 777 ", local.output, "/scatter_plot_of_error_",dataset,"_sw_quantile", reduce.key, ".ps", sep = ""))
   }
)

job$parameters <- list(col=col)
job$input <- rhfmt(file.path(rh.datadir, dataset,"100stations","sharepredict",index,"36.lap.quantile"), type="sequence")
job$output <- rhfmt(file.path(rh.output, dataset), type="sequence")
job$mapred <- list(mapred.reduce.tasks = 72)
job$jobname <- paste(dataset, "abs error quantile")
job$readback <- FALSE
job$mon.sec <- 10
#job.mr <- do.call("rhwatch", job)

