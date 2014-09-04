#################################################################
##100 stations, compare multiple stl2 fitting
#################################################################

#load the data
library(lattice)
library(plyr)
library(stl2)
outputdir <- "~/Projects/Spatial/NCAR/output/"
datadir <- "~/Projects/Spatial/NCAR/RData/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))
load(paste(datadir,"stations.RData", sep=""))
dataset <- "tmax"
fc <- FALSE

if(fc){
    parameter <- list(
        run1=list(
            sw="periodic", 
            tw=1141, 
            td=1, 
            sd=1, 
            inner=10, 
            outer=0, 
            flag=FALSE
        ),
        run2=list(
            sw="periodic", 
            tw=1141, 
            sd=1, 
            td=1, 
            fcw=1141, 
            fcd=1, 
            ssw=825, 
            ssd=2, 
            inner=10, 
            outer=0, 
            flag=TRUE
        )
    )
}else{
    parameter <- expand.grid(sw="periodic", tw=c(241,1141), td=1, sd=1, inner=10, outer=0)
    parameter$sw <- as.character(parameter$sw)
}
if(dataset == "tmax"){
ylab <- "Maximum Temperature (degrees centigrade)"
}else if(dataset == "tmin"){
        ylab <- "Minimum Temperature (degrees centigrade)"
}else {
        ylab <- "Precipitation (millimeters)"
}
if(dataset %in% c("tmax", "tmin")){
        data <- UStemp
        datainfo <- UStinfo
}else{
        data <- USppt
        ddatainfo <- USpinfo
}
rm(list=grep("US", ls(), value=TRUE))

#find the 100 stations with largest observation number for max temperature
stations <- get(paste("stations", dataset, sep="."))
tmp <- data[with(data, station.id %in% stations),]
if(any(with(tmp, is.na(get(dataset))) == TRUE)) stop("The first 100 stations have NA")

month <- tmp$month
levels(month) <- c(1:12)
month <- as.numeric(factor(month, levels=c(1:12)))
date <- paste(tmp$year, month, "01", sep="-")
tmp$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
column <- c("station.id", "year","month","date", dataset)
tmp <- tmp[with(tmp, order(station.id, date)), column]
names(tmp)[dim(tmp)[2]] <- "response"

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

if(!fc){
	rst <- alply(parameter, 1, function(x){
		ddply(tmp, "station.id", function(r){stl2(r$response, r$date, n.p=12, s.window=x$sw, s.degree=x$sd, t.window=x$tw, t.degree=x$td, inner=x$inner, outer=x$outer)$data})
})
}else{
	rst <- llply(parameter, function(x){
	if(x$flag){
	  ddply(tmp, "station.id", function(r){do.call("cbind",stl2(r$response, r$date, n.p=12, s.window=x$sw, s.degree=x$sd, t.window=x$tw, t.degree=x$td, fc.window=c(x$fcw, x$ssw), fc.degree=c(x$fcd, x$ssd), inner=x$inner, outer=x$outer)[c("data","fc")])
		})
	}else{
      ddply(tmp, "station.id", function(r){stl2(r$response, r$date, n.p=12, s.window=x$sw, s.degree=x$sd, t.window=x$tw, t.degree=x$td, inner = x$inner, outer = x$outer)$data})
	}
	})
	rst <- llply(rst, function(x){
	names(x)[grep("fc", names(x))] <- c("fc.trend", "fc.second", "remainder")
	x
})
}
####################################################################
##
####################################################################
remainders <- as.data.frame(do.call(cbind, lapply(rst, function(x){x$remainder})))
remainders$station.id <- rep(stations, each=1236)
names(remainders)[1:2] <- c("model1", "model2")

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmax_models_comparison",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
        b <- xyplot( model1 ~ model2,
             data = subset(remainders, station.id==i),
             xlab = list(label = "model 1.2.4 remainder", cex = 1.2),
             ylab = list(label = "model 1.2.5 remainder", cex = 1.2),
             main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = "p",
             pch = 16,
	    	 col = "red",
             cex = 0.3,
	     	 aspect = 1,
#            scales = list(y=list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             panel = function(...) {
                  panel.abline(a=0, b=1, color="black", lty=1)
                  panel.xyplot(...)
             }
        )
        print(b)
   }
dev.off()

trellis.device(postscript, file = paste(outputdir, "scatterplot_of_tmax_models_comparison_diff",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
        b <- xyplot( (model1-model2) ~ model1,
             data = subset(remainders, station.id==i),
             xlab = list(label = "model 1.2.4 remainder", cex = 1.2),
             ylab = list(label = "model 1.2.4 - model 1.2.6 remainder", cex = 1.2),
             main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = "p",
             pch = 16,
             col = "red",
             cex = 0.3,
#             aspect = 1,
#            scales = list(y = list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             panel = function(...) {
                  panel.abline(h=0,v=0, color="black", lty=1)
                  panel.xyplot(...)

             }
        )
        print(b)
   }
dev.off()



tmp2 <- remainders[,-1]
tmp1 <- remainders[,-2]
names(tmp1)[1] <- names(tmp2)[1] <- "remainder"
tmp1$group <- rep("model1", 123600)
tmp2$group <- rep("model2", 123600)
tmp.remainders <- rbind(tmp1, tmp2)
trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_tmax_models_comparison",".ps", sep = ""), color=TRUE, paper="legal")
     for(i in stations){
        b <- qq( group ~ remainder,
             data = subset(tmp.remainders, station.id==i),
             xlab = list(label = "model 1.2.4 remainder", cex = 1.2),
             ylab = list(label = "model 1.2.6 remainder", cex = 1.2),
             main = list(label = paste("Station ", i, sep=""), cex=1.5),
             type = "p",
             pch = 16,
             col = "red",
             cex = 0.3,
             aspect = 1,
#            scales = list(y=list(relation = 'free', cex=1.5), x=list(relation= 'free',format = "%b %Y", tick.number=10), cex=1.2),
             panel = function(...) {
                  panel.abline(a=0, b=1, color="black", lty=1)
                  panel.xyplot(...)

             }
        )
        print(b)
   }
dev.off()

ACF1 <- ddply(
    .data=remainders,
    .variables="station.id",
    .fun= summarise,
    correlation = c(acf(model1, plot=FALSE)$acf),
    lag = c(acf(model1, plot=FALSE)$lag)
)
ACF1 <- subset(ACF1, lag!=0)
ACF2 <- ddply(
    .data=remainders,
    .variables="station.id",
    .fun= summarise,
    correlation = c(acf(model2, plot=FALSE)$acf),
    lag = c(acf(model2, plot=FALSE)$lag)
)
ACF2 <- subset(ACF2, lag!=0)
ACF2$lag <- ACF2$lag + 0.1
ACF <- rbind(ACF1, ACF2)
ACF$group <- rep(c("model1","model2"), each=3000)
trellis.device(postscript, file = paste(outputdir, "acf_of_", dataset, "_with_stl2_remainder_for_100_stations",".ps", sep = ""), color=TRUE, paper="legal")
   for(i in stations){
        b <- xyplot( correlation ~ lag,
             data = subset(ACF, station.id==i),
             xlab = list(label = "Lag", cex = 1.2),
             ylab = list(label = paste("Station", i, "ACF"), cex = 1.2),
			 key=list(type="l", text=list(label=as.character(parameter$tw)), lines=list(col=col[1:2], lwd=1.5), columns=2),
             type = "h",
			 groups = group,
             panel = function(x,y,...) {
                  panel.abline(h=0)
                  panel.xyplot(x,y,...)
             }
        )
        print(b)
   }
dev.off()


summary <- ddply(
    .data=remainders, 
	.variables = "station.id",
	.fun = summarise,
	mean = c(mean(model1), mean(model2)),
	std = c(sd(model1), sd(model2)),
	group = c("model1", "model2")
)
m <- ddply(
    .data=summary,
	.variables="group",
	.fun=function(r){
		arrange(r, mean)
    }
)
s <- ddply(
    .data=summary,
    .variables="group",
    .fun=function(r){
        arrange(r, std)
    }
)
trellis.device(postscript, file = paste(outputdir, "dotplot_of_tmax_models_comparison",".ps", sep = ""), color=TRUE, paper="legal")
		od.mean <- tail(m$station.id, 100)
		summary$station.id <- factor(summary$station.id, levels=od.mean)
        b <- dotplot( mean ~ station.id,
             data = summary,
             ylab = list(label = "Mean of Remainders", cex = 1.2),
             xlab = list(label = "Stations", cex = 1.2),
			 key=list(type="p", text=list(label=as.character(parameter$tw)),  points=list(col=col[1:2], pch=1), columns=2),
			 groups = group,
             type = "p",
             scales = list(x=list(draw=FALSE), y=list(relation="free")),
             panel = function(x,y,...){
                  panel.abline(h=0, color="black")
                  panel.dotplot(x,y,...)
             }
        )
        od.std <- head(s$station.id, 100)
        summary$station.id <- factor(summary$station.id, levels=od.std)
		a <- dotplot( std ~ station.id,
             data = summary,
             ylab = list(label = "Std of Remainders", cex = 1.2),
             xlab = list(label = "Stations", cex = 1.2),
             key=list(type="p", text=list(label=as.character(parameter$tw)),  points=list(col=col[1:2], pch=1), columns=2),
             groups = group,
             type = "p",
             scales = list(x=list(draw=FALSE), y=list(relation="free")),
             panel = function(x,y,...){
                  panel.abline(h=0, color="black")
                  panel.dotplot(x,y,...)
             }
        )
		print(b)
        print(a)
dev.off()

