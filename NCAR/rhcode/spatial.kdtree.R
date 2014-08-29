###############################################
##Get the kd tree of spatial loess
###############################################

#load the data
library(lattice)
library(maps)
library(plyr)

load(file.path(local.datadir, "USmonthlyMet.RData"))
#station information is in USpinfo and UStinfo object.
save(list = grep("info", ls(), value=T), file = "~/Projects/Spatial/NCAR/RData/info.RData")
dataset <- "tmax"
if(dataset == "precip"){
	data <- USppt
	info <- USpinfo
}else{
	data <- UStemp
	info <- UStinfo
}

#using my kd-tree function to get the kd-tree vertices
source("~/Projects/Spatial/NCAR/code/spatial/kdtree.R")
#alpha here is basically the span/5 in loess
kd001 <- kdtree(sub[,c("station.id","lon","lat")], alpha=0.001)
kd003 <- kdtree(sub[,c("station.id","lon","lat")], alpha=0.003)
rhsave(list=("kd001"), file=file.path(rh.datadir, dataset, "Rdata", paste("kd001", ".RData", sep="")))
rhsave(list=("kd003"), file=file.path(rh.datadir, dataset, "Rdata", paste("kd003", ".RData", sep="")))
