###############################################
##Get the kd tree of spatial loess
###############################################

#load the data
library(lattice)
library(maps)
library(plyr)
datadir <- "~/Projects/Spatial/NCAR/RData/"
outputdir <- "~/Projects/Spatial/NCAR/output/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))
#station information is in USpinfo and UStinfo object.
dataset <- "tmax"
if(dataset == "precip"){
	data <- USppt
	info <- USpinfo
}else{
	data <- UStemp
	info <- UStinfo
}
rm(list=grep("US", ls(), value=T))
sub <- subset(data, year==y & month==m)
source("~/Projects/Spatial/NCAR/code/spatial/kdtree.R")
#alpha here is basically the span/5 in loess
kd001 <- kdtree(sub[,c("station.id","lon","lat")], alpha=0.001)
kd003 <- kdtree(sub[,c("station.id","lon","lat")], alpha=0.003)
rhsave(list=("kd001"), file=file.path(rh.datadir, dataset, "Rdata", paste("kd001", ".RData", sep="")))
rhsave(list=("kd003"), file=file.path(rh.datadir, dataset, "Rdata", paste("kd003", ".RData", sep="")))

