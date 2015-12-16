#Set up the directory and load the data.
library(lattice)
library(plyr)
library(hexbin)
library(maps)
library(magrittr)

options(stringsAsFactors=FALSE)
options(java.parameters = "-Xmx1024m")

#source("~/Rhipe/ross.initial.R")
source("~/Rhipe/rhinitial.R")

par <- list()
par$dataset <- "tmax"
par$Machine <- "adhara"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")


if(par$dataset == "precip") {
  Nstations <- 11918
} else {
  Nstations <- 8125
}

#lattice.theme <- trellis.par.get()
#col <- lattice.theme$superpose.symbol$col
col <- c("#0080ff", "#ff00ff", "darkgreen", "#ff0000", "orange", "#00ff00", "brown") 
if(par$dataset == "tmax") {
  ylab <-  "Maximum Temperature (degrees centigrade)"
}else if(par$dataset == "tmin") {
  ylab <- "Minimum Temperature (degrees centigrade)"
}else {
  ylab <- "Precipitation (millimeters)"
}

#########################################
##    Copy raw data files to HDFS      ##
#########################################
## first copy the station information on to HDFS
USpinfo <- read.table(
  file = file.path(local.root, "Raw", "NCAR_pinfill", "METAinfo"), 
  sep = "", skip = 1,
  col.names = c("station.id", "elev","lon", "lat"),
  stringsAsFactors = FALSE
)
UStinfo <- read.table(
  file = file.path(local.root, "Raw", "NCAR_tinfill", "METAinfo"), 
  sep = "", skip = 1,
  col.names = c("station.id", "elev","lon", "lat"),
  stringsAsFactors = FALSE
)
rhsave(list=c("UStinfo","USpinfo"), file = file.path(rh.root, "stationinfo", "USinfo.RData"))

##then copy the raw data file to HDFS
source("~/Projects/Spatial/NCAR/rhcode/All/RawtoHDFS.R")



#################################################################
##                                                             ##
##               Dataset of the 100 Stations                   ##
##                                                             ##
#################################################################

############################################################################
##        Raw and stl fit for the 100 stations with full obs              ##
############################################################################     
## Rhipe job creating data.frame of 100 stations  ##
s100()
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "aggregated"))[[1]][[2]]
## create plots for raw observations  ##
#source(file.path(local.root, "code", "raw.visual", "100stations.plot.R"))
scatterRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
monthRaw(data=rst, outputdir=file.path(local.root, "output"), target = "tmax", size = "letter", test = F)
monthQQ(data=rst, outputdir=file.path(local.root, "output"), target = "tmax", size = "letter", test = F)
for (i in c("lon","lat","elev")) {
  monthSpatial(data=rst, outputdir=file.path(local.root, "output"), target = "tmax", vari=i, size = "letter")
}

## Rhipe job running stl2 on each station of 100 stations ##
#source(file.path(local.root, "rhcode", "s100", "100stations.plot.R"))
s100.STLfit(sw=77, sd=1, tw=495, td=2, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t495td2_s77sd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)

s100.STLfit(sw=103, sd=1, tw=617, td=2, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t617td2_s77sd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)

s100.STLfit(sw="periodic", sd=1, tw=1141, td=2, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1141td2_speriodicsd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)

s100.STLfit(sw="periodic", sd=1, tw=1141, td=1, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1141td1_speriodicsd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)

s100.STLfit(sw="periodic", sd=1, tw=241, td=1, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t241td1_speriodicsd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)

s100.STLfit(sw="periodic", sd=1, tw=1855, td=1, fcw=121, fcd=2)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1855td1_speriodicsd1_f121fd2"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = TRUE)

s100.STLfit(sw="periodic", sd=1, tw=1855, td=1, fcw=241, fcd=1)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1855td1_speriodicsd1_f241fd1"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = TRUE)

s100.STLfit(sw=103, sd=2, tw=1141, td=1, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1141td1_s103sd2_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)


s100.STLfit(sw=77, sd=2, tw=1855, td=1, fcw=241, fcd=1)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1855td1_s77sd2_f241fd1"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = TRUE)


##################################################################
##    STL tunning fit for the 100 stations with full obs        ##
################################################################## 
FileInput <- file.path(rh.root, par$dataset, "100stations", "bystation")
stationSplit(reduce=200, input=FileInput)
type <- "100stations"
index <- "E1"
parameter <- expand.grid(
  sw = c(51, 73, 93), tw = c(617, 865, 1113), td = 2, sd = 1, fc.flag = FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type, parameter, k, index, valid=600))
}

index <- "E2"
parameter <- expand.grid(
  sw = c(25, 125, "periodic"), tw = c(617, 865, 1113), td = 2, 
  sd = 1, fc.flag = FALSE,stringsAsFactors=FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type, parameter, k, index))
}

index <- "E3"
parameter <- expand.grid(
  sw = c(25, 125, "periodic"), tw = c(121, 361, 1141), 
  td = 2, sd = 1, fc.flag=FALSE, stringsAsFactors=FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type, parameter, k, index))
}
  
index <- "E4"
parameter <- expand.grid(
  sw = c("periodic"), tw = c(121, 241, 361, 751, 1141), 
  td = c(1,2), sd = 1, fc.flag = FALSE, stringsAsFactors = FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type, parameter, k, index))
}

index <- "E5"
parameter <- do.call("rbind", list(
  data.frame(sw="periodic", tw=1855, sd=1, td=1, fcw=1855, fcd=1, scw=241, scd=1, fc.flag=TRUE),
  data.frame(sw="periodic", tw=241, sd=1, td=1, fcw=NA, fcd=NA, scw=NA, scd=NA, fc.flag=FALSE),
  data.frame(sw="periodic", tw=1855, sd=1, td=1, fcw=1855, fcd=1, scw=121, scd=2, fc.flag=TRUE)
))
for(k in 1:nrow(parameter)) {
  try(predict36(type, parameter, k, index))
}

index <- "E6"
parameter <- expand.grid(
  sw = "periodic", tw = 1855, td = 1, sd = 1,fcw = 1855, fcd = 1, 
  scw = c(241, 361, 751), scd = c(1, 2), fc.flag = TRUE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type, parameter, k, index))
}

for(i in paste("E", 1:6, sep="")) {
  
  try(lagResidual(n=nrow(parameter), index=i, type="100stations"))
  try(lagResidQuan(index=i, type="100stations"))
  try(StdMean.group(index=i, type="100stations"))
  try(StdMean.grouplag(index=i, parameter, type="100stations"))

}


#################################################################
##                                                             ##
##                     Dataset After 1950                      ##
##                                                             ##
#################################################################
## Create a1950 by month and by station subsets from All by station subsets
#source(file.path(local.root, "rhcode", "a1950", "rh.a1950.R"))
a1950()

## a1950 by month, interpolate missing obs by spatial loess, cross-validation
for(k in c("direct","interpolate")) {
  for(i in c(1,2)) {
    for(j in seq(0.01, 0.1, 0.005)) {
      try(interpolate(Elev = TRUE, sp=j, deg=2, Edeg=i, surf=k, fam="symmetric"))
    }
    for(j in seq(0.001, 0.009, 0.001)){
      try(interpolate(Elev = TRUE, sp=j, deg=2, Edeg=i, surf=k, fam="symmetric"))
    }
    try(crossValid(fam="symmetric", Edeg=i, surf=k, span = c(seq(0.001, 0.009, 0.001),seq(0.01, 0.1, 0.005))))
  }
}
#intpolat.visual(surf="direct")
intpolat.visual(surf="interpolate", SPsize="median")
intpolat.visual(surf="interpolate", SPsize="large")
intpolat.visual(surf="interpolate", SPsize="small")
intpolat.visual(surf="interpolate", SPsize=NULL, check=list(c(0.008, 1), c(0.015, 2)))


## Try the new code on Dec 9 2015. each month 128 sampled (from a kd-tree) locations 
## is predicted (leave-p-out cross validation, p=128) 
for(i in c(1,2)) {
  for(j in seq(0.1, 0.2, 0.05)) {
    try(newCrossValid(Elev = TRUE, sp=j, deg=2, Edeg=i, surf="direct", fam="symmetric", error="mse"))
  }
  for(j in seq(0.005, 0.1, 0.01)) {
    try(newCrossValid(Elev = TRUE, sp=j, deg=2, Edeg=i, surf="direct", fam="symmetric", error="mse"))
  }
  try(crossValidMerge(fam="symmetric", Edeg=i, surf="direct", span = c(seq(0.005, 0.1, 0.01),seq(0.1,0.2,0.05))))
}



## After found the best interpolation parameter
best <- data.frame(sp=0.015, Edeg=2, deg=2, fam="symmetric", surf="direct")
## switch the key from bymonth to by station
FileInput <- try(
  interpolateStation(sp=best$sp, Edeg=best$Edeg, deg=best$deg, fam=best$fam, surf=best$surf)
)

## tunning the STL+ parameter for stations a1950
## first sample the 128 stations from a1950 for demonstration 
source("~/Projects/Spatial/NCAR/code/kdtree/kdfindcells.R")

rst <- try(a1950.STLfit(input=FileInput, reduce=100, sw=35, sd=1, tw=231, td=2, fcw=NULL, fcd=NULL))[[1]][[2]]

a1950.fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", St.num = 128, test = TRUE)

#####################################################
##   STL experiment for tunning parameters a1950   ##
#####################################################
## Use the previous FileInput which is the bystation division of best spatial
## interpolation model.
stationSplit(reduce=200, input=FileInput, type="a1950", tn=576, valid=270)
###################
## Experiment E1 ##
###################
parameter <- expand.grid(
  sw = c(21, 30, 39), 
  tw = c(231, 313, 451), 
  td = 2, 
  sd = 1,
  fc.flag = FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type="a1950", parameter=parameter, k=k, index="E1", valid=270))
}
try(lagResidual(n=nrow(parameter), index="E1", type="a1950"))
try(lagResidQuan(index="E1", type="a1950", reduce=10))
try(StdMean.group(index="E1", type="a1950", num=36))
try(StdMean.grouplag(index="E1", parameter=parameter, type="a1950"))
try(QQDivFromNormal(index="E1", type="a1950"))
try(QQstationlag(param=parameter, index="E1", type="a1950"))
for(j in 1:nrow(parameter)) {
  rhget(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "meanerrorqqplot", "_outputs", 
      paste("QQ.error", par$dataset, "group", j, "ps", sep= ".")
    ),
    file.path(local.root, "output")
  )
}
## filtering only 128 station for visualization
try(subsetStations(index="E1", type="a1950"))
for(j in c("absmeans","means","std")) {
  errorVsLag(type="a1950", index="E1", var=c("sw","tw"), target=j)
}
for(j in c("mean.absmeans","mean.std")){
  overallErrorVsStation(type="a1950", index="E1", var=c("sw","tw"), target=j, sub=TRUE)
}
###################
## Experiment E2 ##
###################
parameter <- expand.grid(
  sw = c(21, 41, "periodic"), tw = c(123, 313, 451), td = 2, 
  sd = 1, fc.flag = FALSE, stringsAsFactors=FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type="a1950", parameter=parameter, k=k, index="E2", valid=270))
}
try(lagResidual(n=nrow(parameter), index="E2", type="a1950"))
try(lagResidQuan(index="E2", type="a1950", reduce=10))
try(StdMean.group(index="E2", type="a1950", num=36))
try(StdMean.grouplag(index="E2", parameter=parameter, type="a1950"))
try(QQDivFromNormal(index="E2", type="a1950"))
try(QQstationlag(param=parameter, index="E2", type="a1950"))
for(j in 1:nrow(parameter)) {
  rhget(
    file.path(rh.root, par$dataset, "a1950", "STLtuning", "E2", "meanerrorqqplot", "_outputs", 
      paste("QQ.error", par$dataset, "group", j, "ps", sep= ".")
    ),
    file.path(local.root, "output")
  )
}
try(subsetStations(index="E2", type="a1950"))
for(j in c("absmeans","means","std")) {
  errorVsLag(type="a1950", index="E2", var=c("sw","tw"), target=j)
}
for(j in c("mean.absmeans","mean.std")){
  overallErrorVsStation(type="a1950", index="E2", var=c("sw","tw"), target=j, sub=TRUE)
}
###################
## Experiment E3 ##
###################
parameter <- expand.grid(
  sw = "periodic", tw = c(41, 83, 123, 313, 451), td = c(1, 2), 
  sd = 1, fc.flag = FALSE, stringsAsFactors=FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type="a1950", parameter=parameter, k=k, index="E3", valid=270))
}
try(lagResidual(n=nrow(parameter), index="E3", type="a1950"))
try(lagResidQuan(index="E3", type="a1950", reduce=10))
try(StdMean.group(index="E3", type="a1950", num=36))
try(StdMean.grouplag(index="E3", parameter=parameter, type="a1950"))
try(QQDivFromNormal(index="E3", type="a1950"))
try(QQstationlag(param=parameter, index="E3", type="a1950"))
for(j in 1:nrow(parameter)) {
  rhget(
    file.path(rh.root, par$dataset, "a1950", "STLtuning", "E3", "meanerrorqqplot", "_outputs", 
      paste("QQ.error", par$dataset, "group", j, "ps", sep= ".")
    ),
    file.path(local.root, "output")
  )
}
try(subsetStations(index="E3", type="a1950"))
for(j in c("absmeans","means","std")) {
  errorVsLag(type="a1950", index="E3", var=c("td","tw"), target=j)
}
for(j in c("mean.absmeans","mean.std")){
  overallErrorVsStation(type="a1950", index="E3", var=c("td","tw"), target=j, sub=TRUE)
}
###################
## Experiment E4 ##
###################
parameter <- expand.grid(
  sw = c(21,41,"periodic"), tw = 241, td = c(1, 2), 
  sd = 1, fc.flag = FALSE, stringsAsFactors=FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type="a1950", parameter=parameter, k=k, index="E4", valid=270))
}
try(lagResidual(n=nrow(parameter), index="E4", type="a1950"))
try(lagResidQuan(index="E4", type="a1950", reduce=10))
try(StdMean.group(index="E4", type="a1950", num=36))
try(StdMean.grouplag(index="E4", parameter=parameter, type="a1950"))
try(QQDivFromNormal(index="E4", type="a1950"))
try(QQstationlag(param=parameter, index="E4", type="a1950"))
for(j in 1:nrow(parameter)) {
  rhget(
    file.path(rh.root, par$dataset, "a1950", "STLtuning", "E4", "meanerrorqqplot", "_outputs", 
      paste("QQ.error", par$dataset, "group", j, "ps", sep= ".")
    ),
    file.path(local.root, "output")
  )
}
try(subsetStations(index="E4", type="a1950"))
for(j in c("absmeans","means","std")) {
  errorVsLag(type="a1950", index="E4", var=c("sw","td"), target=j)
}
for(j in c("mean.absmeans","mean.std")){
  overallErrorVsStation(type="a1950", index="E4", var=c("sw","td"), target=j, sub=TRUE)
}
###################
## Experiment E5 ##
###################
parameter <- expand.grid(
  sw = c(11, 21, 31, 47,"periodic"), tw = 241, td = 1, 
  sd = c(1, 2), fc.flag = FALSE, stringsAsFactors=FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type="a1950", parameter=parameter, k=k, index="E5", valid=270))
}
try(lagResidual(n=nrow(parameter), index="E5", type="a1950"))
try(lagResidQuan(index="E5", type="a1950", reduce=10))
try(StdMean.group(index="E5", type="a1950", num=36))
try(StdMean.grouplag(index="E5", parameter=parameter, type="a1950"))
try(QQDivFromNormal(index="E5", type="a1950"))
try(QQstationlag(param=parameter, index="E5", type="a1950"))
for(j in 1:nrow(parameter)) {
  rhget(
    file.path(rh.root, par$dataset, "a1950", "STLtuning", "E5", "meanerrorqqplot", "_outputs", 
      paste("QQ.error", par$dataset, "group", j, "ps", sep= ".")
    ),
    file.path(local.root, "output")
  )
}
try(subsetStations(index="E5", type="a1950"))
for(j in c("absmeans","means","std")) {
  errorVsLag(type="a1950", index="E5", var=c("sw","sd"), target=j)
}
for(j in c("mean.absmeans","mean.std")){
  overallErrorVsStation(type="a1950", index="E5", var=c("sw","sd"), target=j, sub=TRUE)
}
###################
## Experiment E6 ##
###################
parameter <- expand.grid(
  sw = "periodic", tw = c(123, 241, 313), td = c(1, 2), 
  sd = 1, fc.flag = FALSE, stringsAsFactors=FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type="a1950", parameter=parameter, k=k, index="E6", valid=270))
}
try(lagResidual(n=nrow(parameter), index="E6", type="a1950"))
try(lagResidQuan(index="E6", type="a1950", reduce=10))
try(StdMean.group(index="E6", type="a1950", num=36))
try(StdMean.grouplag(index="E6", parameter=parameter, type="a1950"))
try(QQDivFromNormal(index="E6", type="a1950"))
try(QQstationlag(param=parameter, index="E6", type="a1950"))
for(j in 1:nrow(parameter)) {
  rhget(
    file.path(rh.root, par$dataset, "a1950", "STLtuning", "E6", "meanerrorqqplot", "_outputs", 
      paste("QQ.error", par$dataset, "group", j, "ps", sep= ".")
    ),
    file.path(local.root, "output")
  )
}
try(subsetStations(index="E6", type="a1950"))
for(j in c("absmeans","means","std")) {
  errorVsLag(type="a1950", index="E6", var=c("tw","td"), target=j)
}
for(j in c("mean.absmeans","mean.std")){
  overallErrorVsStation(type="a1950", index="E6", var=c("tw","td"), target=j, sub=TRUE)
}

###################
## Experiment E7 ##
###################
parameter <- expand.grid(
  sw = c(21, 39, "periodic"), 
  tw = c(123, 231, 313, 451), 
  td = 1, 
  sd = 1,
  fc.flag = FALSE, stringsAsFactors=FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36(type="a1950", parameter=parameter, k=k, index="E7blend", valid=270))
}
try(lagResidual(n=nrow(parameter), index="E7blend", type="a1950"))
try(lagResidQuan(index="E7blend", type="a1950", reduce=10))
try(StdMean.group(index="E7blend", type="a1950", num=36))
try(StdMean.grouplag(index="E7blend", parameter=parameter, type="a1950"))
try(QQDivFromNormal(index="E7blend", type="a1950"))
try(QQstationlag(param=parameter, index="E7blend", type="a1950"))
for(j in 1:nrow(parameter)) {
  rhget(
    file.path(rh.root, par$dataset, type, "STLtuning", index, "meanerrorqqplot", "_outputs", 
      paste("QQ.error", par$dataset, "group", j, "ps", sep= ".")
    ),
    file.path(local.root, "output")
  )
}
## filtering only 128 station for visualization
try(subsetStations(index="E7blend", type="a1950"))
for(j in c("absmeans","means","std")) {
  errorVsLag(type="a1950", index="E7blend", var=c("sw","tw"), target=j)
}
for(j in c("mean.absmeans","mean.std")){
  overallErrorVsStation(type="a1950", index="E7blend", var=c("sw","tw"), target=j, sub=TRUE)
}



##########################################
##      Backfitting for a1950           ##
##########################################

## backfitting iteration for three components
paramt <- data.frame(
  sw = "periodic", tw = 425, sd = 1, td = 1, fcw = 425, fcd = 1, scw = 214, scd = 2, fc.flag = TRUE
) 
backfitStart(span=0.015, family="symmetric", type="interpolate", degree=2)
backfitAll(span=0.015, family="symmetric", type="interpolate", parameter=paramt, index="E1", degree=2, inner=5, outer=1)

backfitComp(family="symmetric", type="interpolate", degree=2, span=0.015, index="E1", fc.flag=T, comp="resid")

backfitComp(family="symmetric", type="interpolate", degree=2, span=0.015, index="E1", fc.flag=T, comp="residfit")

residfit(comp = "residfit", family="symmetric", type="interpolate", degree=2, span=0.015, index="E1")