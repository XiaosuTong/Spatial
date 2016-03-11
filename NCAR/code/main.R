#Set up the directory and load the data.
library(lattice)
library(plyr)
library(hexbin)
library(maps)
library(magrittr)
library(gtools)
us.map <- map('state', plot = FALSE, fill = TRUE)
options(stringsAsFactors=FALSE)
options(java.parameters = "-Xmx1024m")

#source("~/Rhipe/ross.initial.R")
source("~/Rhipe/rhinitial.R")

par <- list()
par$dataset <- "tmax"
par$Machine <- "wsc"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")


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
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t617td2_s103sd1_ffd"))[[1]][[2]]
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
FileInput <- file.path(rh.root, par$dataset, "a1950", "bymonth")
FileInput <- bymonthSplit(input=FileInput, leaf = 100, vari="resp")

#########################################
## Imputation of missing w/o elevation ##
#########################################
## For each month, calculate the range of each spatial factor in the 
## equal cut of 20 intervals.


## E1 span = 0.05, degree=2, without Elevation
para <- list(span=0.05, Edeg=0)
try(interpolate(sp=para$span, deg=2, Edeg=para$Edeg, surf="direct", fam="symmetric"))
## Check the residual of the spatial loess imputing
a1950.spaImputeVisual(family="symmetric", surf="direct", Edeg=para$Edeg, span=para$span)

## E2 span = 0.025, degree=2, without Elevation
para <- list(span=0.025, Edeg=0)
try(interpolate(sp=para$span, deg=2, Edeg=para$Edeg, surf="direct", fam="symmetric"))
## Check the residual of the spatial loess imputing
a1950.spaImputeVisual(family="symmetric", surf="direct", Edeg=para$Edeg, span=para$span)

## E3 span = 0.015, degree=2, without Elevation
para <- list(span=0.015, Edeg=0)
try(interpolate(sp=para$span, deg=2, Edeg=para$Edeg, surf="direct", fam="symmetric"))
## Check the residual of the spatial loess imputing
a1950.spaImputeVisual(family="symmetric", surf="direct", Edeg=para$Edeg, span=para$span)

## E4 span = 0.005, degree=2, without Elevation
para <- list(span=0.005, Edeg=0)
try(interpolate(sp=para$span, deg=2, Edeg=para$Edeg, surf="direct", fam="symmetric"))
## Check the residual of the spatial loess imputing
a1950.spaImputeVisual(family="symmetric", surf="direct", Edeg=para$Edeg, span=para$span)

## Try the new code for Cross Validation. Each month 128 sampled 
## (from a kd-tree) locations is predicted (leave-p-out cross validation, p=128)
FileInput <- file.path(rh.root, par$dataset, "a1950", "bymonthSplit")
for(i in c(0)) {
  for(j in seq(0.1, 0.2, 0.05)) {
    try(newCrossValid(input=FileInput, vari="resp", sp=j, deg=2, Edeg=i, surf="direct", fam="symmetric", error="mse"))
  }
  for(j in c(0.001, seq(0.005, 0.1, 0.01))) {
    try(newCrossValid(input=FileInput, vari="resp", sp=j, deg=2, Edeg=i, surf="direct", fam="symmetric", error="mse"))
  }
  for(j in seq(0.002, 0.004, 0.001)) {
    try(newCrossValid(input=FileInput, vari="resp", sp=j, deg=2, Edeg=i, surf="direct", fam="symmetric", error="mse"))
  }
  try(crossValidMerge(
    input=file.path(rh.root, par$dataset, "a1950", "bymonth.fit.cv", "symmetric", "direct", i), 
    span = c(seq(0.001,0.004, 0.001), seq(0.005, 0.1, 0.01),seq(0.1,0.2,0.05))
  ))
}

imputeCrossValid(input=file.path(rh.root, par$dataset, "a1950", "bymonth.fit.cv", "symmetric", "direct"), Edeg = FALSE)


########################################
## Imputation of missing w/ elevation ##
########################################
for(ii in c(1,2)) {
  ## E1 span = 0.05, degree=2, without Elevation
  para <- list(span=0.05, Edeg=ii)
  try(interpolate(sp=para$span, deg=2, Edeg=para$Edeg, surf="direct", fam="symmetric"))
  ## Check the residual of the spatial loess imputing
  a1950.spaImputeVisual(family="symmetric", surf="direct", Edeg=para$Edeg, span=para$span)  

  ## E2 span = 0.025, degree=2, without Elevation
  para <- list(span=0.025, Edeg=ii)
  try(interpolate(sp=para$span, deg=2, Edeg=para$Edeg, surf="direct", fam="symmetric"))
  ## Check the residual of the spatial loess imputing
  a1950.spaImputeVisual(family="symmetric", surf="direct", Edeg=para$Edeg, span=para$span)  

  ## E3 span = 0.015, degree=2, without Elevation
  para <- list(span=0.015, Edeg=ii)
  try(interpolate(sp=para$span, deg=2, Edeg=para$Edeg, surf="direct", fam="symmetric"))
  ## Check the residual of the spatial loess imputing
  a1950.spaImputeVisual(family="symmetric", surf="direct", Edeg=para$Edeg, span=para$span)  

  ## E4 span = 0.005, degree=2, without Elevation
  para <- list(span=0.005, Edeg=ii)
  try(interpolate(sp=para$span, deg=2, Edeg=para$Edeg, surf="direct", fam="symmetric"))
  ## Check the residual of the spatial loess imputing
  a1950.spaImputeVisual(family="symmetric", surf="direct", Edeg=para$Edeg, span=para$span)
}

## Try the new code for Cross Validation. Each month 128 sampled 
## (from a kd-tree) locations is predicted (leave-p-out cross validation, p=128) 

for(i in c(1, 2)) {
  for(j in seq(0.1, 0.2, 0.05)) {
    try(newCrossValid(input=FileInput, vari="resp", sp=j, deg=2, Edeg=i, surf="direct", fam="symmetric", error="mse"))
  }
  for(j in seq(0.005, 0.1, 0.01)) {
    try(newCrossValid(input=FileInput, vari="resp", sp=j, deg=2, Edeg=i, surf="direct", fam="symmetric", error="mse"))
  }
  try(crossValidMerge(
    input=file.path(rh.root, par$dataset, "a1950", "bymonth.fit.cv", "symmetric", "direct", i), 
    span = c(seq(0.005, 0.1, 0.01),seq(0.1,0.2,0.05))
  ))
}

imputeCrossValid(input=file.path(rh.root, par$dataset, "a1950", "bymonth.fit.cv", "symmetric", "direct"), Edeg = TRUE)


##############################################
##  Visualize the stlplus fit based on the  ##
##  best model                              ##
##############################################
## first sample the 512 stations from a1950 for demonstration 
source("~/Projects/Spatial/NCAR/code/kdtree/kdfindcells.R")
## the best imputed model
FileInput <- file.path(rh.root, par$dataset, "a1950", "bymonth.fit", "symmetric", "direct", 2, "sp0.015")
FileOutput <- file.path(rh.root, par$dataset, "a1950", "bymonth.fit", "symmetric", "direct", 2, "sp0.015.bystation")
## switch the key from bymonth to by station
swapTostation(FileInput, FileOutput)
FileInput <- FileOutput

stlplusFitVisl<- function(fileinput=FileInput, SW, SD, TW, TD, FCW=NULL, FCD=NULL) {

  paras <- list(sw=SW, sd=SD, tw=TW, td=TD, fcw=FCW, fcd=FCD)
  input <- a1950.STLfit(fileinput, reduce=72, tuning=paras)
  a1950.STLvisual(
    paras=paras, input=input, plotEng=plotEng.raw, 
    name="stlraw.vs.time", sample = TRUE, multiple=NULL
  )
  a1950.STLvisual(
    paras=paras, input=input, plotEng=plotEng.trend, 
    name="trend.vs.time", sample = TRUE, multiple=c(4,3)
  )
  a1950.STLvisual(
    paras=paras, input=input, plotEng=plotEng.periodicseason, 
    name="seasonal.vs.year", sample = TRUE, multiple=c(5,4)
  )
  a1950.STLvisual(
    paras=paras, input=input, plotEng=plotEng.searemainder, 
    name="searemainder.vs.year", sample = TRUE, multiple=NULL
  )
  a1950.STLvisual(
    paras=paras, input=input, plotEng=plotEng.QQremaider, 
    name="QQ.remainder", sample = TRUE, multiple=c(5,3)
  )
  a1950.STLvisual(
    paras=paras, input=input, plotEng=plotEng.remainderMonth, 
    name="remainder.vs.year", sample = TRUE, multiple=NULL
  )
  a1950.STLvisual(
    paras=paras, input=input, plotEng=plotEng.remainderMonth2, 
    name="remainder.vs.year2", sample = TRUE, multiple=NULL
  )
  a1950.STLvisual(
    paras=paras, input=input, plotEng=plotEng.remainderACF, 
    name="remainder.acf", sample = TRUE, multiple=c(3,3)
  )
  a1950.STLvisual(
    paras=paras, input=input, plotEng=plotEng.remainderDate, 
    name="remainder.vs.time", sample = TRUE, multiple=NULL
  )

}

stlplusFitVisl(fileinput=FileInput, SW=21, SD=1, TW=231, TD=2, FCW=NULL, FCD=NULL)

stlplusFitVisl(fileinput=FileInput, SW=41, SD=1, TW=241, TD=1, FCW=NULL, FCD=NULL)

stlplusFitVisl(fileinput=FileInput, SW=45, SD=2, TW=451, TD=1, FCW=NULL, FCD=NULL)

stlplusFitVisl(fileinput=FileInput, SW="periodic", SD=1, TW=241, TD=1, FCW=NULL, FCD=NULL)

stlplusFitVisl(fileinput=FileInput, SW="periodic", SD=1, TW=451, TD=1, FCW=NULL, FCD=NULL)

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
    file.path(rh.root, par$dataset, "a1950", "STLtuning", "E1", "meanerrorqqplot", "_outputs", 
      paste("QQ.error", par$dataset, "group", j, "ps", sep= ".")
    ),
    file.path(local.root, "output")
  )
}
## filtering only 512 station for visualization
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
  sw = c(11, 41, "periodic"), tw = c(123, 241, 451), td = 2, 
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
  sw = "periodic", tw = c(41, 83, 123, 241, 451), td = c(1, 2), 
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
#####################################
## Experiment E5, on the web is E4 ##
#####################################
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
  errorVsLag(type="a1950", index="E5", var=c("sw","sd"), target=j, sd.rm=TRUE)
}
for(j in c("mean.absmeans","mean.std")){
  overallErrorVsStation(type="a1950", index="E5", var=c("sw","sd"), target=j, sub=TRUE, sd.rm=TRUE)
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
  try(predict36.blend(type="a1950", parameter=parameter, k=k, index="E7blend", valid=270))
}
try(lagResidual(n=nrow(parameter), index="E7blend", type="a1950"))
try(lagResidQuan(index="E7blend", type="a1950", reduce=10))
try(StdMean.group(index="E7blend", type="a1950", num=36))
try(StdMean.grouplag(index="E7blend", parameter=parameter, type="a1950"))
try(QQDivFromNormal(index="E7blend", type="a1950"))
try(QQstationlag(param=parameter, index="E7blend", type="a1950"))
for(j in 1:nrow(parameter)) {
  rhget(
    file.path(rh.root, par$dataset, "a1950", "STLtuning", "E7blend", "meanerrorqqplot", "_outputs", 
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

###################
## Experiment E8 ##
###################
parameter <- expand.grid(
  sw = c(7, 11, 21, "periodic"), 
  tw = 241, 
  td = 1, 
  sd = c(1, 2),
  fc.flag = FALSE, stringsAsFactors=FALSE
)
for(k in 1:nrow(parameter)) {
  try(predict36.blend(type="a1950", parameter=parameter, k=k, index="E8blend", valid=270))
}
try(lagResidual(n=nrow(parameter), index="E8blend", type="a1950"))
try(lagResidQuan(index="E8blend", type="a1950", reduce=10))
try(StdMean.group(index="E8blend", type="a1950", num=36))
try(StdMean.grouplag(index="E8blend", parameter=parameter, type="a1950"))
try(QQDivFromNormal(index="E8blend", type="a1950"))
try(QQstationlag(param=parameter, index="E8blend", type="a1950"))
for(j in 1:nrow(parameter)) {
  rhget(
    file.path(rh.root, par$dataset, "a1950", "STLtuning", "E8blend", "meanerrorqqplot", "_outputs", 
      paste("QQ.error", par$dataset, "group", j, "ps", sep= ".")
    ),
    file.path(local.root, "output")
  )
}
## filtering only 128 station for visualization
try(subsetStations(index="E8blend", type="a1950"))
for(j in c("absmeans","means","std")) {
  errorVsLag(type="a1950", index="E8blend", var=c("sw","sd"), target=j)
}
for(j in c("mean.absmeans","mean.std")){
  overallErrorVsStation(type="a1950", index="E8blend", var=c("sw","sd"), target=j, sub=TRUE)
}


#############################################################
##  visualize remainder of stlplus against spatial factors ##
##  longitude, latitude, and elevation                     ##
#############################################################
bestStlplus <- "t241td1_speriodicsd1_ffd"
FileInput <- file.path(rh.root, par$dataset, "a1950", "bymonth.fit", "symmetric", "direct", 2, "sp0.015.bystation")
FileOutput <- file.path(rh.root, par$dataset, "a1950", "STL.bymonth", bestStlplus)

swapTomonth(FileInput, FileOutput)
FileInput <- FileOutput

spaPara <- data.frame(permutations(3, 2, c("lon","lat","elev")))
for (k in 1:nrow(spaPara)) {
  vars <- spaPara[k, ]
  a1950.spafitVisualMon(input=FileInput, plotEng.RevsSpa, vars, target="remainder")
}

###################################################
##  Spatial fit on the remainder of stlplus      ##
##  parallelly visualize the residual against    ##
##  three spatial factors for each month.        ##
##  Then visualize the residual against time,    ##
##  and finally visualize the overall quantiles  ##
##  of the residual                              ##
###################################################
## find all stations after 1950 without missing value, save all station.id
## as a1950Nomiss vector in nomiss.a1950.RData
a1950.Nomiss(file.path(rh.root, par$dataset, "a1950", "bystation"))

for (ii in c(0.05,0.025,0.015, 0.005)) {
  for (jj in c(1, 2)) {
    residSpaFitVisl(i=ii, j=jj, bestStlplus="t241td1_speriodicsd1_ffd")
  }
}


#  df <- a1950.residQuant(
#    input=FileInput, target="residual", by="coastdis", 
#    probs=seq(0.005, 0.995, 0.005), nBins = 10000, tails = 100, coast=TRUE, dislim=300
#  )
#  trellis.device(
#    device = postscript, 
#    file = file.path(local.root, "output", "tmp.ps"),
#    color = TRUE, 
#    paper = "letter"
#  )
#    b <- xyplot(q ~ fval
#      , data = df
#      , group = coastdis
#      , xlab = list(label="f-value", cex=1.5)
#      , ylab = list(label="Residual", cex=1.5)
#      , scale = list(cex=1.2)
#    )
#    print(b)
#  dev.off()  

############################################
##  Cross validation for the spatial fit  ##
############################################
bestStlplus <- "t241td1_speriodicsd1_ffd"
FileInput <- file.path(rh.root, par$dataset, "a1950", "STL.bymonth", bestStlplus)
FileInput <- bymonthSplit(input=FileInput, leaf = 100, vari="remainder")

for(i in c(1, 2)) {
  for(j in seq(0.1, 0.2, 0.05)) {
    try(newCrossValid(input=FileInput, vari="remainder", sp=j, deg=2, Edeg=i, surf="direct", fam="symmetric", error="mse"))
  }
  for(j in seq(0.005, 0.1, 0.01)) {
    try(newCrossValid(input=FileInput, vari="remainder", sp=j, deg=2, Edeg=i, surf="direct", fam="symmetric", error="mse"))
  }
  try(crossValidMerge(
    input=file.path(rh.root, par$dataset, "a1950", "STL.bymonth", "t241td1_speriodicsd1_ffd.fit.cv", "symmetric", "direct", i), 
    span = c(seq(0.005, 0.085, 0.01),seq(0.1,0.2,0.05))
  ))
}

imputeCrossValid(
  input=file.path(rh.root, par$dataset, "a1950", "STL.bymonth", "t241td1_speriodicsd1_ffd.fit.cv", "symmetric", "direct"), 
  Edeg = TRUE
)


##############################################
## Visualize all outliers after spatial fit ##
##############################################
bestStlplus <- "t241td1_speriodicsd1_ffd"
FileInput <- file.path(
  rh.root, par$dataset, "a1950", "STL.bymonth.remaindfit", 
  bestStlplus, "symmetric", "direct", "2", "sp0.015"
)
FileOutput <- paste(FileInput, "outliers", sep=".")
## generate the outliers.a1950.RData object which 
## is a vector of station.id
outliersStatAll(FileInput, FileOutput, lim=2)
## generate the dataframe including all outliers
## There are 2,096,822 observations, and 21,471 are outliers
outliers <- outliersTotal(FileInput, FileOutput, lim=2)
outlierBymonth <- outliersCount(FileInput, FileOutput, lim=2, by="month")
outlierBymonth$mFlag <- outlierBymonth$month %in% c("Jan","Feb","Dec")
outliers$mFlag <- outliers$month %in% c("Jan","Feb","Dec")

outlierVisAll(outliers)
outlierVisAll(outliers, byname="Winter", byvari="mFlag")
outlierVisMonth(outlierBymonth)
outlierVisMonth(outlierBymonth, byname="Winter", byvari="mFlag")

## generate the count of outliers in each station
FileInput <- file.path(
  rh.root, par$dataset, "a1950", "STL.bymonth.remaindfit", 
  bestStlplus, "symmetric", "direct", "2", "sp0.015.bystation"
)
FileOutput <- paste(FileInput, "outliers", sep=".")
outlierBystation <- outliersCount(FileInput, FileOutput, lim=2, by="station")

outlierVisStat(outlierBystation)
## generate the stations which has more than 2^5 outliers or has outliers larger than 10
## and save the on HDFS as outliersTop.a1950.RData
outliersTop(FileInput, FileOutput, lim=2, top=2^5, ORtop=10, by="station")
## visualize the raw and outliers against time for each station from outliersTop.a1950.RData
paras <- list(sw="periodic", sd=1, tw=241, td=1, fcw=NULL, fcd=NULL)
a1950.STLvisual(paras, FileInput, plotEng.rawOutlier, name="stlOutlier.vs.time", sample = FALSE, Alloutlier=FALSE, multiple=NULL)

FileInput <- file.path(
  rh.root, par$dataset, "a1950", "STL.bymonth.remaindfit", 
  bestStlplus, "symmetric", "direct", "2", "sp0.015"
)
FileOutput <- paste(FileInput, "outliers", sep=".")
outliersTop(FileInput, FileOutput, lim=2, top=2^6, ORtop=50, by="month")
##visualize the outliers location in space for each month
outlierMonth(FileInput, FileOutput, name="outlierLocbyMonth", plotEng=plotEng.outlierMonth)




## find out the number of outliers that remainder is small but spafit is very large
## totally there are 3805 such outliers from 1960 stations
tmp <- outliersMissCount(FileInput, FileOutput, lim=2)
dim(tmp)


xyplot(radius ~ elev2, data = job.mr[[1]][[2]])






####################################
##    Duplicate each time series  ##
####################################
FileInput <- file.path(rh.root, par$dataset, "a1950", "bymonth.fit", "symmetric", "direct", 2, "sp0.015")
FileOutput <- file.path(rh.root, par$dataset, "simulate", "bystation.orig")
swapTostation(FileInput, FileOutput)

FileInput <- file.path(rh.root, par$dataset, "simulate", "bystation.orig")
FileOutput <- file.path(rh.root, par$dataset, "simulate", "bystation")
system.time(repTime(FileInput, FileOutput, buffSize=10000, Rep=5800))

FileInput <- file.path(rh.root, par$dataset, "simulate", "bystation")
FileOutput <- file.path(rh.root, par$dataset, "simulate", "bymonth")

rst <- data.frame(large=c(0,0,0), small=c(0,0,0), rep = c(1,2,3))
for (i in 1:3) {
  
  rst[i, 1] <- system.time(swapTomonth(FileInput, FileOutput, io_sort=400, spill_percent=0.9, parallelcopies=5, reduce_merge_inmem=1000, reduce_input_buffer=0.0))[3]
  Sys.sleep(60)
  rst[i, 2] <- system.time(swapTomonth(FileInput, FileOutput, io_sort=10, spill_percent=0.6, parallelcopies=5))[3]
  Sys.sleep(60)

}
#   large   small rep
# 295.755 339.719   1
# 309.155 330.587   2
# 309.345 326.252   3

rst <- data.frame(large=c(0,0,0), small=c(0,0,0), rep = c(1,2,3))
for (i in 1:3) {
  
  rst[i, 1] <- system.time(swapTomonth(FileInput, FileOutput, io_sort=400, spill_percent=0.8, parallelcopies=5, reduce_merge_inmem=1000, reduce_input_buffer=0.0))[3]
  Sys.sleep(60)
  rst[i, 2] <- system.time(swapTomonth(FileInput, FileOutput, io_sort=400, spill_percent=0.8, parallelcopies=5, reduce_merge_inmem=0, reduce_input_buffer=1))[3]
  Sys.sleep(60)
  
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