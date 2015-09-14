#Set up the directory and load the data.
library(lattice)
library(plyr)

source("~/Rhipe/ross.initial.R")

par <- list()
par$dataset <- "tmax"
par$Machine <- "rossmann"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

if(par$dataset == "precip") {
  Nstations <- 11918
} else {
  Nstations <- 8125
}

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

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
STLfit(sw=77, sd=1, tw=495, td=2, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t495td2_s77sd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)

STLfit(sw=103, sd=1, tw=617, td=2, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t617td2_s77sd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)

STLfit(sw="periodic", sd=1, tw=1141, td=2, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1141td2_speriodicsd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)

STLfit(sw="periodic", sd=1, tw=1141, td=1, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1141td1_speriodicsd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)

STLfit(sw="periodic", sd=1, tw=241, td=1, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t241td1_speriodicsd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)

STLfit(sw="periodic", sd=1, tw=1855, td=1, fcw=121, fcd=2)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1855td1_speriodicsd1_f121fd2"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = TRUE)

STLfit(sw="periodic", sd=1, tw=1855, td=1, fcw=241, fcd=1)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1855td1_speriodicsd1_f241fd1"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = TRUE)

STLfit(sw=103, sd=2, tw=1141, td=1, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1141td1_s103sd2_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)


STLfit(sw=77, sd=2, tw=1855, td=1, fcw=241, fcd=1)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1855td1_s77sd2_f241fd1"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = TRUE)


##################################################################
##    STL tunning fit for the 100 stations with full obs        ##
################################################################## 
stationSplit(reduce=200)
index <- "E1"
parameter <- expand.grid(
  sw = c(51, 73, 93), tw = c(617, 865, 1113), td = 2, sd = 1
)
for(k in 1:nrow(parameter)) {
  predict36(parameter, k, index)
}
lagResidual(n=nrow(parameter), index="E1")
lagResidQuan(index="E1")
StdMean.group(index="E1")
StdMean.grouplag(index="E1", parameter=parameter)

index <- "E4"
parameter <- expand.grid(
  sw = c("periodic"), tw = c(121, 241, 361, 751, 1141), 
  td = c(1,2), sd = 1, stringsAsFactors = FALSE
)
for(k in 1:nrow(parameter)) {
  predict36(parameter, k, index)
}
lagResidual(n=nrow(parameter), index="E4")
lagResidQuan(index="E4")
StdMean.group(index="E4")
StdMean.grouplag(index="E4", parameter=parameter)
#################################################################
##                     Dataset After 1950                      ##
#################################################################
## Create a1950 by month and by station subsets from All by station subsets
#source(file.path(local.root, "rhcode", "a1950", "rh.a1950.R"))
a1950()
## a1950 by month, interpolate missing obs by spatial loess, cross-validation
for(k in c("interpolate","direct")) {
  for(i in c(1,2)) {
    for(j in seq(0.01, 0.1, 0.005)) {
      interpolate(Elev = TRUE, sp=j, deg=2, Edeg=i, surf=k, fam="symmetric")
    }
    crossValid(fam="symmetric", Edeg=i, surf=k)
  }
}






