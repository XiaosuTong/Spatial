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

## Rhipe job creating data.frame of 100 stations  ##
#source(file.path(local.root, "rhcode", "s100", "rh.100stations.R"))
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "aggregated"))[[1]][[2]]
## create plots for raw observations  ##
#source(file.path(local.root, "code", "raw.visual", "100stations.plot.R"))
scatterRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
monthRaw(data=rst, outputdir=file.path(local.root, "output"), target = "tmax", size = "letter", test = F)
monthQQ(data=rst, outputdir=file.path(local.root, "output"), target = "tmax", size = "letter", test = F)
monthSpatial(data=rst, outputdir=file.path(local.root, "output"), target = "tmax", vari="lon", size = "letter")
monthSpatial(data=rst, outputdir=file.path(local.root, "output"), target = "tmax", vari="lat", size = "letter")
monthSpatial(data=rst, outputdir=file.path(local.root, "output"), target = "tmax", vari="elev", size = "letter")

## Rhipe job running stl2 on each station of 100 stations ##
#source(file.path(local.root, "rhcode", "s100", "rh.100stations.stl.R"))
#source(file.path(local.root, "code", "raw.visual", "100stations.stl2.R"))

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

STLfit(sw=periodic, sd=1, tw=1141, td=2, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1141td2_speriodicsd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)

STLfit(sw=periodic, sd=1, tw=1141, td=1, fcw=NULL, fcd=NULL)
rst <- rhread(file.path(rh.root, par$dataset, "100stations", "STL", "t1141td1_speriodicsd1_ffd"))[[1]][[2]]
fitRaw(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test = F)
remainderDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
seasonalDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F)
trendDiag(data=rst, outputdir=file.path(local.root, "output"), target="tmax", size = "letter", test=F, fc = FALSE)


## Create a1950 by month and by station subsets from All by station subsets
#source(file.path(local.root, "rhcode", "a1950", "rh.a1950.R"))
a1950()
## a1950 by month, interpolate missing obs by spatial loess
for(i in seq(0.01, 0.1, 0.005)) {
  interpolate(Elev = TRUE, sp=i, deg=2, fam="symmetric")
}
crossValid(fam="symmetric")
rst <- rhread("/wsc/tongx/Spatial/tmp/tmax/a1950/bymonth.fit/symmetric/MSE")[[1]][[2]]
