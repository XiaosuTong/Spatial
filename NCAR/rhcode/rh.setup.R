#Set up the directory and load the data.
library(lattice)
library(plyr)
local.datadir <- "/home/tongx/Projects/Spatial/NCAR/RData"
local.output <- "/home/tongx/Projects/Spatial/NCAR/output"
rh.datadir <- "/wsc/tongx/Spatial/tmp"
#The output directory for saving plots should have all permission since the plots are written to the output directory by a user related to hadoop.
#In the reduce step, the permission for the plots should be changed to be all +wrx.
rh.output <- "/wsc/tongx/Spatial/output"
lib.loc <- "/home/tongx/R_LIBS"
