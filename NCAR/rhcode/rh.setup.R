#Set up the directory and load the data.
library(lattice)
library(plyr)
if(par$Machine == "rossmann"){
	local.root <- "/home/tongx/Projects/Spatial/NCAR"
	rh.root <- "/wsc/tongx/Spatial/tmp"
	#The output directory for saving plots should have all permission since the plots are written to 
    #the output directory by a user related to hadoop.
	#In the reduce step, the permission for the plots should be changed to be all +wrx.
	rh.output <- "/wsc/tongx/Spatial/output"
	lib.loc <- "/home/tongx/R_LIBS"
}else if(par$Machine == "gacrux" | par$Machine == "adhara"){
	local.root <- "/home/shaula/u16/tongx/Projects/Spatial/NCAR"
	rh.root <- "/ln/tongx/Spatial/tmp"
	local.raw <- "/home/shaula/u16/tongx/Projects/Spatial/NCAR/Raw"
	#The output directory for saving plots should have all permission since the plots are written to 
    #the output directory by a user related to hadoop.
	#In the reduce step, the permission for the plots should be changed to be all +wrx.
	rh.output <- "/ln/tongx/Spatial/output"
	lib.loc <- "/home/shaula/u16/tongx/R_LIBS"
}
