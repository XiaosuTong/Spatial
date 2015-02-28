#############################################
##spatial loess fit comparasion after 1950 ##
#############################################

## compare the loess fit without elev with fit without elev
## so compare loess02 with loess04
library(maps)
library(lattice)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
#par$family <- "gaussian"
par$span <- 0.05
par$family <- "symmetric"
#par$span <- 0.025
par$loess1 <- "loess02.bystation.all.10pc"
par$loess2 <- "loess04.bystation.all.10pc"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

rst1 <- rhread(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), par$loess1
	)
)
result1 <- do.call("rbind", lapply(rst1, "[[", 2))
result1$residual <- result1$tmax - result1$fitted
result1$station.id <- as.character(result1$station.id)
result1$month <- factor(
	result1$month, 
	levels = c(
		"Jan","Feb","Mar","Apr","May","June",
		"July","Aug", "Sep", "Oct", "Nov", "Dec"
	)
)
rm(rst1)
rst2 <- rhread(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), par$loess2
	)
)
result2 <- do.call("rbind", lapply(rst2, "[[", 2))
result2$residual <- result2$tmax - result2$fitted
result2$station.id <- as.character(result2$station.id)
result2$month <- factor(
	result2$month, 
	levels = c(
		"Jan","Feb","Mar","Apr","May","June",
		"July","Aug", "Sep", "Oct", "Nov", "Dec"
	)
)
rm(rst2)
result <- cbind(
	result1[order(result1$station.id, result1$time),], 
	residual2 = result2[order(result2$station.id, result2$time),]$residual
)
rm(result1, result2)

trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.elev.vs.noelev.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- xyplot( (residual-residual2) ~ (residual+residual2)/2 | month*factor(year),
	data = result,
#	subset = year == "1950",
	pch  = 16,
	layout = c(3,2),
	aspect = 1,
	cex  = 0.3,
	scale = list(
		y = list(relation = "free"),
		x = list(relation = "free")
	),
	main = paste("Residuals from Robust Fit span=", par$span, sep=""),
	xlab = list(label = "Mean"),
	ylab = list(label = "Difference"),
	panel = function(x, y,...) {
			panel.abline(h=0, lwd=0.5)
			panel.xyplot(x,y,...)
	}
)
print(a)
dev.off()

load(file.path(local.datadir, "samplestation.a1950.RData"))
sub <- result[with(result, station.id %in% tmax.sample.a1950$station.id),]
sub$station.id <- factor(
	sub$station.id, 
	levels = tmax.sample.a1950[order(tmax.sample.a1950$leaf),1]
)


trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.elev.vs.noelev.bystation.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- xyplot( (residual-residual2) ~ (residual+residual2)/2 | factor(station.id, labels = 1:128),
	data = sub,
	pch  = 16,
	layout = c(4,3),
#	aspect = 1,
	cex  = 0.3,
	main = paste("Residuals from Robust Fit span=", par$span, sep=""),
	scale = list(y=list(relation="free"),x=list(relation="free")),
	xlab = list(label = "Mean"),
	ylab = list(label = "Difference"),
	panel = function(x, y,...) {
			panel.abline(h=0, lwd=0.5)
			panel.xyplot(x,y,...)
	}
)
print(a)
dev.off()