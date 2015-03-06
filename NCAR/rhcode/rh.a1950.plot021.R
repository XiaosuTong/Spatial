#################################################
## Plotting at sampled stations after 1950
## 
## - plot021:
##     dignostic plots conditional on year*month
##     all stations 7,738 are included
## - plot02:
##     dignostic plots condtional on stations for both loess02/loess04
##     sampled stations from kd-tree built based on 4,978 stations
##     have over 300 obs
##
## - loess02: is set to be without elevation
## - loess04: is set to be with elevation
#################################################
library(maps)
library(lattice)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
#par$family <- "gaussian"
#par$span <- 0.05
par$family <- "symmetric"
par$span <- 0.025
#par$loess <- "loess02.bystation.all.10pc"
par$loess <- "loess04.bystation.all.10pc"
if(pmatch(par$loess, "loess02.bystation.all.10pc", nomatch = 0)){
	title <- paste("from Robust Fit without Elev span=", par$span, sep="")
}else{
	title <- paste("from Robust Fit with Elev span=", par$span, sep="")
}
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

rst <- rhread(
	file.path(
		rh.datadir, par$dataset, "spatial", "a1950", 
		par$family, paste("sp", par$span, sep=""), par$loess
	)
)
result <- do.call("rbind", lapply(rst, "[[", 2))
result$residual <- result$tmax - result$fitted
result$station.id <- as.character(result$station.id)
result$month <- factor(
	result$month, 
	levels = c(
		"Jan","Feb","Mar","Apr","May","June",
		"July","Aug", "Sep", "Oct", "Nov", "Dec"
	)
)
## QQ plot of residuals overall
or.resid <- result[order(result$residual), "residual"]
or.resid <- or.resid[!is.na(or.resid)]
idx <- round(seq(1, length(or.resid), length.out = 1000))
resid <- data.frame(residual = or.resid[idx])
resid$f.value <- (idx - 0.5) / length(or.resid)
resid$qnorm <- qnorm(resid$f.value)

trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.qq.overall.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- xyplot( residual ~ qnorm ,
	data = resid,
	pch  = 16,
	aspect = 1,
	cex  = 0.5,
	xlab = list(label = "Unit normal quantile"),
	ylab = list(label = "Loess residuals"),
	panel = function(x, y,...) {
			panel.grid(h=-1,v=-1)
			panel.xyplot(x,y,...)
			panel.qqmathline(y, y=y,...)
	}
)
print(a)
dev.off()

##quantile plot of residuals overall
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.dist.overall.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- xyplot( residual ~ f.value ,
	data = resid,
	pch  = 16,
	aspect = 1,
	cex  = 0.5,
	xlab = list(label = "f-value"),
	ylab = list(label = "Loess residuals"),
	panel = function(x, y,...) {
			panel.grid(h=-1,v=-1)
			panel.xyplot(x,y,...)
	}
)
print(a)
dev.off()

## QQ plot of residual for each time point
## observations outside of 0.015 and 0.985 quantiles are included
Qrst <- ddply(
	.data = result,
	.variable = c("year", "month"),
	.fun = function(r) {
		r <- r[!is.na(r$residual),]
		a <- sort(r$residual)
		idx <- round(seq(1, nrow(r), length.out = 100))
		n.idx <- c(1:idx[2], idx[3:98], idx[99]:length(a))
		f.value <- (n.idx - 0.5) / length(a)
		qnorm <- qnorm(f.value)
		value <- data.frame(
			residual = a[n.idx], 
			qnorm = qnorm, 
			group = as.numeric(n.idx<idx[2] | n.idx>idx[98]),
			fv = f.value
		)
	}
)
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.dist.bytime.",
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in 1950:1997){
a <- xyplot(residual ~ qnorm | month*as.factor(year),
	data = Qrst, 
	subset = year == paste(i),
	layout = c(12,1),
	group = group,
  col = col[1:2],
  distribute.type= TRUE,
	pch  = 16,
	#aspect = 1,
	cex  = 0.3,
#	scale = list(y=list(relation="free")),
	xlab = list(label = "Unit normal quantile"),
	ylab = list(label = "Loess Residuals"),
	main = list(
		label = paste("Normal Quantiles of Residuals", title)
	),
#	prepanel = prepanel.qqmathline,
	panel = function(x, y, ...) {
		#panel.grid()
		panel.abline(h=0, lwd=0.5, col="black")
		panel.xyplot(x,y,...)
	}
)
print(a)
}
dev.off()
## bytime2 is only the quantiles between 0.015 and 0.985 quantiles
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.dist.bytime2.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- xyplot(residual ~ qnorm | month*as.factor(year),
	data = Qrst, 
	subset = group == 0,
	layout = c(6,2),
	pch  = 16,
#	aspect = 1,
	cex  = 0.3,
	xlab = list(label = "Unit normal quantile"),
	ylab = list(label = "Loess Residuals"),
	main = list(
		label = paste("Normal Quantiles of Residuals", title)
	),
	panel = function(x,y, ...) {
		panel.qqmathline(y,y=y,...)
		panel.abline(h=0, lwd=0.5, col="black")
		panel.xyplot(x,y,...)
	}
)
print(a)
dev.off()
## scatter plot of residual against fited value conditional on month
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.vs.fit.bytime.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
a <- xyplot(residual ~ fitted | month*as.factor(year),
	data = result,
	layout = c(4,3),
  key=list(
    text = list(label=c(
      "residuals",
      "loess smoothing: span=0.75, degree=2 "
    )),
    lines = list(
      pch=c(".",""), 
      cex=4, 
      lwd=1.5, 
      type=c("p","l"), 
      col=col[1:2]
    ), 
    columns = 2
  ),
	pch  = 16,
	cex  = 0.2,
	scale = list(
		x=list(relation = "free"),
		y=list(relation = "free")
	),
	ylab = list(label = "Loess Residuals"),
	xlab = list(label = "Loess Fitted Value"),
	panel = function(x, y,...) {
		panel.abline(h=0, lwd=0.5, col="black")
		panel.xyplot(x, y, ...)
		#panel.abline(h=0, col="black", lwd = 0.5, lty=1)
		panel.loess(x,y,span=0.75,degree=2,col=col[2],family="symmetric",...)
	}
)
print(a)
dev.off()

## scatter plot of residual against lon/lat conditional on interval of lat/lon
trellis.device(
	postscript, 
	file = paste(
		local.output, "/a1950.loess.resid.vs.lon.bytime.", 
		par$dataset, ".ps", sep = ""
	), 
	color = TRUE,
	paper = "legal"
)
for(i in sort(unique(result$year))) {
#	for(j in levels(result$month)) {
		a <- xyplot(residual ~ lon | equal.count(lat, 20, overlap=0),
			data = subset(result, year == i),
			layout = c(10,2), #for vs.lon layout is c(10,2)
			strip=strip.custom(var.name = "Latitude", strip.levels=rep(FALSE, 2)),
			pch  = 16,
			cex  = 0.3,
			scale = list(
				y = list(
					relation = "same", 
					alternating = TRUE
				),
				x = list(
					relation = "free",
					tick.number = 3
				)
			),
    	key=list(
        text = list(label=c(
        	"residuals",
          "loess smoothing: span=0.75, degree=2 "
        )),
        lines = list(
            pch=c(".",""), 
            cex=4, 
            lwd=1.5, 
            type=c("p","l"), 
            col=col[1:2]
        ), 
        columns = 2
    	),
  		main = list(label = paste("Year", i)),
			ylab = list(label = "Loess Residuals"),
			xlab = list(label = "Latitude"),
			prepanel = function(x,y,...) {
				prepanel.loess(x, y, span = 3/4, degree=2)
			},
			panel = function(x, y,...) {
				panel.abline(h=0, col="black")
				panel.xyplot(x, y,...)
				panel.loess(x,y,span = 3/4, degree=2, col=col[2],...)
			}
		)
		print(a)
#	}
}
dev.off()
