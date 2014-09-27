#########################################################
##ploting code for spatial loess fit at the 100 stations
#########################################################
##loess01 is set to be span=0.125, degree=2
##loess02 is set to be span=0.06, degree=2
source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$loess <- "loess02"
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
load(paste(
	"~/Projects/Spatial/NCAR/RData/", 
	par$dataset,
	".100stations.stl.RData", 
	sep="")
)
rst <- rhread(file.path(rh.datadir, par$dataset, "spatial", "100stations", par$loess))
result <- do.call("rbind", lapply(rst, "[[", 2))
result$factor <- factor(
	rep(rep(paste("Period", 1:9), c(rep(144,8),84)), times = 100), 
	levels = paste("Period", c(9:1))
)
if(par$dataset == "precip") {
	ylab <- "Precipitation (millimeters)"
}else if(par$dataset == "tmax") {
	ylab <- "Maximum Temperature (degree centigrade)"
}else {
	ylab <- "Minimum Temperature (degree centigrade)"
}
result$time2 <- c(rep(0:143,8), 0:83)

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

order <- ddply(
	.data = result,
  .variables = "station.id",
  .fun = summarise,
  mean = mean(fitted)
)
order.st <- as.character(
	order[order(order$mean, decreasing = TRUE), ]$station.id
)
result$station.id <- factor(
	result$station.id, 
	levels = order.st
)
result <- result[with(result, order(station.id, time)),]
stations100.stl <- stations100.stl[with(stations100.stl, order(station.id, date)),]
result$remainder <- stations100.stl$fc.remainder
trellis.device(
	postscript, 
	file = paste(local.output, "/scatterplot_of_resid.vs.remainder_", par$dataset, ".ps", sep = ""), 
	color = TRUE,
	paper = "legal"
)
	b <- xyplot( 
		get(par$dataset) - fitted ~ remainder | station.id,
    data = result,
    xlab = list(label = "Residual of Spatial Loess Fit"),
    ylab = list(label = "Remainder of STL+ Fit"),
		type = "p",
		aspect = 1,
		pch = 16,
		cex = 0.2,
		layout = c(5,2),
		panel = function(x, y,...) {
      panel.abline(a=0, b=1, color="black", lty=1)
			panel.xyplot(x, y, col = "red",...)
		}
	)
	print(b)
dev.off()

trellis.device(
	postscript, 
	file = paste(local.output, "/scatterplot_of_loess.fit_", par$dataset, ".ps", sep = ""), 
	color = TRUE,
	paper = "legal"
)
  for(i in levels(result$station.id)) {
	b <- xyplot( get(par$dataset) ~ time2 | factor,
    	data = subset(result, station.id == i),
        xlab = list(label = "Month"),
        ylab = list(label = ylab),
		main = list(label = paste("Station", i, sep=" ")),
		type = "b",
		pch = 16,
		cex = 0.5,
		layout = c(1,9),
		strip = FALSE,
		xlim = c(0, 143),
		ylim = c(
			min(
				min(subset(result, station.id == i)$fitted), 
				min(with(subset(result, station.id == i), get(par$dataset)))
			), 
			max(
				max(subset(result, station.id == i)$fitted), 
				max(with(subset(result, station.id == i), get(par$dataset)))
			)
		),
		scales = list(
			y = list(relation = 'same', alternating=TRUE), 
			x = list(at=seq(0,143,by=12), relation='same')
		),
		key = list(
			text = list(
				label = c("observed","fitted")
			), 
			lines = list(
				cex = 4, 
				lwd = 1.5, 
				type = c("l","l"), 
				col = col[1:2]
			),
			columns = 2
		),		
		panel = function(x, y, subscripts,...) {
			tmp <- subset(result, station.id == i)[subscripts, ]
			panel.abline(v=seq(0,145,by=12), color="lightgrey", lty=3, lwd=0.5)
			panel.xyplot(x, y, col = col[1],...)
			panel.xyplot(x, tmp[, "fitted"], col = col[2], ...)
		}
	)
	print(b)
  }
dev.off()

trellis.device(
postscript,
file = paste(local.output, "/scatterplot_of_loess.fit.diff_", par$dataset, ".ps", sep = ""),
color = TRUE,
paper = "legal"
)
  for(i in levels(result$station.id)) {
	b <- xyplot( get(par$dataset)-fitted ~ time2 | factor,
		data = subset(result, station.id == i),
		xlab = list(label = "Month"),
		ylab = list(label = ylab),
		main = list(label = paste("Station", i, sep=" ")),
		type = "b",
		pch = 16,
		cex = 0.5,
		layout = c(1,9),
		strip = FALSE,
		xlim = c(0, 143),
		scales = list(
			y = list(relation = 'same', alternating=TRUE),
			x = list(at=seq(0,143,by=12), relation='same')
		),
		panel = function(...) {
			panel.abline(v=seq(0,145,by=12), color="lightgrey", lty=3, lwd=0.5)
			panel.xyplot(...)
		}
	)
	print(b)
  }
dev.off()

trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/QQ_plot_diff_conditional_month_", 
		par$dataset,
		".ps", 
		sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
  for(i in levels(result$station.id)){
    a <- qqmath(~ (get(par$dataset)-fitted) | month,
			data = subset(result, station.id == i),
			distribution = qnorm,
			aspect = "xy",
			layout = c(12,1),
			pch  = 16,
			cex  = 0.5,
			main = list(label = paste("Station", i, sep=" ")),
			xlab = list(label = "Unit normal quantile"),
			ylab = list(label = ylab, cex=1.2),
			prepanel = prepanel.qqmathline,
			panel = function(x, y,...) {
				panel.grid()
				panel.qqmathline(x, y=x)
				panel.qqmath(x, y, ...)
			}
		)
		print(a)
  }
dev.off()

div.norm <- function(data){
	resid <- with(data, get(par$dataset)) - data$fitted
	yy <- quantile(resid, c(0.25, 0.75))
	xx <- qnorm(c(0.25, 0.75))
	r <- diff(yy)/diff(xx)
	x <- qnorm(ppoints(length(resid)))
	y <- r*x + yy[1] - xx[1]*r
	div <- sum(abs(sort(resid) - y))
}
order.norm <- ddply(
	.data = result,
	.variables = c("station.id"),
	.fun = div.norm
)
order.norm <- as.character(order.norm[order(order.norm$V1), ]$station.id)
trellis.device(
	device = postscript, 
	file = paste(
		local.output, 
		"/QQ_plot_diff_", 
		par$dataset, 
		".ps", 
		sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
  a <- qqmath(~ (get(par$dataset)-fitted) | factor(station.id, levels = order.norm),
    data = result,
    distribution = qnorm,
    aspect = "xy",
    layout = c(10,1),
    pch  = 16,
    cex  = 0.5,
		xlab = list(label = "Unit normal quantile"),
		ylab = list(label = ylab, cex=1.2),
    prepanel = prepanel.qqmathline,
    panel = function(x, y,...) {
      panel.grid()
      panel.qqmathline(x, y=x)
      panel.qqmath(x, y, ...)
    }
  )
  print(a)
dev.off()


dd <- ddply(
	.data = result,
	.variables = c("station.id", "month"),
	.fun = summarise,
	mean = mean(fitted),
	lon = lon[1],
	lat = lat[1],
	elev = elev[1]
)
mm <- dd[rep(row.names(dd), each = 103), ]
result <- result[with(result, order(station.id, month, year)), ]
result$central <- result$fitted - mm$mean
trellis.device(
	device = postscript, 
	file = paste(local.output, "/", par$dataset, "loess.fit_conditional_month.ps", sep = ""), 
	color = TRUE, 
	paper = "legal"
)
  for(i in levels(result$station.id)) {
    b <- xyplot( central ~ year | month,
		data = subset(result, station.id == i),
		xlab = list(label = "Year"),
		ylab = list(label = ylab),
		main = list(label = paste("Station", i, sep=" ")),
		type = "p",
		pch = 16,
		cex = 0.5,
		layout = c(4, 3),
		strip = TRUE,
		key = list(
			text = list(
				label = c("centralized fitted","loess smoothing")
			), 
			lines = list(
				pch = c(".", ""), 
				cex = 4, 
				lwd = 1.5, 
				type = c("p","l"), 
				col = col[1:2]
			),
			columns = 2
		),		
		scales = list(
			y = list(relation = 'same', alternating = TRUE), 
			x = list(tick.number = 10, relation = 'same')
		),
		panel = function(x, y, ...) {
			panel.xyplot(x, y, ...)
			panel.loess(x, y, span = 2/3, degree = 1, col = col[2],...)
		}
	)
	print(b)
  }
dev.off()

trellis.device(
	device = postscript, 
	file = paste(local.output, "/", par$dataset, "loess.fit_diff_conditional_month.ps", sep = ""), 
	color = TRUE, 
	paper = "legal"
)
  for(i in levels(result$station.id)) {
    b <- xyplot( get(par$dataset)-fitted ~ year | month,
		data = subset(result, station.id == i),
		xlab = list(label = "Year"),
		ylab = list(label = ylab),
		main = list(label = paste("Station", i, sep=" ")),
		type = "p",
		pch = 16,
		cex = 0.5,
		layout = c(4, 3),
		strip = TRUE,
		key = list(
			text = list(
				label = c("residual","loess smoothing")
			), 
			lines = list(
				pch = c(".", ""), 
				cex = 4, 
				lwd = 1.5, 
				type = c("p","l"), 
				col = col[1:2]
			),
			columns = 2
		),
		scales = list(
			y = list(relation = 'same', alternating = TRUE), 
			x = list(tick.number = 10, relation = 'same')
		),
		panel = function(x, y, ...) {
			panel.xyplot(x, y, ...)
			panel.loess(x, y, span = 2/3, degree = 1, col = col[2],...)
		}
	)
	print(b)
  }
dev.off()


##########################################################
##monthly mean and std for the residual from spatial loess
##########################################################
resid.mean <- ddply(
	.data = result,
	.variables = c("station.id", "month"),
	.fun = summarise,
	mean = mean(get(par$dataset) - fitted),
	std = sd(get(par$dataset) - fitted)
)
resid.std <- ddply(
	.data = result,
	.variables = "station.id",
	.fun = summarise,
	std = sd(get(par$dataset) - fitted)
)
od <- ddply(
	.data = resid.mean,
	.variables = "station.id",
	.fun = summarise,
	overall = mean(mean)
)
resid.mean$station.id <- factor(
	resid.mean$station.id, 
	levels = od[with(od, order(overall)), ]$station.id
)
trellis.device(
	device = postscript, 
	file = paste(
		local.output, "/", par$dataset, 
		"_loess.resid_month.std_vs_month.ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
xyplot( std ~ month | station.id,
	data = resid.mean,
	layout = c(10,1),
	aspect = "xy",
	xlab = list(label = "Month", cex = 1.2),
  ylab = list(label = "Standard Diviation of Residual", cex = 1.2),	
  scales = list(
		y = list(relation = 'same', alternating = TRUE), 
		x = list(at=c(1, 7, 12), relation = 'same')
  ),
  panel = function(x, y, subscripts, ...){
		panel.xyplot(x,y,col=col[1],type="b",...)
	}
)
dev.off()

trellis.device(
	device = postscript, 
	file = paste(
		local.output, "/", par$dataset, 
		"_loess.resid_month.mean_vs_month.ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
xyplot( mean ~ month | station.id,
	data = resid.mean,
	layout = c(5,1),
	col = col[1],
	type = "b",
	cex = 0.8,
	xlab = list(label = "Month", cex = 1.2),
  ylab = list(label = "Mean of Residual", cex = 1.2),	
  scales = list(
		y = list(relation = 'same', alternating = TRUE), 
		x = list(at=c(1, 7, 12), relation = 'same')
  )
)
dev.off()

trellis.device(
	device = postscript, 
	file = paste(
		local.output, "/", par$dataset, 
		"_dotplot_overall_mean_of_loess.resid.ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
qqmath( ~ overall,
	data = od,
	distribution = qunif,
	xlab = "f-value",
	ylab = "Mean of Residuals"
)
dev.off()

trellis.device(
	device = postscript, 
	file = paste(
		local.output, "/", par$dataset, 
		"_dotplot_overall_std_of_loess.resid.ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
qqmath( ~ std,
	data = resid.std,
	distribution = qunif,
	xlab = "f-value",
	ylab = "SD of Residuals"
)
dev.off()

##########################################################
##Auto correlation ACF for the residual from spatial loess
##########################################################
ACF <- ddply(
	.data = result,
	.variables = "station.id",
	.fun = summarise,
	correlation = c(acf(get(par$dataset)-fitted, plot=FALSE)$acf),
	lag = c(acf(get(par$dataset)-fitted, plot=FALSE)$lag) 
)
ACF$station.id <- factor(
	ACF$station.id, 
	levels = od[with(od, order(overall)), ]$station.id
)
trellis.device(
	device = postscript, 
	file = paste(
		local.output, "/acf_of_", par$dataset, 
		"_loess.fit_remainder_for_100_stations",
		".ps", sep = ""
	), 
	color = TRUE, 
	paper = "legal"
)
for(i in levels(ACF$station.id)){
	b <- xyplot( correlation ~ lag,
		data = subset(ACF, station.id == i & lag != 0),
		xlab = list(label = "Lag", cex = 1.2),
    ylab = list(label = paste("Station", i, "ACF"), cex = 1.2),
    type = "h",
    panel = function(x,y,...) {
      panel.abline(h=0)
      panel.xyplot(x,y,...)
    }
  )
  print(b)
}
dev.off()