source("~/Rhipe/rhinitial.R")
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
par <- list()
par$dataset <- "tmax"
par$loess <- "loess01"

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
    .variables = "fac",
    .fun = summarise,
    mean = mean(fitted)
)
order.st <- as.character(
	order[order(order$mean, decreasing=TRUE), ]$fac
)
result$fac <- factor(
	result$fac, 
	levels = order.st
)
