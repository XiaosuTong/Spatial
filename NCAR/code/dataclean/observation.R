###############################################################
##Count the observation number in precipitation and temperature
##############################################################

#Set up the directory and load the data.
library(lattice)
library(plyr)
library(maps)
datadir <- "~/Projects/Spatial/NCAR/RData/"
outputdir <- "~/Projects/Spatial/NCAR/output/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))

sum(unique(USppt$station.id) %in% unique(UStemp$station.id))
##There are 4267 stations only for precipitation, 
##7651 stations for both, 
##and 474 stations only for temperature,
##totally, there are 12,392 stations.

dim(USppt[!is.na(USppt[,8]),])
dim(UStemp[!is.na(UStemp[,8]),])
##There are 6,204,442 observations for precipitation
##and 4,285,841 observations for temperature.
##Totally, there should be 11918*12*103=14,730,648 observations for precipitation
##and should be 8125*12*103=10,042,500 observations for temperature.


#Count the observation number for each station.
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)
TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)




########################################################################
##QQ plot of observation count of each station for three types of response.
######################################################################
#count the observation number for each station at least have one obs after 1950
tmp <- subset(UStemp, year>= 1950)
tmp <- tmp[!is.na(tmp$tmax),]
tmp$station.id <- as.character(tmp$station.id)
TMaxcount <- tapply(!is.na(tmp$tmax), tmp$station.id, sum)
#data <- c(Pcount, TMaxcount, TMincount)
#names(data) <- NULL

#data <- data.frame(
#    response = data, 
#    type= rep(c("Precip","tmax","tmin"), c(length(Pcount), length(TMaxcount), length(TMincount)))
#)
data <- data.frame(as.vector(TMaxcount))
names(data) <- c("response")
trellis.device(postscript, 
    file = paste(outputdir, 
        "Distn_of_Observation_count_a1950_tmax", 
        ".ps", sep = ""
    ), 
    color=TRUE, 
    paper="legal"
)
a <- qqmath(
    ~ response,	
	data = data,
	distribution = qunif,
	grid = TRUE,
	pch = 16,
	cex = 0.3,
    main = "Maximum Temperature",
    xlab = list(label="f-value", cex=1.2),
	ylab = list(label="Number of Observations", cex=1.2),
    panel = function(...){
        panel.qqmath(...)
        panel.abline(h = 576, col = "red", lty=2, lwd = 1)
    }
)
print(a)
dev.off()

trellis.device(
    device = postscript, 
    file = paste(outputdir, "Rate_of_Observation_count", ".ps", sep = ""), 
    color=TRUE,
    paper="legal"
)
a <- qqmath(
    ~ response/1236 | type,
    data = data,
    distribution = qunif,
    grid = TRUE,
    pch = 16,
    cex = 0.5,
    xlab = list(label="f-value", cex=1.5),
    ylab = list(label="Rate of Valid Observation", cex=1.5)
)
print(a)
dev.off()

Pcount <- ddply(
    .data = USppt[!is.na(USppt$precip),],
    .variable = "station.id",
    .fun = summarize,
    zero = sum(precip == 0),
    total = length(precip)
)
trellis.device(postscript, 
    file = paste(outputdir, 
        "Distn_of_zero_count_precip", 
        ".ps", sep = ""
    ), 
    color=TRUE, 
    paper="legal"
)
a <- qqmath(
    ~ zero/1236, 
    data = Pcount,
    distribution = qunif,
    grid = TRUE,
    pch = 16,
    cex = 0.4,
    main = "Precipitation",
    xlab = list(label="f-value", cex=1.2),
    ylab = list(label="Rate of Zero over All", cex=1.2),
)
b <- qqmath(
    ~ zero/total, 
    data = Pcount,
    distribution = qunif,
    grid = TRUE,
    pch = 16,
    cex = 0.4,
    main = "Precipitation",
    xlab = list(label="f-value", cex=1.2),
    ylab = list(label="Rate of Zero over Valid Observations", cex=1.2),
)
print(a)
print(b)
dev.off()

################################################
##Line plot of observation count vs. time
###############################################
pcount <- ddply(
    .data = USppt,
	.variable = c("year","month"),
	.fun = summarize,
	count = sum(!is.na(precip))
)
pcount$group <- "precip"
tmincount <- ddply(
    .data = UStemp,
	.variable = c("year", "month"),
	.fun = summarize,
	count = sum(!is.na(tmin))
)
tmincount$group <- "tmin"
tmaxcount <- ddply(
    .data = UStemp,
    .variable = c("year", "month"),
    .fun = summarize,
    count = sum(!is.na(tmax))
)
tmaxcount$group <- "tmax"
#data <- do.call(rbind, list(pcount, tmaxcount, tmincount))
data <- tmaxcount
data$month <- as.numeric(data$month)
date <- paste(data$year, data$month, "01", sep="-")
data$date <- as.POSIXct(strptime(date, format = "%Y-%m-%d"), format='%Y%m%d', tz="")
data <- data[order(data$date),]
#data$total <- c(11918, 8125, 8125)
#numstation <- length(
#    unique(
#        USppt[!is.na(USppt$precip) & USppt$year >= 1955,]$station.id
#    )
#)
numstation <- length(
    unique(
        UStemp[!is.na(UStemp$tmax) & UStemp$year >= 1950,]$station.id
    )
)
#start <- as.POSIXct(
#    strptime(paste(head(data$year,1)-1, "01", "01", sep="-"), format = "%Y-%m-%d"), 
#    format='%Y%m%d', 
#    tz=""
#)
start <- as.POSIXct(
    strptime(paste("1950", "01", "01", sep="-"), format = "%Y-%m-%d"), 
    format='%Y%m%d', 
    tz=""
)
end <- as.POSIXct(
    strptime(paste(tail(data$year,1)+1, "12", "01", sep="-"), format = "%Y-%m-%d"), 
    format='%Y%m%d', 
    tz=""
)

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

trellis.device(
    postscript, 
    file = paste(
        outputdir, 
        "tmax.lineplot_of_observation_count",
        ".ps", 
        sep = ""
    ), 
    color=TRUE, 
    paper="legal"
)
##comment out part is the plot function for overall obs counts for temp and precip
##no comment out part is for the precip after year 1955.
b <- xyplot( 
    count ~ date,
	data = subset(data, year >= 1950),
	type = c("l"),
	grib = TRUE,
#	groups = group,
	xlab = list(label = "Time", cex = 1.2),
	xlim = c(start, end),
    ylim = c(4000, 8000),
	ylab = list(label="Number of Observation",cex=1.2),
    main = "Maximum Temperature",
	scales = list(
        y = list(relation = 'free'), 
        x = list(format = "%b %Y", tick.number = 14)
    ),                
#    key = list(
#        columns = 3, 
#        text = list(c("Precip", "Tmax", "Tmin"), cex = 1.2), 
#        lines=list(lwd = 3, col = col[1:3])
#    ),
	panel = function(x,y,...) {
		panel.abline(
#            h=seq(0,13,by=1), 
            v=seq(start, end, by="60 month"), 
            color="lightgrey", 
            lty=3, 
            lwd=0.5
        )
        panel.abline(h = 7738, col = "red", lty=1, lwd = 1)
		panel.xyplot(x,y,...)
	}
)
print(b)
dev.off()

trellis.device(
  device = postscript, 
  file = paste(outputdir, "na_rate_of_observation_count", ".ps", sep = ""), 
  color = TRUE, 
  paper = "legal"
)
b <- xyplot( 
    (1 - count/total) ~ date,
    data = data,
    type = c("l"),
    grib = TRUE,
    groups = group,
    xlab = list(label = "Time", cex = 1.5),
    xlim = c(start, end),
    ylab = list(label="Rate of Missing Value",cex=1.5),
    scales = list(
        y = list(relation = 'free', cex = 1.5), 
        x = list(format = "%b %Y", tick.number = 14), 
        cex = 1.2
    ),
    key = list(
        columns = 3, 
        text = list(c("Precip", "Tmax", "Tmin"), cex = 1.2), 
        lines=list(lwd = 3, col = col[1:3])
    ),
    panel = function(x,y,...) {
        panel.abline(
            h=seq(0,13,by=1), 
            v=seq(start, end, by="48 month"), 
            color="lightgrey", 
            lty=3, 
            lwd=0.5
        )
        panel.xyplot(x,y,...)
        panel.text(x[which.max(y)], max(y), round(max(y),2), adj=c(0,0))
		panel.text(x[which.min(y)], min(y), round(min(y),2), adj=c(1,1))
    }
)
print(b)
dev.off()

###################################################################
##Quantile plot for max length of continuous missing value for tmax
###################################################################
data <- subset(UStemp, year >= 1950)
data <- data[order(data$station.id,data$year, data$month),]

## seqle function calculate the length of consecutive integer of a vector.
## the input to this function for our example, is the index of NA in tmax.
seqle <- function(x,incr=1) { 
  if(!is.numeric(x)) x <- as.numeric(x) 
  n <- length(x)  
  y <- x[-1L] != x[-n] + incr 
  i <- c(which(y|is.na(y)),n) 
  list(lengths = diff(c(0L,i)),
       values = x[head(c(0L,i)+1L,-1L)]) 
} 

missing <- ddply(
  .data = data,
  .variable = "station.id",
  .fun = function(r) {
    tmp <- r[order(r$year, r$month),]
    me <- is.na(tmp$tmax) ## what is the index for NA
    max(seqle(which(me))$lengths)
  }
)
names(missing)[2] <- "consec.miss"
missing <- subset(missing, consec.miss != 572)
trellis.device(postscript, 
    file = paste(outputdir, 
        "Distn_of_max.consecutive.missing_a1950_tmax", 
        ".ps", sep = ""
    ), 
    color=TRUE, 
    paper="legal"
)
a <- qqmath(
  ~ consec.miss, 
  data = subset(missing, consec.miss != 576),
  distribution = qunif,
  grid = TRUE,
  pch = 16,
  cex = 0.3,
  main = "Maximum Temperature",
  xlab = list(
    label="f-value", 
    cex=1.2
  ),
  ylab = list(
    label = "Maximum Length of Consecutive Missing Value", 
    cex = 1.2
  ),
  panel = function(x,...){
    panel.qqmath(x,...)
    panel.text(1, max(x), max(x), adj=c(0,0))
  }
)
print(a)
dev.off()
