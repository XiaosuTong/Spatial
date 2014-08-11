###############################################
##Plot the location of stations in a US map
###############################################

#load the data
library(lattice)
library(maps)
library(plyr)
datadir <- "~/Projects/Spatial/NCAR/RData/"
outputdir <- "~/Projects/Spatial/NCAR/output/"

load(paste(datadir,"USmonthlyMet.RData", sep=""))
load(file.path(datadir, "stations.RData"))

##Find out all unique stations
temp <- UStinfo
precip <- USpinfo
tmp <- temp[!(temp$station.id %in% precip$station.id),]
data <- rbind(precip, tmp)
data <- data[order(data$elev),]
data$group <- rep(paste("level", 1:8, sep=" "), each = 1549)
data$group <- factor(data$group, levels=unique(data$group))
rm(temp)
rm(precip)

#index <- c(seq(-100,3700,by=100),3810)
#data$group <- cut(data$elev, index)

us.map <- map('state', plot = FALSE, fill = TRUE)

label <- ddply(
    .data = data,
    .variables = "group",
    .fun = summarise,
    label = paste(min(elev), "m", " ~ ", max(elev), "m", sep="")
)
levels(data$group) <- label$label

##############################
##Stations location on US map
##############################
trellis.device(postscript, file = paste(outputdir, "spatial", ".ps", sep = ""), color=TRUE, paper="legal")
b <- xyplot(
    lat ~ lon | group,
	data = data,
	xlab = list(label="Longitude", cex=1.5),
	ylab = list(label="Latitude", cex=1.5),
	layout = c(1,1),
    pch = 16,
    cex = 0.25,
    col = "red",
	strip = strip.custom(par.strip.text= list(cex = 1.5)),
	par.settings = list(layout.heights = list(strip = 1.5)),
    panel = function(...) {
        panel.polygon(us.map$x,us.map$y)   
        panel.xyplot(...)
    }
)
print(b)
dev.off()

#############################################################
##dot plot of max/min elevation for each level in spatial.ps
#############################################################

dr.spatial <- ddply(
    .data = data,
	.variables = "group",
	.fun = summarise,
	type = c("min","max"),
	elev = range(elev),
	label = rep(paste(min(elev),max(elev), sep="~"),2)
)

dr.spatial$label <- factor(dr.spatial$label, levels=unique(dr.spatial$label))

lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

trellis.device(postscript, file = paste(outputdir,"elevation_range", ".ps", sep = ""), color=TRUE, paper="legal")
b <- dotplot(
    label ~ elev,
    data = dr.spatial,
    xlab = list(label="Elevation (meter)", cex=1.5),
    ylab = list(label="Level", cex=1.5),
    groups = type,
	grib = TRUE,
	key = list(columns=2, points= list(pch=16, col = col[1:2]), text= list(label=c("Max","Min"), cex=1.2)),
	scales = list(x = list(cex=1.5), y = list(cex=1.5)),          
	pch = 16,
    panel = function(...) {
        panel.dotplot(...)
    }
)
print(b)
dev.off()



#######################################
##QQ plot of elevation of all stations
#######################################
trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_elevation", ".ps", sep = ""), color=TRUE, paper="legal")
a <- qqmath(
    ~ elev,
    data = data,
    distribution = qunif,
    aspect = 1,
	pch = 16,
	cex = 0.3,
	grib = FALSE,
    xlab = list(label="f-value", cex=1.5),
    ylab = list(label="Elevation (meter)", cex=1.5),
    scales = list(x = list(cex=1.5), y = list(cex=1.5)),
    prepanel = prepanel.qqmathline,
    panel = function(x, y,...) {
		panel.abline( 
            h = dr.spatial$elev, 
            color = "lightgrey", 
            lty = 3, 
            lwd = 0.5
        )
        panel.qqmath(x, y,...)
    }
)
print(a)
dev.off()

######################################################
##QQ plot of elevation of stations whoes above 1 meter
######################################################
trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_elevation_log", ".ps", sep = ""), color=TRUE, paper="legal")
a <- qqmath(
    ~ log2(elev),
    data = subset(data, elev>=1),
    distribution = qunif,
    aspect = 1,
    pch = 16,
    cex = 0.3,
    grib = FALSE,
    xlab = list(label="f-value", cex=1.5),
    ylab = list(label="Log of Elevation (log base 2 meter)", cex=1.5),
    scales = list(
        x = list(cex=1.5), 
        y = list(cex=1.5)
    ),
    prepanel = prepanel.qqmathline,
    panel = function(x, y,...) {
        panel.abline( 
            h= dr.spatial$elev, 
            color="lightgrey", 
            lty=3, 
            lwd=0.5
        )
        panel.qqmath(x, y,...)
    }
)
print(a)
dev.off()


#####################################################################
##QQ plot of elevation of stations conditional on level of elevation
#####################################################################
trellis.device(postscript, file = paste(outputdir, "QQ_plot_of_elevation_conditional_on_level", ".ps", sep = ""), color=TRUE, paper="legal")
a <- qqmath(
    ~ elev | group,
    data = data,
    distribution = qunif,
    aspect = 1,
	pch = 16,
	cex = 0.3,
	grib = FALSE, 
	layout = c(4,2),
    xlab = list(label="f-value", cex=1.5),
    ylab = list(label="Elevation (meter)", cex=1.5),
    scales = list(x = list(cex=1.5), y = list(cex=1.5)),
    prepanel = prepanel.qqmathline,
    panel = function(x, y,...) {
        panel.qqmath(x, y,...)
    }
)
print(a)
dev.off()



#####################################################
##Location of 100 stations in US map
######################################################

us.map <- map('state', plot = FALSE, fill = TRUE)
tmp1 <- UStinfo[with(UStinfo, station.id %in% stations.tmax), ]
tmp1$group <- "max temperature"
tmp2 <- UStinfo[with(UStinfo, station.id %in% stations.tmin), ]
tmp2$group <- "min temperature"
tmp3 <- USpinfo[with(USpinfo, station.id %in% stations.precip), ]
tmp3$group <- "precipitation"
data <- rbind(tmp1, tmp2, tmp3)

trellis.device(postscript, file = paste(outputdir, "spatial_100_stations", ".ps", sep = ""), color=TRUE, paper="legal")
b <- xyplot(
    lat ~ lon | group,
    data = data,
    xlab = list(label="Longitude", cex=1.5),
    ylab = list(label="Latitude", cex=1.5),
    layout = c(1,1),
    pch = 16,
    cex = 0.3,
    col = "red",
    strip = strip.custom(par.strip.text= list(cex = 1.5)),
    par.settings = list(layout.heights = list(strip = 1.5)),
    panel = function(x,y,...) {
        panel.polygon(us.map$x,us.map$y)
		for(i in 1:length(x)){
			panel.text(
                x = x[i], 
                y = y[i],
                labels = data[which(data$lon==x[i]&data$lat==y[i]),]$station.id, 
                col = "blue", 
                adj = c(0,0), 
                cex=0.5
            )
		}
        panel.xyplot(x,y,...)
    }
)
print(b)
dev.off()


#######################################################################################################
##Or find the  stations which have full observation # for both temperature and precipitation
########################################################################################################
tmp <- UStemp[(UStemp[1:8125,]$station.id %in% USppt[1:11918,]$station.id),]

station <- as.character(unique(tmp$station.id))

count1 <- ddply(
    .data = subset(USppt, station.id %in% station),
    .variables = "station.id",
    .fun = summarise,
     count.p = sum(!is.na(precip))
)
#full <- subset(count1, count.p == 1236)
#there are 287 stations have full observations for precipitation

count2 <- ddply(
    .data = subset(UStemp, station.id %in% station),
    .variables = "station.id",
    .fun = summarise,
    count.t = sum(!is.na(tmax)|!is.na(tmin))
)
#full <- subset(count2, count.p == 1236)
#there are 513 stations have full observations for temperature
count.pt <- merge(count1,count2)
count.pt <- count.pt[order(count.pt$count.p, count.pt$count.t, decreasing=TRUE),]
count.pt[1:100,]

data.tmp <- USpinfo[USpinfo$station.id %in% count.pt[1:100,]$station.id,]
##################################################