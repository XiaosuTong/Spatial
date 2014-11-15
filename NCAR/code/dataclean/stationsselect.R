#load the data
library(lattice)
library(plyr)
datadir <- "~/Projects/Spatial/NCAR/RData/"
outputdir <- "~/Projects/Spatial/NCAR/output/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))

TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)


#find the 100 stations with largest observation number
Pcount <- sort(Pcount, decreasing=TRUE)
stations.precip <- head(names(Pcount),100)

TMincount <- sort(TMincount, decreasing=TRUE)
stations.tmin <- head(names(TMincount),100)

TMaxcount <- sort(TMaxcount, decreasing=TRUE)
stations.tmax <- head(names(TMaxcount),100)


##find the 100 stations with largest observation number and close with each other
stations.tmax.close <- names(
	TMaxcount[which(TMaxcount >= 1236)]
)
stations.tmax.close <- as.character(
	subset(
		x = UStinfo, 
		subset = station.id %in% stations.tmax.close & lon>=-85.5 & lon<=-70 & lat>=38 & lat<= 43
	)$station.id
)[1:100]

stations.tmin.close <- names(
	TMincount[which(TMincount >= 1236)]
)
stations.tmin.close <- as.character(
	subset(
		x = UStinfo, 
		subset = station.id %in% stations.tmin.close & lon>=-100 & lat <= 41 & lat >= 29
	)$station.id
)[1:100]

stations.precip.close <- names(
	Pcount[which(Pcount >= 1236)]
)
stations.precip.close <- as.character(
	subset(
		x = USpinfo, 
		subset = station.id %in% stations.precip.close & lon>=-100 & lat>=38.5 &lat<=42.5
	)$station.id
)[1:100]


save(
	list=c(
		"stations.tmax", "stations.tmin", "stations.precip", 
		"stations.tmax.close", "stations.tmin.close", "stations.precip.close"
	), 
	file=paste(datadir, "stations.100.RData", sep="")
)

## find the stations at list have one obs after given year for given variable
find.stations <- function(dataset, variable, y = 1950, na = TRUE) {
	tmp <- subset(dataset, year >= y)
	if(na){
		tmp <- tmp[!is.na(tmp[[variable]]),]
	}
	tmp$station.id <- as.character(tmp$station.id)
	stations <- unique(tmp$station.id)
}

stations.a1950.tmax <- find.stations(UStemp, "tmax")
stations.a1950.tmin <- find.stations(UStemp, "tmin")
save(
	list = c(
		"stations.a1950.tmax", "stations.a1950.tmin"
	),
	file = paste(datadir, "stations.a1950.RData", sep="")
)