load(file.path(local.datadir,"USmonthlyMet.RData"))
load(file.path(local.datadir, "stations.RData"))

#save the stations information on HDFS
rhsave(
    list = c("USpinfo", "UStinfo"), 
    file = file.path(rh.datadir, "stationinfo", "USinfo.RData")
)

#Count the observation number for each station.
TMaxcount <- tapply(!is.na(UStemp$tmax), UStemp$station.id, sum)

#Order the stations by the observation counts
TMaxcount <- sort(TMaxcount, decreasing=TRUE)
tmp.tmax <- UStemp
tmp.tmax$station.id <- factor(
    tmp.tmax$station.id, 
    levels = names(TMaxcount)
)
tmp.tmax <- tmp.tmax[order(tmp.tmax$station.id),]

#The reason for using dlply instead of using lapply is because in lapply, every time when calling subset() function,
#it will iterative over all observations. However in dlply, the subset will be done once when we specify the .variables
tmp <- dlply(
    .data = tmp.tmax,
    .variables = "station.id",
    .fun = function(r) {
        tmaxSub <- r[,c(6,7,9)]
        attributes(tmaxSub)$lat  <- unique(r$lat)
        attributes(tmaxSub)$elev  <- unique(r$elev)
        attributes(tmaxSub)$lon <- unique(r$lon)
        attributes(tmaxSub)$station.name <- as.character((unique(r$station.name)))
        list(as.character(r$station.id)[1], tmaxSub)
    }
)

#Write the data which is a list with # of stations of elements, each element is a list of key-value pair. 
#Key is the station.id, value is dataframe for that station.
rhwrite(tmp, file.path(rh.datadir,"tmax","data"))
tmax <- tmp.tmax[,-8]
rhsave(
    list = ("tmax"), 
    file = file.path(rh.datadir, "tmax","Rdata","tmax.RData")
)
tmax.100stations <- subset(tmax, station.id %in% stations.tmax)
rhsave(
    list = ("tmax.100stations"), 
    file = file.path(rh.datadir, "tmax","100stations","Rdata","tmax.100stations.RData")
)
rm("tmax", "tmp", "tmax.100stations")

#Count the observation number for each station.
TMincount <- tapply(!is.na(UStemp$tmin), UStemp$station.id, sum)

#Order the stations by the observation counts
TMincount <- sort(TMincount, decreasing=TRUE)
tmp.tmin <- UStemp
tmp.tmin$station.id <- factor(
    tmp.tmin$station.id, 
    levels = names(TMincount)
)
tmp.tmin <- tmp.tmin[order(tmp.tmin$station.id),]

##The reason for using dlply instead of using lapply is because in lapply, every time when calling subset() function,
##it will iterative over all observations. However in dlply, the subset will be done once when we specify the .variables
tmp <- dlply(
    .data = tmp.tmin,
    .variables = "station.id",
    .fun = function(r) {
        tminSub <- r[,c(6,7,8)]
        attributes(tminSub)$lat  <- unique(r$lat)
        attributes(tminSub)$elev  <- unique(r$elev)
        attributes(tminSub)$lon <- unique(r$lon)
        attributes(tminSub)$station.name <- as.character(unique(r$station.name))
        list(as.character(r$station.id)[1], tminSub)
    }
)

#Write the data which is a list with # of stations of elements, each element is a list of key-value pair. 
#Key is the station.id, value is the dataframe of that station.
rhwrite(tmp, file.path(rh.datadir,"tmin","data"))
tmin <- tmp.tmin[,-9]
rhsave(
   list = ("tmin"), 
   file = file.path(rh.datadir, "tmin","Rdata","tmin.RData")
)
tmin.100stations <- subset(tmin, station.id %in% stations.tmin)
rhsave(
    list = ("tmin.100stations"), 
    file = file.path(rh.datadir, "tmin","100stations","Rdata","tmin.100stations.RData")
)
rm("tmin", "tmp", "tmin.100stations")


#Count the observation number for each station.
Pcount <- tapply(!is.na(USppt$precip), USppt$station.id, sum)

#Order the stations by the observation counts
Pcount <- sort(Pcount, decreasing=TRUE)
tmp.precip <- USppt
tmp.precip$station.id <- factor(
    tmp.precip$station.id, 
    levels = names(Pcount)
)
tmp.precip <- tmp.precip[order(tmp.precip$station.id),]


##The reason for using dlply instead of using lapply is because in lapply, every time when calling subset() function,
##it will iterative over all observations. However in dlply, the subset will be done once when we specify the .variables
tmp <- dlply(
    .data = tmp.precip,
    .variables = "station.id",
    .fun = function(r) {
        precipSub <- r[,c(6,7,8)]
        attributes(precipSub)$lat  <- unique(r$lat)
        attributes(precipSub)$elev  <- unique(r$elev)
        attributes(precipSub)$lon <- unique(r$lon)
        attributes(precipSub)$station.name <- as.character(unique(r$station.name))
        list(as.character(r$station.id)[1], precipSub)
    }
)


#Write the data which is a list with # of stations of elements, each element is a list of key-value pair. 
#Key is the station.id, value is the dataframe of that station.
rhwrite(tmp, file.path(rh.datadir, "precip","data"))
precip <- tmp.precip
rhsave(
    list = ("precip"), 
    file = file.path(rh.datadir, "precip","Rdata","precip.RData")
)
precip.100stations <- subset(precip, station.id %in% stations.precip)
rhsave(
    list = ("precip.100stations"), 
    file = file.path(rh.datadir, "precip","100stations","Rdata","precip.100stations.RData")
)
rm("precip", "tmp", "precip.100stations")

#################################
##Create the RData for after1950
#################################
load(file.path(local.datadir,"USmonthlyMet.RData"))
load(file.path(local.datadir, "stations.a1950.RData"))
tmax.a1950 <- subset(UStemp, station.id %in% stations.a1950.tmax & year >= 1950)[, -8]
rhsave(
  list = ("tmax.a1950"),
  file = file.path(rh.datadir, "tmax", "a1950", "Rdata", "tmax.a1950.RData")
)