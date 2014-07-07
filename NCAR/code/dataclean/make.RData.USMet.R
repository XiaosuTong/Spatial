################################################
##Write the data into txt file and RData object
################################################

datadir<- "~/Projects/Spatial/NCAR/Raw"

######################################Temperature#################################################
#Scan META information of all stations into list of three vectors.
#In scan() function, "what" statement is used to specify the type of data to be read.
#Scan() can also be used to scan data from standard input.
#look is a list object with length=4, which is the spatial information of each station. 
look <- scan(paste(datadir, "/NCAR_tinfill/METAinfo", sep=""), skip=1, what=list( "a", 1,1,1))
names(look) <- c("station.id", "elev", "lon", "lat")

UStinfo <- look

#Read a table of fixed width formatted data into a data.frame.
#look2 is a dataframe that is the station name of each station.
look2 <- read.fwf(paste(datadir, "/USmonthly.names.txt", sep=""), widths=c(2,6,1,12))[,c(2,4)]

#Match the stations spatial information and stations name.
#Find the position of station.id of look$station.id in look2(station.id). look3 is the postions in look2[,1].
look3 <- match(look$station.id, look2[,1])

#Assign the station name to the each station. 
UStinfo$station.name <- as.character(look2[,2])[look3]
UStinfo <- as.data.frame(UStinfo)

N <- length(look$station.id)


# file names for tmax and tmin saved to be two vectors named "fnames"
cmd<-paste( "ls ",datadir,  "/NCAR_tinfill/tmax*", sep="")
fnames.tmax<- system( cmd, intern=TRUE)
cmd<-paste( "ls ",datadir,  "/NCAR_tinfill/tmin*", sep="")
fnames.tmin<- system( cmd, intern=TRUE)

# now accumulate data

# loop over data
# initial the dataframe.
UStmax <- UStmin <- data.frame()

#UStmax and UStmin are two data.frame each of which has 13 columns, # of stations*103 rows. Sorted by year. 
for ( k in 1:103){
	print( k)
	# read in funny formated file for year k
	dat <- read.fwf( fnames.tmax[k], width = c(6, 7, rep(5,11), 3, rep(1,11)))
	#miss indicates which observation is missing. Even though all time point have observation measurement, if miss == 1, it is actually missing,
	#filled by the prediction value.
	miss <- dat[, (1:12) + 13]
	tmp <- dat[, (1:12) +  1]
	name <- dat[, 1]
	tmp <- tmp/10
	tmp[miss == 1] <- NA
	names(tmp) <- c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec")
	tmp <- cbind(station.id=name, tmp, year=rep((k + 1894)))
	UStmax <- rbind(UStmax, tmp)

	dat <- read.fwf( fnames.tmin[k], width = c(6, 7, rep(5,11), 3, rep(1,11)))
	miss <- dat[,(1:12) +13]
	tmp <- dat[,(1:12) +1 ]
	name <- dat[, 1]
	tmp <- tmp/10
	tmp[miss==1]<- NA
	names(tmp) <- c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec")
	tmp <- cbind(station.id=name, tmp, year=rep((k + 1894)))
	UStmin <- rbind(UStmin, tmp)
}

#UStemp sorted by month, year, station.id
UStemp <- data.frame(station.id = rep(UStmin$station.id, 12), 
		     elev = rep(UStinfo$elev, 12*103),
		     lon = rep(UStinfo$lon, 12*103),
		     lat = rep(UStinfo$lat, 12*103),
		     station.name = rep(UStinfo$station.name, 12*103),
		     year=rep(UStmin$year,12), 
		     month=rep(c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"), each=dim(UStmin)[1]), 
		     tmin=c(UStmin[,2], UStmin[,3], UStmin[,4], UStmin[,5], UStmin[,6], UStmin[,7], UStmin[,8], UStmin[,9], UStmin[,10],UStmin[,11], UStmin[,12], UStmin[,13]),
		     tmax=c(UStmax[,2], UStmax[,3], UStmax[,4], UStmax[,5], UStmax[,6], UStmax[,7], UStmax[,8], UStmax[,9], UStmax[,10],UStmax[,11], UStmax[,12], UStmax[,13]))

write.table(UStemp, file="UStemp.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE)



#############################Precipitation#######################################################

#look <- scan(paste(datadir,"/NCAR_pinfill/METAinfo", sep="") ,skip=1, what=list( "a", 1,1,1))
look <- read.table( paste(datadir,"/NCAR_pinfill/METAinfo",sep=""), header=TRUE, row.names=NULL)
look <- look[, c( 1,4,2,3)]
names(look) <- c("station.id", "elev", "lon", "lat")

USpinfo <- look

#Read a table of fixed width formatted data into a data.frame.
#look2 is a dataframe that is the station name of each station.
look2 <-read.fwf(paste(datadir, "/USmonthly.names.txt", sep=""), widths=c(2,6,1,12))[,c(2,4)]

#Match the stations spatial information and stations name.
#Find the position of station.id of look$station.id in look2(station.id). look3 is the postions in look2[,1].
look3<- match(USpinfo$station.id, look2[,1])

USpinfo$station.name<- as.character(look2[,2])[look3]
USpinfo <- as.data.frame(USpinfo)

N<- length(look$station.id)

# file names for ppt

cmd <- paste( "ls ",datadir,  "/NCAR_pinfill/ppt*", sep="")
fnames.ppt <- system( cmd, intern=TRUE)

### now accumulate data

# loop over data

USppt <- data.frame()
 
for ( k in 1:103){
	print(k)
	# read in funny formated file for year k
	dat <- read.fwf( fnames.ppt[k], width= c(6, 7,rep(5,11), 3,rep(1,11)))
	miss <- dat[,(1:12) +13]
	tmp <- dat[,(1:12) +1 ]
	name <- dat[,1]
	tmp <-  tmp/10
	tmp[miss==1] <- NA
	names(tmp) <- c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec")
        tmp <- cbind(station.id=name, tmp, year=rep((k + 1894)))
        USppt <- rbind(USppt, tmp)
}


USppt <- data.frame(station.id = rep(USppt$station.id, 12),
                     elev = rep(USpinfo$elev, 12*103),
                     lon = rep(USpinfo$lon, 12*103),
                     lat = rep(USpinfo$lat, 12*103),
                     station.name = rep(USpinfo$station.name, 12*103),
                     year = rep(USppt$year,12),
                     month = rep(c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"), each=dim(USppt)[1]),
                     precip = c(USppt[,2], USppt[,3], USppt[,4], USppt[,5], USppt[,6], USppt[,7], USppt[,8], USppt[,9], USppt[,10],USppt[,11], USppt[,12], USppt[,13]))


write.table(USppt, file="USppt.txt", quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE, col.names = FALSE)

UStemp$month <- factor(UStemp$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))
USppt$month <- factor(USppt$month, levels=c("Jan", "Feb", "Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"))

save( list=c("UStemp","UStinfo", "USppt", "USpinfo"), file= "~/Projects/Spatial/NCAR/RData/USmonthlyMet.RData")

print( "all done!")

