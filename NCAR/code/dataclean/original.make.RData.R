datadir<- "~/Projects/Spatial"

######## temperature
scan(paste(datadir,"/NCAR_tinfill/METAinfo", sep="") ,
skip=1, what=list( "a", 1,1,1))-> look
names( look)<-c("station.id", "elev", "lon", "lat")

UStinfo<- look
look2<-read.fwf("USmonthly.names.txt",
                     widths=c(2,6,1,12))[,c(2,4)]
look3<- match(look$station.id, look2[,1])

UStinfo$station.name<- as.character(look2[,2])[look3]
UStinfo <- as.data.frame( UStinfo)

N<- length( look$station.id)


# file names for tmax and tmin

cmd<-paste( "ls ",datadir,  "/NCAR_tinfill/tmax*", sep="")
fnames.tmax<- system( cmd, intern=TRUE)
cmd<-paste( "ls ",datadir,  "/NCAR_tinfill/tmin*", sep="")
fnames.tmin<- system( cmd, intern=TRUE)

# now accumulate data
#

# loop over data

UStmax<- UStmin<- array(NA, c( 103, 12,N) )
 
for ( k in 1:103){
print( k)
# read in funny formated file for year k
read.fwf( fnames.tmax[k],
       width= c(6,   7,rep(5,11), 3,rep(1,11))) -> dat
 miss<- as.matrix( dat[,(1:12) +13])
 dat<-  as.matrix(  dat[,(1:12) +1 ])
 dat<-  matrix( as.single( dat/10), ncol=12)
 dat[miss==1]<- NA
UStmax[k,,] <- dat

read.fwf( fnames.tmin[k],
       width= c(6,   7,rep(5,11), 3,rep(1,11))) -> dat
 miss<- as.matrix( dat[,(1:12) +13])
 dat<-  as.matrix(  dat[,(1:12) +1 ])
 dat<-  matrix( as.single( dat/10), ncol=12)
 dat[miss==1]<- NA
UStmin[k,,] <- dat

  }

############################ precip
scan(paste(datadir,"/NCAR_pinfill/METAinfo", sep="") ,
  skip=1, what=list( "a", 1,1,1))-> look

read.table( paste(datadir,"/NCAR_pinfill/METAinfo",sep=""), header=TRUE,
row.names=NULL)-> look

look<- look[, c( 1,4,2,3)]
names( look)<-c("station.id", "elev", "lon", "lat")

USpinfo<- look
  look2<-read.fwf("USmonthly.names.txt",
                     widths=c(2,6,1,12))[,c(2,4)]
  look3<- match(USpinfo$station.id, look2[,1])
USpinfo$station.name<- as.character(look2[,2])[look3]
USpinfo <- as.data.frame( USpinfo)

N<- length( look$station.id)

# file names for ppt

cmd<-paste( "ls ",datadir,  "/NCAR_pinfill/ppt*", sep="")
fnames.ppt<- system( cmd, intern=TRUE)

# now accumulate data
#

# loop over data

USppt<- array(NA, c( 103, 12,N) )
 
for ( k in 1:103){
print( k)
# read in funny formated file for year k
read.fwf( fnames.ppt[k],
       width= c(6,   7,rep(5,11), 3,rep(1,11))) -> dat
 miss<- as.matrix( dat[,(1:12) +13])
 dat<-  as.matrix(  dat[,(1:12) +1 ])
 dat<-  matrix( as.single( dat/10), ncol=12)
 dat[miss==1]<- NA
USppt[k,,] <- dat
  }






save( list=c("UStmax","UStmin","UStinfo", "USppt", "USpinfo"),
               file= "RData.USmonthlyMet.bin")

print( "all done!")

