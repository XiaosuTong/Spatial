###################################################
##QQ plot of year for each station
###################################################

#load the data
library(lattice)
library(maps)
datadir <- "~/Projects/Spatial/NCAR/RData/"
load(paste(datadir,"USmonthlyMet.RData", sep=""))
#precip <- USppt[!(is.na(USppt$precip)),1:7]

tmax <- UStemp[!(is.na(UStemp$tmax)),1:7]
tmin <- UStemp[!(is.na(UStemp$tmin)),1:7]


#month <- precip$month
#levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
#month <- as.numeric(factor(month, levels=c(1:12)))
#precip$time <- precip$year-1895 + month/12

month <- tmax$month
levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
month <- as.numeric(factor(month, levels=c(1:12)))
tmax$time <- tmax$year-1895 + month/12

month <- tmin$month
levels(month) <- c(4,8,12,2,1,7,6,3,5,11,10,9)
month <- as.numeric(factor(month, levels=c(1:12)))
tmin$time <- tmin$year-1895 + month/12

#Pcount <- sort(tapply(!is.na(USppt$precip), USppt$station.id, sum), decreasing=TRUE)
TMaxcount <- sort(tapply(!(is.na(UStemp$tmax)), UStemp$station.id, sum), decreasing=TRUE)
TMincount <- sort(tapply(!(is.na(UStemp$tmin)), UStemp$station.id, sum), decreasing=TRUE)


#precip$station.id <- factor(precip$station.id, levels=names(Pcount))
tmax$station.id <- factor(tmax$station.id, levels=names(TMaxcount))
tmin$station.id <- factor(tmin$station.id, levels=names(TMincount))


trellis.device(postscript, file = paste("QQ_plot_of_year_for_precipitation", ".ps", sep = ""), color=TRUE, paper="legal")
        a <- qqmath(~ time | station.id,
                data = precip,
                distribution = qunif,
                aspect = 1,
		pch = 16,
		cex = 0.5,
             	strip = FALSE,		
		layout = c(10,6),
                xlab = list(label="f-value", cex=1.5),
                ylab = list(label="Year", cex=1.5),
                scales = list(x = list(cex=1.5), y = list(cex=1.5)),
                prepanel = prepanel.qqmathline,
                panel = function(x, y,...) {
                        panel.grid()
                        panel.qqmath(x, y,...)
                }

        )
        print(a)
dev.off()


trellis.device(postscript, file = paste("QQ_plot_of_year_for_tmax", ".ps", sep = ""), color=TRUE, paper="legal")
        a <- qqmath(~ time | station.id,
                data = tmax,
                distribution = qunif,
                aspect = 1,
		pch = 16,
		cex = 0.5,
                strip = FALSE,
                layout = c(10,6),
                xlab = list(label="f-value", cex=1.5),
                ylab = list(label="Year", cex=1.5),
                scales = list(x = list(cex=1.5), y = list(cex=1.5)),
                prepanel = prepanel.qqmathline,
                panel = function(x, y,...) {
                        panel.grid()
                        panel.qqmath(x, y,...)
                }

        )
        print(a)
dev.off()

trellis.device(postscript, file = paste("QQ_plot_of_year_for_tmin", ".ps", sep = ""), color=TRUE, paper="legal")
        a <- qqmath(~ time | station.id,
                data = tmin,
                distribution = qunif,
                aspect = 1,
                pch = 16,
                cex = 0.5,
                strip = FALSE,
                layout = c(10,6),
                xlab = list(label="f-value", cex=1.5),
                ylab = list(label="Year", cex=1.5),
                scales = list(x = list(cex=1.5), y = list(cex=1.5)),
                prepanel = prepanel.qqmathline,
                panel = function(x, y,...) {
                        panel.grid()
                        panel.qqmath(x, y,...)
                }

        )
        print(a)
dev.off()


 
