####################################
## stl dignostic plotting for E28 ##
####################################

#############
## vs month #
#############
library(lattice)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

rst <- rhread("/ln/tongx/Spatial/tmp/tmax/stl/a1950/digno/E28")
result <- rst[[1]][[2]]

if(par$dataset == "tmax"){
    ylab <- "Maximum Temperature (degrees centigrade)"
}else if(par$dataset == "tmin"){
    ylab <- "Minimum Temperature (degrees centigrade)"
}else {
    ylab <- "Precipitation (millimeters)"
}
##  trend + seasonal component with original obs vs. month
trellis.device(
  device = postscript, 
  file = file.path(
    local.output, paste(par$dataset, "trend+seasonal.vs.month", "ps", sep = ".")
  ), 
  color=TRUE, 
  paper="legal"
)
for(i in 1:128) {
  data <- subset(result, leaf == i)
  a <- xyplot((trend+seasonal) ~ time1 | factor
    , data = data
    , layout = c(1,5)
    , strip = FALSE
    , aspect = "xy"
    , scale = list(
        y = list(
          relation = 'same', 
          tick.number=4, 
          alternating=TRUE
        ), 
        x = list(
          at = seq(0, 143, by=12), 
          relation = 'same'
        )
      )
    , key = list(
        text = list(label=c(
          "raw obs",
          "trend+seasonal component"
        )),
        lines = list(
          pch = 1, 
          cex = 0.5, 
          lwd = 1.5, 
          type = c("p","l"), 
          col = col[1:2]
        ), 
        columns = 2
      )
    , main = list(label = "Trend and Seasonal component vs. Month")
    , sub = paste(
        "Station from cell",
        i
      )
    , xlim = c(0,119)
    , ylim = c(
        min(min(data$resp, na.rm = TRUE), min(data$seasonal+data$trend)) - 1,
        max(max(data$resp, na.rm = TRUE), max(data$seasonal+data$trend)) + 1
      )
    , ylab = list(label = ylab)
    , xlab = list(label = "Month")
    , panel = function(x, y, subscripts ,...) {
        panel.xyplot(
          x = sort(x), 
          y = y[order(x)], 
          type="b", col=col[2], pch=16,cex =0.5,lwd=1, ...
        )
        panel.xyplot(
          x = data[subscripts,]$time1, 
          y = data[subscripts,]$resp,
          type="p", col=col[1], cex=0.5, ...
        )
        panel.abline(
          v=seq(0,119, by=12), color="lightgrey", lty=3, lwd=0.5
        )
      }
  )
  print(a)
}
dev.off()

##  remainder vs. month 
trellis.device(
  device = postscript, 
  file = file.path(
    local.output, paste(par$dataset, "remainder.vs.month", "ps", sep = ".")
  ), 
  color = TRUE,
  paper = "legal"
)
for(i in 1:128) {
  a <- xyplot(remainder ~ time1 | factor,
    , data = result
    , subset = leaf == i
    , layout = c(1,5)
    , strip = FALSE
    , pch  = 16
#   , aspect = "xy"
    , cex  = 0.4
    , xlim = c(0, 120)
    , main = list(label = "Remainder vs. Month")
    , ylab = list(label = ylab)
    , sub = paste("Station from cell", i)
    , xlab = list(label = "Month")
    , panel = function(x, y,...) {
        panel.abline(h=0, col="red", lty=1, lwd=0.5)
        panel.xyplot(x, y,...)
      }
  )
  print(a)
}
dev.off()
trellis.device(
  device = postscript, 
  file = file.path(
    local.output, paste(par$dataset, "remainder2.vs.month", "ps", sep = ".")
  ), 
  color = TRUE,
  paper = "legal"
)
for(i in 1:128) {
  a <- xyplot(remainder ~ time,
    , data = result
    , subset = leaf == i & !is.na(remainder)
    , layout = c(1,1)
    , strip = FALSE
    , pch  = 16
#   , aspect = "xy"
    , cex  = 0.4
    , key = list(
        text = list(label=c(
          "remainder",
          "loess soomthing: span=0.15, degree=2",
          "loess soomthing: span=0.35, degree=1"
        )),
        lines = list(
          pch = 16, 
          cex = 0.5, 
          lwd = 1.5, 
          type = c("p","l","l"), 
          col = col[1:3]
        ), 
        columns = 3
      )
    , scale = list(x = list(relation="free"))
    , main = list(label = "Remainder vs. Month")
    , ylab = list(label = ylab)
    , sub = paste("Station from cell", i)
    , xlab = list(label = "Month")
    , panel = function(x, y,...) {
        panel.abline(h=0, col="black", lty=1, lwd=0.5)
        panel.xyplot(x, y,...)
        panel.loess(x,y,degree=2,span=0.15, col=col[2], evaluation=200,...)
        panel.loess(x,y,degree=1,span=0.35, col=col[3], evaluation=200,...)
      }
  )
  print(a)
}
dev.off()

## QQ plot of remainder for each station
trellis.device(
  device = postscript, 
  file = file.path(
    local.output, paste(par$dataset, "QQ.remainder", "ps", sep = ".")
  ), 
  color = TRUE,
  paper = "legal"
)
a <- qqmath(~ remainder | factor(leaf)
  , data = result
  , distribution = qnorm
  , layout = c(5,3)
  , pch  = 16
  , aspect = 1
  , cex  = 0.3
  , main = list(label = "Normal Quantiles of Remainder")
  , scale = list(y=list(relation="free"))
  , xlab = list(label = "Unit normal quantile")
  , ylab = list(label = "Remainder")
# , prepanel = prepanel.qqmathline
  , panel = function(x, ...) {
      panel.abline(v=seq(-4,4,2), h=seq(-9,9,3), col="lightgrey", lty=1, lwd=0.5)
      panel.qqmathline(x, y=x)
      panel.qqmath(x, ...)
    }
)
print(a)
dev.off()

## trend component and yearly mean vs.month
dr <- ddply(
  .data = result,
  .variable = c("leaf", "year"),
  .fun = summarise,
  mean = mean(resp, na.rm=TRUE)
)
result<- result[order(result$leaf),]
result <- dr[rep(row.names(dr), each=12), "mean", drop=FALSE] %>% cbind(result)
trellis.device(
  device = postscript, 
  file = file.path(
    local.output, paste(par$dataset, "trend.yealymean.vs.month", "ps", sep = ".")
  ),
  color = TRUE, 
  paper = "legal"
)
b <- xyplot( mean ~ time | factor(leaf) 
  , data = result
  , xlab = list(label = "Month", cex = 1.2)
  , ylab = list(label = ylab, cex = 1.2)
  , xlim = c(0, 576)
  , pch = 16
  , layout = c(5,2)
#  , aspect="xy"
  , key=list(
      text = list(label=c("trend component","yearly mean")), 
      lines = list(pch=16, cex=0.7, lwd=1.5, type=c("l","p"), col=col[1:2]),
      columns=2
    )
  , scales = list(
      y = list(relation = 'free'), 
      x=list(at=seq(0, 576, by=120), relation = 'same')
    )
  , prepanel = function(x,y,subscripts,...){
      v <- result[subscripts,]
      ylim <- range(v$mean)
      ans <- prepanel.default.xyplot(v$time, v$trend, ...)
      ans$ylim <- range(ans$ylim, ylim)
      ans
    }
  , panel = function(x, y, subscripts, ...) {
      panel.xyplot(
        x = x[seq(1, 576, by=12)],
        y = y[seq(1, 576, by=12)],
        type="p", col=col[2], cex = 0.5, ...
      )
      panel.xyplot(x, result$trend[subscripts], type="l", col=col[1], ...)
    }
)
print(b)
dev.off()

##########################
## conditional on month ##
##########################

## seasonal component vs month for each station
## seasonal is a constant for each month because of periodic
trellis.device(
  device = postscript, 
  file = file.path(
    local.output, paste(par$dataset, "seasonal.vs.month", "ps", sep = ".")
  ),  
  color = TRUE, 
  paper = "legal"
)
  b <- xyplot( seasonal ~ month | factor(leaf)
    , data = result
    , subset = year == 1950
    , xlab = list(label = "Month in year", cex = 1.2)
    , ylab = list(label = ylab, cex = 1.2)
    , main = "Seasonal Component vs Month in Year"
    , type = "b"
    , pch = 16
    , aspect = "xy"
    , cex = 0.5
    , layout = c(5, 4),
    , scales = list(
        y = list(relation = 'same', alternating=TRUE), 
        x = list(at=c(1, 3, 5, 7, 9, 11), relation='same')
      )
    , panel = function(x,y,...) {
        panel.xyplot(x,y,...)
      }
  )
  print(b)
dev.off()

## remainder vs year condtional on month in year
trellis.device(
  device = postscript, 
  file = file.path(
    local.output, paste(par$dataset, "remainder.vs.year", "ps", sep = ".")
  ),    
  color = TRUE, 
  paper = "legal"
)
for(i in 1:128){
  b <- xyplot( remainder ~ (year-1949) | month
    , data = result
    , subset = leaf == i 
    , xlab = list(label = "Year", cex = 1.2)
    , ylab = list(label = ylab, cex = 1.2)
    , main = "Remainder vs Year"
    , sub = paste("Station from cell", i)
    , pch = 16
    , cex = 0.5
    , key=list(
        text = list(label=c("remainder","loess smoothing: span=0.85, degree=1")), 
        lines = list(pch=16, cex=0.7, lwd=1.5, type=c("l","p"), col=col[1:2]),
        columns=2
      )
    , layout = c(12,1)
    , strip = TRUE
    , scales = list(
        y = list(relation = 'same', alternating=TRUE), 
        x = list(tick.number=10, relation='same')
      )
    , panel = function(x,y,...){
        panel.abline(h=0, color="black", lty=1, lwd=0.5)
        panel.xyplot(x,y,...)
        panel.loess(x,y,span=3/4, degree=1, col=col[2],...)
      }
  )
  print(b)
}
dev.off()

## Normal quantiles of remainder conditional on month in year
trellis.device(
  device = postscript, 
  file = file.path(
    local.output, paste(par$dataset, "QQ.remainder.month", "ps", sep = ".")
  ),     
  color = TRUE, 
  paper = "legal"
)
for(i in 1:128){
  a <- qqmath(~ remainder | month,
    , data = result
    , subset = leaf == i
    , distribution = qnorm
    , aspect = 1
    , pch = 16
    , cex = 0.4
    , layout = c(4,3)
    , main = paste("Normal Quantiles of Remainder for Station from Cell", i)
    , xlab = list(label="Unit normal quantile")
    , ylab = list(label = "Remainder")
    , prepanel = prepanel.qqmathline,
    , panel = function(x, y,...) {
        panel.grid(lty=3, lwd=0.5, col="black",...)
        panel.qqmathline(x, y=x)
        panel.qqmath(x, y,...)
      }
  )
  print(a)
}
dev.off()
