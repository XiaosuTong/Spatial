######################################
## stl dignostic plotting for fc.E1 
## after spatial loess fitting, check the fitting or remainder
########################################

##############
## vs month ##
##############

library(lattice)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

rst <- rhread("/ln/tongx/Spatial/tmp/tmax/stl/a1950/digno/fc.E1")
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
    local.output, paste(par$dataset, "fitted.vs.month", "ps", sep = ".")
  ), 
  color=TRUE, 
  paper="legal"
)
for(i in 1:128) {
  data <- subset(result, leaf == i)
  a <- xyplot((low+seasonal+middle) ~ time1 | factor
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
          "fitted value"
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
    , main = list(label = "Fitted Value vs. Month")
    , sub = paste(
        "Station from cell",
        i
      )
    , xlim = c(0,119)
    , ylim = c(
        min(min(data$resp, na.rm = TRUE), with(data, min(seasonal+low+middle))) - 1,
        max(max(data$resp, na.rm = TRUE), with(data, max(seasonal+low+middle))) + 1
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

## lowfc component and yearly mean vs.month
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
    local.output, paste(par$dataset, "low.yealymean.vs.month", "ps", sep = ".")
  ),
  color = TRUE, 
  paper = "legal"
)
b <- xyplot( mean ~ time | factor(leaf) 
  , data = result
  , xlab = list(label = "Month", cex = 1.2)
  , ylab = list(label = ylab, cex = 1.2)
  , main = list(label = "Low Frequency Component vs Month")
  , xlim = c(0, 576)
  , pch = 16
  , layout = c(8,1)
  , key=list(
      text = list(label=c("low-frequency component","yearly mean")), 
      lines = list(pch=16, cex=0.7, lwd=1.5, type=c("l","p"), col=col[1:2]),
      columns=2
    )
  , scales = list(
      y = list(relation = 'same'), 
      x=list(at=seq(0, 576, by=120), relation = 'same')
    )
  , prepanel = function(x,y,subscripts,...){
      v <- result[subscripts,]
      ylim <- range(v$mean)
      ans <- prepanel.default.xyplot(v$time, v$low, ...)
      ans$ylim <- range(ans$ylim, ylim)
      ans
    }
  , panel = function(x, y, subscripts, ...) {
      panel.xyplot(
        x = x[seq(1, 576, by=12)],
        y = y[seq(1, 576, by=12)],
        type="p", col=col[2], cex = 0.5, ...
      )
      panel.xyplot(x, result$low[subscripts], type="l", col=col[1], ...)
    }
)
print(b)
dev.off()

## middlefc component and yearly mean-lowfc vs.month
trellis.device(
  device = postscript, 
  file = file.path(
    local.output, paste(par$dataset, "middle.yealymean.vs.month", "ps", sep = ".")
  ),
  color = TRUE, 
  paper = "legal"
)
b <- xyplot( (mean-low) ~ time | factor(leaf)
  , data = result
  , subset = !is.na(mean)
  , xlab = list(label = "Month")
  , ylab = list(label = ylab)
  , main = list(label = "Middle Frequency Component vs Month")
  , xlim = c(0, 576)
  , pch = 16
  , aspect = "xy"
  , layout = c(2,4)
  , key=list(
      text = list(label=c(
        "middle frequency component",
        "yearly mean-low frequency component"
      )),
      lines = list(pch=16, cex=0.7, lwd=1.5, type=c("l","p"), col=col[1:2]), 
      columns = 2
    )
  , scales = list(
      y = list(relation = 'free'), 
      x=list(at=seq(0, 576, by=120), relation = 'same')
    )
  , prepanel = function(x,y,subscripts,...){
      v <- result[subscripts,]
      ans <- prepanel.default.xyplot(v$time, v$middle, ...)
    }
  , panel = function(x, y, subscripts, ...) {
      panel.abline(h=0, col="black", lwd=0.5)
      panel.xyplot(
        x = x[seq(1, 576, by=12)],
        y = y[seq(1, 576, by=12)],
        type="p", col=col[2], cex = 0.5, ...
      )
      panel.xyplot(x, result$middle[subscripts], type="l", col=col[1], ...)
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

###########################################
## Component spatial loess fit dignostic ##
###########################################

library(magrittr)

loess.comp.plot <- function(loess = "fc.E1.loess", comp = "remainder", f = "symmetric", s=0.025) {

  file <- substr(comp, 1, 3) %>% paste("loess.fit", sep=".")

  rst <- rhread(
    file.path(rh.datadir, par$dataset, "stl", "a1950", loess, file, f, paste("sp", s, sep=""))
  )

  tmp <- do.call(rbind, lapply(rst,"[[",1)) %>%
    data.frame(stringsAsFactors=FALSE) 
  tmp$X2 <- tmp$X2 %>% substr(1,3) %>% match(month.abb)
  mod <- tmp %>% with(tmp[order(X1, X2),]) %>% row.names() %>% as.numeric()

  ## a simple function to Capitalize the first character of a string vector
  simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
  }

  if (comp == "remainder") {
    ylab <- simpleCap(comp)
  } else {
    ylab <- paste(simpleCap(comp), "component")
  } 

  mainlab <- paste(simpleCap(comp), "vs. Latitude")

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.lat.lon", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ lat | equal.count(lon, 20, overlap=0) 
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Longitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
           relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Latitude"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.lat.elev", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ lat | equal.count(elev, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Elevation", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Latitude"
      , ylab = ylab 
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()

  mainlab <- paste(simpleCap(comp), "vs. Longitude")

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.lon.lat", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ lon | equal.count(lat, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Latitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Longitude"
      , ylab = ylab    
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.lon.elev", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ lon | equal.count(elev, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Elevation", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Longitude"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()
 
  mainlab <- paste(simpleCap(comp), "vs. Elevation")

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.elev.lat", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ elev2 | equal.count(lat, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Latitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Log (Elevation + 128) (log base 2 meter)"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "vs.elev.lon", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp) ~ elev2 | equal.count(lon, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Longitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(5,4)
      , xlab = "Log (Elevation + 128) (log base 2 meter)"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
      }
    )
  })
  dev.off()

}

loess.resid.plot <- function(loess = "fc.E1.loess", comp = "remainder", f = "symmetric", s=0.025) {

  file <- substr(comp, 1, 3) %>% paste("loess.fit", sep=".")

  rst <- rhread(
    file.path(rh.datadir, par$dataset, "stl", "a1950", loess, file, f, paste("sp", s, sep=""))
  )

  tmp <- do.call(rbind, lapply(rst,"[[",1)) %>%
    data.frame(stringsAsFactors=FALSE) 
  tmp$X2 <- tmp$X2 %>% substr(1,3) %>% match(month.abb)
  mod <- tmp %>% with(tmp[order(X1, X2),]) %>% row.names() %>% as.numeric()

  ## a simple function to Capitalize the first character of a string vector
  simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
  }

  ylab <- "Residuals"

  mainlab <- paste("Residuals vs. Loess Fitted Value of", simpleCap(comp))

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "fit.vs.resid", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp)-fit ~ fit
      , data = data[[r]][[2]]
      , pch = 16
      , cex = 0.3
      , xlab = "Fitted value"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2],evaluation=100,...)
      }
    )
  })
  dev.off()
  
  mainlab <- "Residuals vs. Latitude"

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "resid.vs.lat.lon", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp)-fit ~ lat | equal.count(lon, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Longitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
           relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Latitude"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2],evaluation=100,...)
      }
    )
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "resid.vs.lat.elev", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp)-fit ~ lat | equal.count(elev, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Elevation", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Latitude"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2],evaluation=100,...)
      }
    )
  })
  dev.off()

  mainlab <- "Residuals vs. Longitude"

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "resid.vs.lon.lat", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp)-fit ~ lon | equal.count(lat, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Latitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Longitude"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2], evaluation=100,...)
      }
    )
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "resid.vs.lon.elev", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp)-fit ~ lon | equal.count(elev, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Elevation", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Longitude"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2],evaluation=100,...)
      }
    )
  })
  dev.off()
 
  mainlab <- "Residuals vs. Elevation"

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "resid.vs.elev.lat", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp)-fit ~ elev2 | equal.count(lat, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Latitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Log (Elevation + 128) (log base 2 meter)"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2], evaluation=100,...)
      }
    )
  })
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "resid.vs.elev.lon", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- xyplot( get(comp)-fit ~ elev2 | equal.count(lon, 20, overlap=0)
      , data = data[[r]][[2]]
      , strip=strip.custom(var.name = "Longitude", strip.levels=rep(FALSE, 2))
      , pch = 16
      , cex = 0.3
      , scale = list(
          y = list(
            relation = "free", 
            alternating = TRUE
          ),
          x = list(
            relation = "free",
            tick.number = 3
          )
        )
      , layout = c(10,2)
      , xlab = "Log (Elevation + 128) (log base 2 meter)"
      , ylab = ylab
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x,y,...) {
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
          panel.loess(x,y, span=0.75, degree=1, col=col[2],evaluation=100,...)
      }
    )
  })
  dev.off()

  tmp <- do.call(rbind, lapply(rst,"[[",1)) %>%
    data.frame(stringsAsFactors=FALSE) 
  names(tmp) <- c("year", "month")
  data <- tmp[rep(1:nrow(tmp), each = 4978),] %>% 
    cbind(do.call("rbind", lapply(rst, "[[", 2)))

  Qrst <- ddply(
    .data = data,
    .variable = c("year", "month"),
    .fun = function(r) {
      r <- r[!is.na(r$remainder - r$fit),]
      a <- sort(r$remainder - r$fit)
      idx <- round(seq(1, nrow(r), length.out = 100))
      n.idx <- c(1:idx[2], idx[3:98], idx[99]:length(a))
      f.value <- (n.idx - 0.5) / length(a)
      qnorm <- qnorm(f.value)
      value <- data.frame(
        residual = a[n.idx], 
        qnorm = qnorm, 
        group = as.numeric(n.idx<idx[2] | n.idx>idx[98]),
        fv = f.value
      )
    }
  )
  Qrst$month <- factor(
  Qrst$month, 
    levels = c(
      "Jan","Feb","Mar","Apr","May","June",
      "July","Aug", "Sep", "Oct", "Nov", "Dec"
    )
  )
  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "QQ.resid.month", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  for(i in 1950:1997){
    a <- xyplot(residual ~ qnorm | month*as.factor(year)
      , data = Qrst
      , subset = year == paste(i)
      , layout = c(6,2)
      , pch  = 16
      , cex  = 0.3
      , key=list(
          text = list(label=c("quantiles in [1.5%, 98.5%]","quantiles outside [1.5%, 98.5%]")), 
          lines = list(pch=16, cex=0.7, lwd=1.5, type=c("p","p"), col=col[1:2]),
          columns=2
        )
      , xlab = list(label = "Unit normal quantile")
      , ylab = list(label = "Residuals")
      , main = "Normal Quantiles of Loess Residuals of Remainder"
      , panel = function(x, y, subscripts, ...) {
          tmp <- subset(Qrst[subscripts, ], year == paste(i))
          panel.qqmathline(tmp[tmp$group == 0, "residual"], ...)
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(tmp[tmp$group == 0, "qnorm"],tmp[tmp$group == 0, "residual"], col=col[1],...)
          panel.xyplot(tmp[tmp$group == 1, "qnorm"],tmp[tmp$group == 1, "residual"], col=col[2],...)
        }
    )
    print(a)
  }
  dev.off()

  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "QQ.resid.month2", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
    a <- xyplot(residual ~ qnorm | month*as.factor(year),
      , data = Qrst, 
      , subset = group == 0,
      , layout = c(6,2),
      , pch  = 16,
      , cex  = 0.3,
      , xlab = list(label = "Unit normal quantile"),
      , ylab = list(label = "Residuals"),
      , main = "Normal Quantiles of Loess Residuals of Remainder"
      , panel = function(x,y, ...) {
          panel.qqmathline(y,y=y,...)
          panel.abline(h=0, lwd=0.5, col="black")
          panel.xyplot(x,y,...)
        }
    )
    print(a)
  dev.off()

}