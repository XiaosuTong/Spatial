q########################################################
## spatial loess fit of stl remainder at stations after 1950 
##
## - a1950/Exp.fc/E1-E3
##     different setting for stl with second fc, read in bystation subsets
##     only includes 4,978 stations that have 300+ obs after 1950
##     key is the station.id, value is the data.frame includes stl fit
##     lon, lat, and elev have been saved as attributes in data.frame
##
## - a1950/fc.quantiles
##     read in all E1-E30 files, for each station calculate the quantiles
##     key is the quantiles, value is a data.frame has station.id, idx which
##     identifies which E it is, and resp
##
## - plot remainders for each quantile
##     best setting: sw="periodic", tw=865, td=1, fcw=c(865,109), fcd=c(1,2)
##
## - a1950/digno/fc.E1
##     read in a1950/Exp.fc/E1, kept only 128 stations
##
## - stl+ dignostic plot
##     source rh.a1950.stlplot02 file
##  
## - comp.smooth(), same as the one in stlloess01.R
##     input file is a1950/Exp.fc/E1 which is the best fitting from stl w/ fc component.
##     output file are a1950/fc.E1.loess/comp.E1-E2, which is the evalutation of component
##     at grid points, then will be used for levelplot only.
##     
## - check.smooth(), same as the one in stlloess01.R
##     input is a1950/E28.loess/comp.E1-E2, plot the levelplot for given component.
##
## - loess.fit(), same as the one in stlloess01.R
##     The loess fitting for Remainder. input file is a1950/E28 and key is station.id
##     output is a1950/fc.E1.loess/comp.loess.fit/family/span and key is month
## 
########################################################

source("~/Rhipe/rhinitial.R")
par <- list()
par$machine <- "gacrux"
par$dataset <- "tmax"
par$N <- 576
par$family <- "symmetric"
par$arg <- list( 
  list(
    sw = "periodic", tw = 865, td = 1, sd = 1, 
    ffcw = 865, ffcd = 1, sfcw = 109, sfcd = 2,
    multiple = TRUE
  ),
  list(
    sw = "periodic", tw = 865, td = 1, sd = 1,
    ffcw = 865, ffcd = 1, sfcw = 241, sfcd = 2,
    multiple = TRUE
  ),
  list(
    sw = "periodic",
    tw = 109, 
    td = 2, 
    sd = 1,
    multiple = FALSE
  )
)
source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")

# load the stations vector after 1950, there are 7,738 stations
load(file.path(local.datadir, "stations.a1950.RData"))
month.or <- c(
    "Jan","Feb","Mar","Apr","May","June",
    "July","Aug", "Sep", "Oct", "Nov", "Dec"
)
source("~/Projects/Spatial/NCAR/myloess/my.loess02.R")
source("~/Projects/Spatial/NCAR/myloess/my.predloess.R")


##########################################################
## Read in bystation subsets, Experiment w/ fc, E1-E3   ##
##########################################################

for (i in 1:length(par$arg)) {
job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
    if (map.keys[[r]] %in% stations) {         # only includes 7,738 stations
      rhcounter("stations", "_a1950_", 1)
      value <- subset(map.values[[r]], year >= 1950)
      if (sum(!is.na(value$resp)) >= 300) {   # only includes 4,978 stations each of which has 300+ obs
        rhcounter("stations", "_300_", 1)
        value$month <- factor(value$month, levels = month.or)
        value <- value[with(value, order(year, month)),]
        if (nrow(value) == 576) rhcounter("stations", "_576_", 1)
        value$time <- 1:nrow(value)
        if(flag) {
          fit <- do.call("cbind",
            stl2(
              x = value$resp, 
              t = value$time, 
              n.p = 12, 
              s.window = sw, 
              s.degree = sd, 
              t.window = tw, 
              t.degree = td,
              fc.window = c(ffcw, sfcw), 
              fc.degree = c(ffcd, sfcd), 
              inner = inner, 
              outer = outer
            )[c("data","fc")]
          )
          names(fit)[grep("fc.fc", names(fit))] <- c("low", "middle")
          names(fit)[grep("data.seasonal", names(fit))] <- "seasonal"
          names(fit)[grep("fc.remainder", names(fit))] <- "remainder"
          outvalue <- cbind(
            subset(value,  select = c(month, year, resp, time)),
            subset(fit, select = -c(data.trend, data.remainder, data.weights, data.sub.labels, data.raw))
          )
        } else {
          fit <- stl2(
            x = value$resp, 
            t = value$time, 
            n.p = 12, 
            s.window = sw, 
            s.degree = sd, 
            t.window = tw, 
            t.degree = td, 
            inner = inner, 
            outer = outer
          )$data
          outvalue <- cbind(
            subset(value,  select = c(month, year, resp, time)),
            subset(fit, select = -c(weights, sub.labels, raw))
          )
        }
        attributes(outvalue)$elev  <- unique(as.character(value$elev))
        attributes(outvalue)$lon <- unique(as.character(value$lon))
        attributes(outvalue)$lat <- unique(as.character(value$lat))
        rhcollect(map.keys[[r]], outvalue)
      }
    }
  })
})
job$setup <- expression(
  map = {
    library(lattice)
    library(yaImpute, lib.loc = lib.loc)
    library(stl2, lib.loc = lib.loc)
  }
)
job$cleanup <- expression(
  map = {
    rm(list=ls())
  }
)
if (par$arg[[i]]$multiple) { 
  job$parameters <- list(
    stations = stations.a1950.tmax,
    sw = par$arg[[i]]$sw,
    sd = par$arg[[i]]$sd,
    tw = par$arg[[i]]$tw,
    td = par$arg[[i]]$td,
    ffcw = par$arg[[i]]$ffcw,
    sfcw = par$arg[[i]]$sfcw,
    ffcd = par$arg[[i]]$ffcd,
    sfcd = par$arg[[i]]$sfcd,
    inner = 10,
    outer = 0,
    flag = par$arg[[i]]$multiple
  )
} else {
  job$parameters <- list(
    stations = stations.a1950.tmax,
    sw = par$arg[[i]]$sw,
    sd = par$arg[[i]]$sd,
    tw = par$arg[[i]]$tw,
    td = par$arg[[i]]$td,
    inner = 10,
    outer = 0,
    flag = par$arg[[i]]$multiple
  )
}
job$input <- rhfmt(
  file.path(rh.datadir, par$dataset, "data", "bystation"),
  type = "sequence"
) 
job$output <- rhfmt(
  file.path(rh.datadir, par$dataset, "stl", "a1950", "Exp.fc", paste("E", i, sep="")), 
  type = "sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 72,
  mapred.tasktimeout = 0
)
job$mon.sec <- 10
job$jobname <- file.path(
  rh.datadir, par$dataset, "stl", "a1950", "Exp.fc", paste("E", i, sep="")
)
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)
}

##################################################
## get the quantiles for each setting in Exp.fc ##
##################################################

job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
    file <- Sys.getenv("mapred.input.file")
    idx <- substr(strsplit(tail(strsplit(file, "/")[[1]],2), "[.]")[[1]], 2, 3)
    sumy <- quantile(
      x = map.values[[r]]$remainder, 
      probs = c(0.05, 0.25, 0.5, 0.75, 0.95), 
      na.rm = TRUE
    )
    mean <- mean(map.values[[r]]$remainder, na.rm = TRUE)
    sumy <- c(sumy, mean)
    lapply(seq_along(sumy), function(k, ind = idx) {
      value <- data.frame(
        station = map.keys[[r]],  
        resp = sumy[k]
      )
      value <- cbind(value, id = ind)
      rhcollect(k, value)
    })
  })
})
job$reduce <- expression(
  pre = {
    combine <- data.frame()
  },
  reduce = {
    combine <- rbind(combine, do.call(rbind, reduce.values))
  },
  post = {
    rhcollect(reduce.key, combine)
  }
)
job$input <- rhfmt(
  file.path(rh.datadir, par$dataset, "stl", "a1950", "Exp.fc"), type="sequence"
)
job$output <- rhfmt(
  file.path(rh.datadir, par$dataset, "stl", "a1950", "fc.quantiles"), type="sequence"
)
job$mapred <- list(
  mapred.reduce.tasks = 6,
  mapred.tasktimeout = 0,
  rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 10
job$jobname <- file.path(rh.datadir, par$dataset, "stl", "a1950", "fc.quantiles")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)

##############################################
## tunning parameters for different setting ##
##############################################

## have to vary the group variable in the plotting function
## then condition on the rest of variables

library(lattice)
lattice.theme <- trellis.par.get()
col <- lattice.theme$superpose.symbol$col

rst <- rhread("/ln/tongx/Spatial/tmp/tmax/stl/a1950/fc.quantiles")
myfun <- function(sub) {
trellis.device(
  postscript, 
  file = paste(
    local.output, "/a1950.quant.", sub[[1]], ".", #with elev at the end of the plot name 
    par$dataset, ".ps", sep = ""
  ), 
  color = TRUE,
  paper = "legal"
)
  b <- qqmath(~ resp
    , data = sub[[2]]
    , pch = 16
    , distribution = qunif
    , cex = 0.2
    , aspect = 1
    , group = factor(id, levels=c(1,2,3))
    , xlab = list(label = "f-value")
    , ylab = list(label = "Remainder")
    , key=list(
        text = list(label=c("middle-span=109", "middle-span=241", "no middel")),
        lines = list(
          pch=16, 
          cex=0.7, 
          lwd=1.5, 
          type=c("p","p","p"), 
          col=col[1:3]
        ), 
        columns = 3
      )
    , panel = function(x,...) {
      panel.qqmath(x,...)
      panel.abline(h=0, col="black", lwd=0.5)
    }
  )
  print(b)
dev.off()
}
lapply(rst, myfun)

#########################################
## only keep 128 stations in digno/E28 ##
#########################################

##  load which 128 stations 
load(file.path(local.datadir, "samplestation.a1950.RData"))
##  read in the best parameter setting only for 128 stations
job <- list()
job$map <- expression({
  lapply(seq_along(map.values), function(r) {
    if(map.keys[[r]] %in% stations$station.id) {
      rhcounter("stations", "_300_", 1)
      value <- map.values[[r]]
      #value$station.id <- map.keys[[r]]
      value$leaf <- stations$leaf[which(stations$station.id == map.keys[[r]])]
      value$time1 <- c(rep(0:119, 4), 0:95)
      value$factor <- factor(
        rep(1:5, c(rep(120, 4), 96)),
        levels = c(5:1)
      )
      rhcollect(1, value)
    }
  })
})
job$reduce <- expression(
  pre = {
    combine <- data.frame()
  },
  reduce = {
    combine <- rbind(combine, do.call(rbind, reduce.values))
  },
  post = {
    rhcollect(reduce.key, combine)
  }
)
job$input <- rhfmt(
  file.path(rh.datadir, par$dataset, "stl", "a1950", "Exp.fc", "E1"), type="sequence"
)
job$output <- rhfmt(
  file.path(rh.datadir, par$dataset, "stl", "a1950", "digno", "fc.E1"), type="sequence"
)
job$parameters <- list(
  stations = tmax.sample.a1950
)
job$mapred <- list(
  mapred.reduce.tasks = 1,
  mapred.tasktimeout = 0,
  rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 10
job$jobname <- file.path(rh.datadir, par$dataset, "stl", "a1950", "digno", "fc.E1")
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)


###########################################################
## read in Exp.fc/E1 and change to be key=c(month, year) ##
###########################################################

library(maps)
library(magrittr)
def.grid <- expand.grid(
  lon = seq(-124, -67, by = 0.25),
  lat = seq(25, 49, by = 0.25)
)
instate <- !is.na(map.where("state", def.grid$lon, def.grid$lat))
def.grid <- def.grid[instate, ]

## E1 is span=0.025, E2 is span=0.01
par$loess <- list(
  family = "symmetric",
  span = 0.01,
  degree = 2
)
comp.smooth <- function(component = "middle", input.stl = "Exp.fc", input.Exp = "E1", out.loess = "fc.E1.loess", out.Exp = "E2", arg = par$loess, new.grid = def.grid) { 

  file <- substr(component, 1, 3) %>% paste(out.Exp, sep=".")

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      if (comp == "seasonal") { ## for seasonal we only need one year info
        d_ply(
        .data = map.values[[r]],
        .variable = c("year","month"),
        .fun = function(k, station=map.keys[[r]]) {
          key <- c(unique(k$year), unique(as.character(k$month)))
          value <- k[, c(comp, "remainder"), drop=FALSE]
          value$station.id <- station
          value$lon <- as.numeric(attributes(map.values[[r]])$lon)
          value$lat <- as.numeric(attributes(map.values[[r]])$lat)
          value$elev <- as.numeric(attributes(map.values[[r]])$elev)
          if (key[1] == "1950") {
            rhcounter("year", "_1950_", 1)    
            rhcollect(key, value)
          }
        })
      } else {  ## for the rest of component we neen 576 month all
        d_ply(
        .data = map.values[[r]],
        .variable = c("year","month"),
        .fun = function(k, station=map.keys[[r]]) {
          key <- c(unique(k$year), unique(as.character(k$month)))
          value <- k[, c(comp, "remainder"), drop=FALSE]
          value$station.id <- station
          value$lon <- as.numeric(attributes(map.values[[r]])$lon)
          value$lat <- as.numeric(attributes(map.values[[r]])$lat)
          value$elev <- as.numeric(attributes(map.values[[r]])$elev)
          rhcollect(key, value)
        })
      }
    })
  })

  job$reduce <- expression(
    pre = {
      combine <- data.frame()
    },
    reduce = {
      combine <- rbind(combine, do.call(rbind, reduce.values))
    },
    post = {
        lo.fit <- my.loess2( get(comp) ~ lon + lat, 
          data    = subset(combine, !is.na(combine$remainder)), 
          degree  = argumt$degree, 
          span    = argumt$span,
          family  = argumt$family,
          control = loess.control(surface = "direct")
        )
        fit <- my.predict.loess(
          object = lo.fit, 
          newdata = data.frame(
            lon = grid.in$lon,
            lat = grid.in$lat
          )
        )
        grid.in$fit <- fit
        rhcollect(reduce.key, grid.in)
    }
  )

  job$parameters <- list(
    argumt = arg,
    my.loess2 = my.loess2,
    my.simple2 = my.simple2,
    my.predict.loess = my.predict.loess,
    my.predLoess = my.predLoess,
    grid.in = new.grid,
    comp = component
  )

  job$setup <- expression(
    map = {
      library(plyr)
    },
    reduce = {
      library(maps, lib.loc = lib.loc)
      system("chmod 777 myloess2.so")
      dyn.load("myloess2.so")
    }
  )

  job$shared <- c(
    file.path(
      rh.datadir, par$dataset, "shareRLib", "myloess2.so"
    )
  )

  job$mapred <- list(
    mapred.reduce.tasks = 144,
    mapred.tasktimeout = 0,
    rhipe_reduce_buff_size = 10000
  )

  job$input <- rhfmt(
    file.path(rh.datadir, par$dataset, "stl", "a1950", input.stl, input.Exp), type="sequence"
  )

  job$output <- rhfmt(
    file.path(rh.datadir, par$dataset, "stl", "a1950", out.loess, file), type="sequence"
  )

  job$mon.sec <- 10

  job$jobname <- file.path(
    rh.datadir, par$dataset, "stl", "a1950", out.loess, file
  )

  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)
}

########################################
## check the smoothness of components ##
########################################

library(magrittr)
library(maps)

check.smooth <- function(loess = "fc.E1.loess", comp = "remainder", Exp = "E2", n.cut = 100) {

  us.map <- map('state', plot = FALSE, fill = TRUE)

  file <- substr(comp, 1, 3) %>% paste(Exp, sep=".")
  rst <- rhread(
    file.path(rh.datadir, par$dataset, "stl", "a1950", loess, file)
  )

  if (comp == "seasonal") {
    mod <- unlist(lapply(rst, function(r) r[[1]][2])) %>%
      substr(1,3) %>%
      match(month.abb) %>%
      order()
  } else {
    tmp <- do.call(rbind, lapply(rst,"[[",1)) %>%
      data.frame(stringsAsFactors=FALSE) 
    tmp$X2 <- tmp$X2 %>% substr(1,3) %>% match(month.abb)
    mod <- tmp %>% with(tmp[order(X1, X2),]) %>% row.names() %>% as.numeric()
  }

  ## a simple function to Capitalize the first character of a string vector
  simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
      sep="", collapse=" ")
  }
  
  if (comp == "remainder") {
    mainlab <- paste("Smoothing", simpleCap(comp), "for")
  } else if (comp == "seasonal") {
    mainlab <- paste("Smoothing", simpleCap(comp), "Component for")
  } else {
    mainlab <- paste("Smoothing", simpleCap(comp), "Frequency Component for")
  }
  trellis.device(
    device = postscript, 
    file = file.path(
      local.output, paste(par$dataset, comp, "smooth", "ps", sep = ".")
    ),
    color = TRUE, 
    paper = "legal"
  )
  lapply(mod, function(r, data = rst) {
    b <- levelplot(fit ~ lon * lat
      , data = data[[r]][[2]]
      , region = TRUE
 #    , at = c(seq(-4,-1,0.2),seq(-0.8,0.8,0.1), seq(1,4,0.2))
      , cuts = n.cut
      , col.regions = colorRampPalette(c("blue", "yellow","red"))
 #    , colorkey = list(at = seq(-1, 1, 0.1))
      , xlab = "Longitude"
      , ylab = "Latitude"
      , main = paste(mainlab, data[[r]][[1]][1], data[[r]][[1]][2])
      , panel = function(x, y, z, ...) {
          panel.levelplot(x,y,z,...)
          panel.polygon(us.map$x,us.map$y, border = "black") 
        }
    )
    print(b)
  })
  dev.off()

}


#####################################
## Spatial loess fit for component ##
#####################################

loess.fit <- function(component = "seasonal", input.stl = "Exp.fc", input.Exp = "E1", out.loess = "fc.E1.loess", f = "symmetric", s = 0.025, d = 2) {

  us.map <- map('state', plot = FALSE, fill = TRUE)
  file <- substr(component, 1, 3) %>% paste("loess.fit", sep=".")
  arg <- list(span=s, degree=d, family = f)

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      if (comp == "seasonal") { ## for seasonal we only need one year info
        d_ply(
        .data = map.values[[r]],
        .variable = c("year","month"),
        .fun = function(k, station=map.keys[[r]]) {
          key <- c(unique(k$year), unique(as.character(k$month)))
          value <- k[, c(comp, "remainder"), drop=FALSE]
          value$station.id <- station
          value$lon <- as.numeric(attributes(map.values[[r]])$lon)
          value$lat <- as.numeric(attributes(map.values[[r]])$lat)
          value$elev <- as.numeric(attributes(map.values[[r]])$elev)
          if (key[1] == "1950") {
            rhcounter("year", "_1950_", 1)    
            rhcollect(key, value)
          }
        })
      } else {  ## for the rest of component we need 576 month all
        d_ply(
        .data = map.values[[r]],
        .variable = c("year","month"),
        .fun = function(k, station=map.keys[[r]]) {
          key <- c(unique(k$year), unique(as.character(k$month)))
          value <- k[, c(comp, "remainder"), drop=FALSE]
          value$station.id <- station
          value$lon <- as.numeric(attributes(map.values[[r]])$lon)
          value$lat <- as.numeric(attributes(map.values[[r]])$lat)
          value$elev <- as.numeric(attributes(map.values[[r]])$elev)
          rhcollect(key, value)
        })
      }
    })
  })

  job$reduce <- expression(
    pre = {
      combine <- data.frame()
    },
    reduce = {
      combine <- rbind(combine, do.call(rbind, reduce.values))
    },
    post = {
      combine$elev2 <- log2(combine$elev+128)
      lo.fit <- my.loess2( get(comp) ~ lon + lat + elev2, 
        data    = subset(combine, !is.na(remainder)), 
        degree  = argumt$degree, 
        span    = argumt$span,
        family  = argumt$family,
        parametric = "elev2",
        control = loess.control(surface = "direct")
      )
      fit <- my.predict.loess(
        object = lo.fit, 
        newdata = data.frame(
          lon = combine$lon, 
          lat = combine$lat,
          elev2 = combine$elev2
        )
      )
      combine$fit <- fit
      rhcollect(reduce.key, combine)
    }
  )

  job$parameters <- list(
    argumt = arg,
    my.loess2 = my.loess2,
    my.simple2 = my.simple2,
    my.predict.loess = my.predict.loess,
    my.predLoess = my.predLoess,
    comp = component
  )

  job$setup <- expression(
    map = {
      library(plyr)
    },
    reduce = {
      library(maps, lib.loc = lib.loc)
      system("chmod 777 myloess2.so")
      dyn.load("myloess2.so")
    }
  )

  job$shared <- c(
    file.path(
      rh.datadir, par$dataset, "shareRLib", "myloess2.so"
    )
  )

  job$mapred <- list(
    mapred.reduce.tasks = 144,
    mapred.tasktimeout = 0,
    rhipe_reduce_buff_size = 2500
  )

  job$input <- rhfmt(
    file.path(rh.datadir, par$dataset, "stl", "a1950", input.stl, input.Exp), type="sequence"
  )

  job$output <- rhfmt(
    file.path(rh.datadir, par$dataset, "stl", "a1950", out.loess, file, f, paste("sp", s, sep="")), type="sequence"
  )

  job$mon.sec <- 10

  job$jobname <- file.path(
    rh.datadir, par$dataset, "stl", "a1950", out.loess, file, f, paste("sp", s, sep="")
  )

  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)  
}

