###############################################################
##  Input is the All bystation key-value pairs, stations     ##
##  which does not have observations after 1950 are filtered ##
##  out. Output is a1950 by stations key-value pairs         ##
###############################################################
a1950 <- function() {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      sub <- subset(map.values[[r]], year >= 1950)
      if(sum(is.na(sub$resp)) != 576) {
        attr(sub, "location") <- attributes(map.values[[r]])$location
        rhcollect(map.keys[[r]], sub)
      }
    })
  })
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "All", "bystation"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bystation"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 100,  #cdh3,4
    mapreduce.job.reduces = 100  #cdh5
  )
  job$readback <- FALSE
  job$combiner <- TRUE
  job$jobname <- file.path(rh.root, par$dataset, "a1950", "bystation")
  job.mr <- do.call("rhwatch", job)
  
  
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      v <- map.values[[r]]
      v$station.id <- as.character(map.keys[[r]])
      v$elev <- as.numeric(attributes(v)$location["elev"])
      v$lon <- as.numeric(attributes(v)$location["lon"])
      v$lat <- as.numeric(attributes(v)$location["lat"])
      lapply(1:nrow(v), function(k){
        value <- subset(v[k,], select = -c(year, month))
        rhcollect(c(v$year[k], v$month[k]), value)
      })
    })
  })
  job$reduce <- expression(
    pre = {
      combined <- data.frame()
    },
    reduce = {
      combined <- rbind(combined, do.call(rbind, reduce.values))
    },
    post = {
      row.names(combined) <- NULL
      rhcollect(reduce.key, combined)
    }
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bystation"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bymonth"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 100,  #cdh3,4
    mapreduce.job.reduces = 100  #cdh5
  )
  job$readback <- FALSE
  job$combiner <- TRUE
  job$jobname <- file.path(rh.root, par$dataset, "a1950", "bymonth")
  job.mr <- do.call("rhwatch", job)

}

#############################################################
##  Spatial loess interpolation function apply spaloess    ##
##  function to each month.                                ##
##  Use bymonth.fit as input since bymonth.fit.cv  is only ##
##  used as cross validation purpose                       ##
#############################################################
interpolate <- function(sp, Edeg, deg=2, fam="symmetric", surf="direct") {
  
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      v <- map.values[[r]]
      if(Edeg == 2) {
        v$elev2 <- log2(v$elev + 128)
        lo.fit <- spaloess( resp ~ lon + lat + elev2, 
          data    = v, 
          degree  = degree, 
          span    = span,
          para    = "elev2",
          family  = family,
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = surf),
          napred = TRUE
        )
      } else if(Edeg == 1) {
        v$elev2 <- log2(v$elev + 128)
        lo.fit <- spaloess( resp ~ lon + lat + elev2, 
          data    = v, 
          degree  = degree, 
          span    = span,
          drop    = "elev2",
          para    = "elev2",
          family  = family,
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = surf),
          napred = TRUE
        )
      } else if (Edeg == 0) {
        lo.fit <- spaloess( resp ~ lon + lat, 
          data    = v, 
          degree  = degree, 
          span    = span,
          family  = family,
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = surf),
          napred = TRUE
        )
      }
      value <- subset(v, !is.na(resp))
      if(Edeg == 0) {
        value.na <- subset(v, is.na(resp))
      } else {
        value.na <- subset(v, is.na(resp), select=-c(elev2))
      }
      value$fitted <- lo.fit$fitted
      value.na <- merge(value.na, lo.fit$pred, by= c("lon","lat"))
      value <- rbind(value, value.na)
      if(nrow(value) == nrow(v)) {
        rhcounter("MERGE","CORRECT",1)
      }
      rhcollect(map.keys[[r]], value)
    })
  })
  job$setup <- expression(
    map = {library(Spaloess, lib.loc=lib.loc)}
  )
  job$parameters <- list(
    span = sp,
    degree = deg,
    family = fam,
    Edeg = Edeg,
    surf = surf
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bymonth"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bymonth.fit", fam, surf, Edeg, paste("sp",sp, sep="")), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 20,  #cdh3,4
    mapreduce.job.reduces = 20  #cdh5
  )
  job$readback <- FALSE
  job$combiner <- TRUE
  job$jobname <- file.path(rh.root, par$dataset, "a1950", "bymonth.fit", fam, surf, Edeg, paste("sp",sp, sep=""))
  job.mr <- do.call("rhwatch", job)

}


########################################################
##  The old corss validation function, did not remove ##
##  the point for the fitting (not leave-128-out)     ##
########################################################
crossValid <- function(fam, Edeg, surf, first = FALSE, span) {
  
  FileInput <- paste(
    file.path(rh.root, par$dataset, "a1950", "bymonth.fit", fam, surf, Edeg), "/sp", span, sep=""
  )
  FileOutput <- file.path(rh.root, par$dataset, "a1950", "bymonth.fit", fam, surf, Edeg, "MSE")  

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      v <- map.values[[r]]
      mse <- mean((v$resp-v$fitted)^2, na.rm=TRUE)
      file <- Sys.getenv("mapred.input.file")
      span <- substr(tail(strsplit(file, "/")[[1]],3)[2], 3, 7)
      value <- data.frame(
        span = span, 
        mse = mse, 
        na = sum(!is.na(v$resp-v$fitted)),
        year = map.keys[[r]][1], 
        month = map.keys[[r]][2],
        stringsAsFactors = FALSE
      )
      rhcollect(1, value)
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
  job$input <- rhfmt(FileInput, type = "sequence")
  job$output <- rhfmt(FileOutput, type = "sequence")
  job$mapred <- list(
    mapred.reduce.tasks = 1,  #cdh3,4
    mapreduce.job.reduces = 1  #cdh5 
    #rhipe_reduce_buff_size = 10000
  )
  job$mon.sec <- 10
  job$jobname <- FileOutput  
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)

}


cppkdtree <- function(data, nb) { ## data matrix, no.of.leaves

  D <- ncol(data)  ## number of columns in dm
  ND <- nrow(data)  ## number of rows in dm
  bucketSize <- ND/nb  ## number of rows in each leaf
  
  ## void getkdtree(double *data, int *D, int *ND, int *bucketSize, int* pidx, int *owner)
  res <- .C("getkdtree"
          , as.numeric(data)  ## data matrix as a vector, by column first
          , as.integer(D)  ## number of columns in dm
          , as.integer(ND)  ## number of rows in dm
          , as.integer(bucketSize)  ## number of rows in each leaf
          , idx = integer(ND)  ## pre-allocate index of rows
          , leaf = integer(ND)  ## pre-allocate leaf
  )
  ## output is a list of 6 elements:
  ## 1. dm as a vector
  ## 2. number of columns in dm
  ## 3. number of rows in dm
  ## 4. number of rows in each leaf
  ## 5. index of rows
  ## 6. index of leaves, can be matched up with index of rows to find leaves
  
  ## return(tapply(res$idx, res$leaf, function(r) data[r,], simplify=FALSE))
  ## return(data.frame(idx=res[[5]], leaf=res[[6]]))
  return(res)
  ## return the index of rows and index of leaf

}

##########################################################################
##  Input is the a1950 bymonth division, output is a1950 bymonthSplit.  ##
##  Key is changed from c(year, month) to c(year, month, rep), where    ##
##  rep is the index of observations in a leaf of a kd-tree built       ##
##  based on the non-NA observations in the given month.                ##
##  Input argument leaf control the number of leafs in that kd-tree.    ##
##  The value of the new key-value pairs are data.frame including all   ##
##  non-NA obs in that given month but with a new column named flag     ##
##  which is to indicate for this given month which stations are used   ##
##  as prediction of cross validation.                                  ##
##  In other words, each monthly data.frames have been duplicated       ##
##  multiple times . Totally, there are 22,916 new key-value pairs.     ## 
##########################################################################
bymonthSplit <- function(input, leaf = 100, vari) {
  
  output <- paste(input, "Split", sep="")
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      v <- subset(map.values[[r]], !is.na(map.values[[r]][, vari]))
      rhcounter("Map", "input", 1)
      v$flag <- 0
      row.names(v) <- NULL
      rst <- cppkdtree(as.matrix(v[,c(4,5)]), leaf)
      error <- 0
      kdtree <- data.frame(idx = rst[[5]], leaf = rst[[6]]) ## grab the row index and the leaf index 
      replicates <- ddply(
        .data = kdtree, 
        .vari = "leaf",
        .fun = function(dr) {
          dr$rep <- c(1:nrow(dr))
          dr
        }
      )
      d_ply(
        .data = replicates,
        .vari = "rep",
        .fun = function(ii) {
          rhcounter("Map","keyvalue", 1)
          value <- v
          value$flag[ii$idx] <- 1
          if(sum(value$flag) == 128) {
            rhcounter("Map", "128", 1)
          } else {
            rhcounter("Map","non128", 1)
          }
          rhcollect(c(map.keys[[r]], unique(ii$rep)), value)
        }
      )
    })
  })
  job$setup <- expression(
    map = {
      library(Spaloess, lib.loc=lib.loc)
      library(plyr)
      system("chmod 777 cppkdtree.so")
      dyn.load("cppkdtree.so")
    }
  )
  job$shared <- c(
    #file.path(file.path(rh.root, par$dataset, "shareRLib", "cppkdtree.so"))
    file.path(file.path(rh.root, "shareRLib", "cppkdtree.so"))
  )
  job$parameters <- list(
    leaf = leaf,
    cppkdtree = cppkdtree,
    vari = vari
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(output, type = "sequence")
  job$mapred <- list(
    mapred.tasktimeout = 0,
    mapred.reduce.tasks = 100,  #cdh3,4
    mapreduce.job.reduces = 100  #cdh5
  )
  job$readback <- FALSE
  job$combiner <- TRUE
  job$jobname <- output
  job.mr <- do.call("rhwatch", job)

  return(output)

}

#############################################################################
##  Input is the a1950 bymonthSplit file, for each key-value pairs, 128    ##
##  flagged locations are predicted using the rest of locations. Then      ##
##  the keys are changed from c(year, month, rep) to c(year, month). Value ##
##  is the error of the prediction. In the Reduce, errors of each month    ##
##  is accumulated.                                                        ##
#############################################################################
newCrossValid <- function(input, vari, sp, Edeg, deg=2, fam="symmetric", surf="direct", error = "mse") {

  output <- gsub("Split", ".fit.cv", input)
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      v <- map.values[[r]]
      orig <- subset(v, flag == 1)
      v[v$flag == 1, vari]<- NA
      if (!("elev2" %in% colnames(v))) {
        v$elev2 <- log2(v$elev + 128)
      }
      if (Edeg == 2) {
        fml <- paste(vari, "~ lon + lat + elev2")
        lo.fit <- spaloess( fml, 
          data    = v, 
          degree  = degree, 
          span    = span,
          para    = "elev2",
          family  = family,
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = surf),
          napred = TRUE
        )
      } else if (Edeg == 1) {
        fml <- paste(vari, "~ lon + lat + elev2")
        lo.fit <- spaloess( fml, 
          data    = v, 
          degree  = degree, 
          span    = span,
          drop    = "elev2",
          para    = "elev2",
          family  = family,
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = surf),
          napred = TRUE
        )
      } else if (Edeg == 0) {
        fml <- paste(vari, "~ lon + lat")
        lo.fit <- spaloess( fml, 
          data    = v, 
          degree  = degree, 
          span    = span,
          family  = family,
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = surf),
          napred = TRUE
        )
      }
      value <- merge(orig, lo.fit$pred, by= c("lon","lat"))
      if (vari == "remainder") {
        names(value)[ncol(value)] <- "fitted"
      }
      if(er == "abs") {
        error <- abs(value[, vari] - value[, "fitted"])
      } else if (er == "mse") {
        error <- (value[, vari] - value[, "fitted"])^2
      }
      rhcollect(map.keys[[r]][1:2], error)
    })
  })
  job$reduce <- expression(
    pre = {
      all <- 0
      len <- 0
    },
    reduce = {
      all <- all + sum(unlist(reduce.values))
      len <- len + length(unlist(reduce.values))
    },
    post = {
      rhcollect(reduce.key, all/len)
    }
  )
  job$setup <- expression(
    map = {library(Spaloess, lib.loc=lib.loc)}
  )
  job$parameters <- list(
    vari = vari,
    span = sp,
    degree = deg,
    family = fam,
    Edeg = Edeg,
    surf = surf,
    er = error
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(
    file.path(output, fam, surf, Edeg, paste("sp",sp, sep="")),
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 20, #cdh3,4
    mapreduce.job.reduces = 20 #cdh5
  )
  job$readback <- FALSE
  job$jobname <- file.path(output, fam, surf, Edeg, paste("sp",sp, sep=""))
  job$mon.sec <- 20
  job.mr <- do.call("rhwatch", job)

}

#########################################
##  Read in all sp files for a given fam, surf, and Edeg
##
#########################################
crossValidMerge <- function(input, span) {
  
  Input <- paste(input, "/sp", span, sep="")
  Output <- file.path(input, "MABSE")  

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      v <- map.values[[r]]
      file <- Sys.getenv("mapred.input.file")
      span <- substr(tail(strsplit(file, "/")[[1]],3)[2], 3, 7)
      value <- data.frame(
        span = span, 
        mse = map.values[[r]], 
        year = map.keys[[r]][1], 
        month = map.keys[[r]][2],
        stringsAsFactors = FALSE
      )
      rhcollect(1, value)
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
  job$input <- rhfmt(Input, type = "sequence")
  job$output <- rhfmt(Output, type = "sequence")
  job$mapred <- list(
    mapred.reduce.tasks = 1,  #cdh3,4
    mapreduce.job.reduces = 1  #cdh5 
    #rhipe_reduce_buff_size = 10000
  )
  job$mon.sec <- 20
  job$jobname <- Output  
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)

}

###################################################
##  swap the input key-value pairs from by month ##
##  to by station.id                             ##
###################################################
swapTostation <- function(input, output, elevFlag=TRUE) {

    # job is just changing the key from month and year to station.id from job22
    job <- list()
    job$map <- expression({
      lapply(seq_along(map.values), function(r) {
        map.values[[r]]$year <- map.keys[[r]][1]
        map.values[[r]]$month <- map.keys[[r]][2]
        lapply(1:nrow(map.values[[r]]), function(i){
          value <- map.values[[r]][i, ]
          rhcollect(as.character(value$station.id), value)
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
        lon <- unique(combine$lon)
        lat <- unique(combine$lat)
        elev2 <- unique(combine$elev2)
        if(elevFlag) {
          combine <- subset(combine, select = -c(lon, lat, elev, elev2, station.id))
        } else {
          combine <- subset(combine, select = -c(lon, lat, elev2, station.id))
        }
        attr(combine, "loc") <- c(lon, lat, elev2)
        rhcollect(reduce.key, combine)
      }
    )
    job$parameters <- list(
      elevFlag = elevFlag
    )
    job$combiner <- FALSE
    job$input <- rhfmt(input , type = "sequence")
    job$output <- rhfmt(output, type = "sequence")
    job$mapred <- list(
      mapred.tasktimeout = 0,
      mapred.reduce.tasks = 100,
      rhipe_reduce_buff_size = 10000
    )
    job$mon.sec <- 10
    job$jobname <- output  
    job$readback <- FALSE  

    job.mr <- do.call("rhwatch", job)

}

a1950.STLfit <- function(input, reduce, tuning) {
  
  output <- file.path(
    rh.root, par$dataset, "a1950", "STL", paste("t",tuning$tw, "td", tuning$td, "_s", tuning$sw, 
    "sd", tuning$sd, "_f", tuning$fcw, "fd", tuning$fcd, sep="")
  )

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      loc <- attributes(map.values[[r]])$loc
      value <- arrange(map.values[[r]], year, match(month, month.abb))
      #value$station.id <- map.keys[[r]]
      value$date <- 1:nrow(value)
      Index <- which(is.na(value$resp))
      Resp <- value$resp
      value$flag <- 1
      value$flag[Index] <- 0
      Resp[Index] <- value$fitted[Index]
      if (is.null(par$fcw)) {
        
        fit <- stlplus(
          x=Resp, t=value$date, n.p=12, s.window=par$sw, s.degree=par$sd, 
          t.window=par$tw, t.degree=par$td, inner=10, outer=0
        )$data
        value <- cbind(value, subset(fit, select = c(seasonal, trend, remainder)))

      } else {

        fit <- do.call("cbind", stlplus(
          x=Resp, t=value$date, n.p=12, s.window=par$sw, s.degree=par$sd, t.window=par$tw, 
          t.degree=par$td, fc.window=c(par$tw,par$fcw), fc.degree=c(par$td,par$fcd), inner=10, outer=0
        )[c("data","fc")])
        value <- cbind(value, subset(fit, select = -c(data.raw, data.trend, data.remainder, data.weights, data.sub.labels)))
        names(value)[grep("fc.fc", names(value))] <- c("fc.first", "fc.second")

      }
      attr(value, "loc") <- loc
      rhcollect(map.keys[[r]], value)
    })
  })
  job$reduce <- expression(
    pre = {
      combined <- data.frame()
    },
    reduce = {
      combined <- rbind(combined, do.call("rbind", reduce.values))
    },
    post = { 
      rhcollect(reduce.key, combined)
    }
  )
  job$parameters <- list(
    par = tuning
  )
  job$setup <- expression(
    map = {
      suppressMessages(library(stlplus, lib.loc=lib.loc))
      library(plyr)
    }
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(output, type = "sequence")
  job$mapred <- list(
    mapred.tasktimeout = 0,
    mapred.reduce.tasks = reduce,  #cdh3,4
    mapreduce.job.reduces = reduce  #cdh5
  )
  job$readback <- FALSE
  job$jobname <- output
  job.mr <- do.call("rhwatch", job)

  return(output)

}

swapTomonth <- function(input, output) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      d_ply(
        .data = map.values[[r]],
        .variable = c("year","month"),
        .fun = function(k, station = map.keys[[r]]) {
          key <- c(unique(k$year), unique(as.character(k$month)))
          value <- subset(k, select = -c(year, month))
          value$station.id <- station
          value$lon <- as.numeric(attributes(map.values[[r]])$loc[1])
          value$lat <- as.numeric(attributes(map.values[[r]])$loc[2])
          value$elev2 <- as.numeric(attributes(map.values[[r]])$loc[3])
          rhcollect(key, value)
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
  job$setup <- expression(
    map = {
      library(plyr)
    }
  )
  job$mapred <- list(
    mapred.reduce.tasks = 72,
    mapred.tasktimeout = 0,
    rhipe_reduce_buff_size = 10000
  )
  job$combiner <- TRUE
  job$input <- rhfmt(input, type="sequence")
  job$output <- rhfmt(output, type="sequence")
  job$mon.sec <- 20
  job$jobname <- output
  job$readback <- FALSE  

  job.mr <- do.call("rhwatch", job)  

}

a1950.Spatialfit <- function(input, output, argumt) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      v <- map.values[[r]]
      if(argumt$Edeg == 2) {
        lo.fit <- spaloess( remainder ~ lon + lat + elev2, 
          data    = v, 
          degree  = argumt$degree, 
          span    = argumt$span,
          para    = "elev2",
          family  = "symmetric",
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = argumt$surf),
          napred = FALSE
        )
      } else if(argumt$Edeg == 1) {
        lo.fit <- spaloess( remainder ~ lon + lat + elev2, 
          data    = v, 
          degree  = argumt$degree, 
          span    = argumt$span,
          drop    = "elev2",
          para    = "elev2",
          family  = "symmetric",
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = argumt$surf),
          napred = FALSE
        )
      } else if (argumt$Edeg == 0) {
        lo.fit <- spaloess( remainder ~ lon + lat, 
          data    = v, 
          degree  = argumt$degree, 
          span    = argumt$span,
          family  = "symmetric",
          normalize = FALSE,
          distance = "Latlong",
          control = loess.control(surface = argumt$surf),
          napred = FALSE
        )
      }
      v$spafit <- lo.fit$fitted
      rhcollect(map.keys[[r]], v)
    })
  })
  job$parameters <- list(
    argumt = argumt
  )
  job$setup <- expression(
    map = {
      library(Spaloess, lib.loc=lib.loc)
    }
  )
  job$mapred <- list(
    mapred.reduce.tasks = 72,
    mapred.tasktimeout = 0
  )
  job$input <- rhfmt(input, type="sequence")
  job$output <- rhfmt(output, type="sequence")
  job$mon.sec <- 10
  job$jobname <- output
  job$readback <- FALSE  

  job.mr <- do.call("rhwatch", job)  

}


constructQuants <- function(obj, probs, tails, mids) {
  keys <- lapply(obj, function(x) x[[1]][[2]])
  intKeys <- sapply(keys, is.integer)
  vals <- lapply(obj, "[[", 2)

  quants <- data.frame(
    idx = unlist(keys[intKeys]),
    freq = unlist(vals[intKeys])
  )

  tot <- sum(as.numeric(quants$freq))
  quants <- quants[order(quants$idx),]
  quants$pct <- quants$freq / tot
  quants$cpct <- cumsum(quants$pct)
  quants$q <- mids[quants$idx]

  fn <- approxfun(quants$cpct, quants$q, method = "constant", f = 1, rule = 2)
  res <- data.frame(
    fval = probs,
    q = fn(probs)
  )

  if(tails > 0) {
    # now append top and bottom
    tailKeys <- unlist(keys[!intKeys])
    tailVals <- vals[!intKeys]

    top <- tailVals[tailKeys == "top"][[1]]
    bot <- tailVals[tailKeys == "bot"][[1]]

    botDf <- data.frame(
      fval = (seq_len(tails) - 1) / tot,
      q = bot
    )

    topDf <- data.frame(
      fval = (seq_len(tails) + tot - tails) / tot,
      q = top
    )

    res <- res[res$fval > max(botDf$fval) & res$fval < min(topDf$fval),]
    res <- rbind(botDf, res, topDf)
  }
  res
}

coastDist <- function(data, lim) {

  pts <- as.matrix(data[,c("lon","lat")], ncol = 2)
  mp <- map("usa", plot = FALSE)
  data$coastdis <- NA
  xy.coast <- cbind(mp$x, mp$y)[!is.na(mp$x), ]
  for (i in 1:nrow(pts)) {
    data$coastdis[i] <- as.numeric(min(spDistsN1(xy.coast, pts[i,], longlat = TRUE)) <= lim)
  }
  
  return(data)

}

getCondCuts <- function(df, splitVars) {
  apply(
    X = do.call("cbind", lapply(df[,splitVars,drop = FALSE], function(x) format(x, scientific = FALSE, trim = TRUE, justify = "none"))), 
    MARGIN = 1,
    FUN = function(x) paste(paste(splitVars, "=", x, sep = ""), collapse = "|")
  )
}

a1950.residQuant <- function(input, target="residual", by=NULL, probs=seq(0, 1, 0.005), nBins = 10000, tails = 100, coast=FALSE, dislim=NULL) {

  ## Get the range of residual
  jobRng <- list()
  jobRng$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      v <- subset(map.values[[r]], flag==1)
      if(target == "residual") (
        value <- with(v, remainder - spafit)
      ) else {
        value <- v[, target]
      }
      rhcollect(1, data.frame(bot=min(value), top=max(value)))
    })
  })
  jobRng$reduce <- expression(
    pre = {
      combine <- data.frame()
    },
    reduce = {
      combine <- rbind(combine, do.call("rbind", reduce.values))
    },
    post = {
      rhcollect(reduce.key, data.frame(bot=min(combine$bot), top=max(combine$top)))
    }
  )
  jobRng$mapred <- list(
    mapred.reduce.tasks = 1,
    mapred.tasktimeout = 0
  )
  jobRng$parameters <- list(
    target = target
  )
  jobRng$input <- rhfmt(input, type="sequence")
  jobRng$output <- rhfmt(paste(input, target, "range", sep="."), type="sequence")
  jobRng$mon.sec <- 10
  jobRng$jobname <- paste(input, target, "range", sep=".")
  jobRng$readback <- TRUE
  varRange <- as.numeric(do.call("rhwatch", jobRng)[[1]][[2]])

  delta <- abs(diff(varRange)) / (nBins - 1)
  cuts <- seq(varRange[1] - delta / 2, varRange[2] + delta / 2, by = delta)
  mids <- seq(varRange[1], varRange[2], by = delta)

  jobQuant <- list()
  jobQuant$map <- expression({
    dat <- do.call("rbind", lapply(seq_along(map.keys), function(r) {
      value <- map.values[[r]]
      if (coast) {
        value <- coastDist(value, dislim)
      }
      data.frame(
        v = with(subset(value, flag==1), remainder - spafit),
        subset(value, flag==1)[, by, drop = FALSE],
        stringsAsFactors = FALSE
      )
    }))
    if(is.null(by)) {
      inds <- list("1" = seq_len(nrow(dat)))
    } else {
      splits <- getCondCuts(dat[, by, drop = FALSE], by)
      inds <- split(seq_along(splits), splits)
    }
    ## inds is a list of vector, each element of list is a vector of row index for the given level
    ## of conditional variable by. Next we loop over all different levels of by variable
    ## If by is NULL, then inds length is 1.
    for (ii in seq_along(inds)) {
      vv <- dat$v[inds[[ii]]]  ## get the target variable for the given level of by variable
      vv <- vv[!is.na(vv)]
      if (length(vv) > 0) {
        ord <- order(vv)
        cutTab <- as.data.frame(
          table(cut(vv, cuts, labels = FALSE)), responseName = "Freq", stringsAsFactors = FALSE
        )
        cutTab$Var1 <- as.integer(cutTab$Var1)
        ## for each bin we collect the key-value pair
        ## key is the bin index cutTab$Var1, and value is the frequency count in that bin
        for (jj in 1:nrow(cutTab)) {
          rhcollect(list(as.list(dat[inds[[ii]][1], by, drop = FALSE]), cutTab$Var1[jj]), cutTab$Freq[jj])
        } 
        rhcollect(list(as.list(dat[inds[[ii]][1], by, drop = FALSE]), "bot"), vv[head(ord, tails)])
        rhcollect(list(as.list(dat[inds[[ii]][1], by, drop = FALSE]), "top"), vv[tail(ord, tails)])
      }
    }    
  })
  jobQuant$reduce <- expression(
    pre = {
      bot <- NULL
      top <- NULL
      sum <- 0
    }, 
    reduce = {
      if(reduce.key[[2]] == "bot") {
        bot <- head(sort(c(bot, do.call(c, reduce.values))), tails)
      } else if(reduce.key[[2]] == "top") {
        top <- tail(sort(c(top, do.call(c, reduce.values))), tails)
      } else {
        sum <- sum + sum(unlist(reduce.values))
      }
    }, 
    post = {
      if(reduce.key[[2]] == "bot") {
        rhcollect(reduce.key, bot)
      } else if(reduce.key[[2]] == "top") {
        rhcollect(reduce.key, top)
      } else {
        rhcollect(reduce.key, sum)
      }
    }
  )
  jobQuant$parameters <- list(
    by = by,
    delta = delta,
    cuts = cuts,
    mids = mids,
    tails = tails,
    coast = coast,
    coastDist = coastDist,
    getCondCuts = getCondCuts,
    dislim = dislim
  )
  jobQuant$setup <- expression(
    map = {
      suppressMessages(library(sp, lib.loc=lib.loc))
      library(maps, lib.loc=lib.loc)
    }
  )
  jobQuant$mapred <- list(
    mapred.reduce.tasks = 50,
    mapred.tasktimeout = 0
  # rhipe_map_buff_size = 5
  )
  jobQuant$input <- rhfmt(input, type="sequence")
  jobQuant$output <- rhfmt(paste(input, target, "quant", sep="."), type="sequence")
  jobQuant$mon.sec <- 10
  jobQuant$jobname <- paste(input, target, "quant", sep=".")
  jobQuant$readback <- TRUE
  mrRes <- do.call("rhwatch", jobQuant)
  
  if(is.null(by)) {
    res <- constructQuants(mrRes, probs, tails, mids)
  } else {
    groups <- sapply(mrRes, function(x) {
      do.call(paste, c(as.list(x[[1]][[1]]), sep = "|"))
    })
    ind <- split(seq_along(groups), groups)

    res <- lapply(seq_along(ind), function(i) {
      data.frame(
        constructQuants(mrRes[ind[[i]]], probs, tails, mids),
        mrRes[[ind[[i]][1]]][[1]][[1]],
        stringsAsFactors = FALSE
      )
    })
    res <- do.call("rbind", res)
  }

  res

}

a1950.Nomiss <- function(input) {

  output <- paste(input, "noMissStations", sep=".") 

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      if (sum(is.na(map.values[[r]]$resp)) == 0) {
        rhcollect(1, map.keys[[r]])
      }
    })
  })
  job$reduce <- expression(
    pre = {
      combine <- character()
    },
    reduce = {
      combine <- c(combine, unlist(reduce.values))
    },
    post = {
      rhcollect(reduce.key, combine)
    }
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(output, type = "sequence")
  job$mapred <- list(
    mapred.reduce.tasks = 1,  #cdh3,4
    mapreduce.job.reduces = 1  #cdh5
  )
  job$readback <- TRUE
  job$jobname <- output
  job.mr <- do.call("rhwatch", job)  

  a1950Nomiss <- job.mr[[1]][[2]]
  rhsave(a1950Nomiss, file=file.path(rh.root, par$dataset, "a1950","Rdata","nomiss.a1950.RData"))

}
