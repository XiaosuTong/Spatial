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

interpolate <- function(Elev = TRUE, sp, Edeg, deg=2, fam="symmetric", surf="direct") {
  
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      v <- map.values[[r]]
      if(Elev) {
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
        }
      } else {
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
      value <- merge(v, lo.fit$pred, by=c("lon","lat"))
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
    Elev = Elev,
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
    file.path(rh.root, par$dataset, "a1950", "bymonth.fit.new", fam, surf, Edeg, paste("sp",sp, sep="")), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 20,  #cdh3,4
    mapreduce.job.reduces = 20  #cdh5
  )
  job$readback <- FALSE
  job$combiner <- TRUE
  job$jobname <- file.path(rh.root, par$dataset, "a1950", "bymonth.fit.new", fam, surf, Edeg, paste("sp",sp, sep=""))
  job.mr <- do.call("rhwatch", job)

}


interpolateStation <- function( sp, Edeg, deg=2, fam="symmetric", surf="direct"){

  FileInput <- file.path(
    rh.root, par$dataset, "a1950", "bymonth.fit", fam, surf, Edeg, paste("sp",sp, sep="")
  )
  FileOutput <- file.path(
    rh.root, par$dataset, "a1950", "bymonth.fit", fam, surf, Edeg, paste("sp",sp, ".bystation", sep="")
  )
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
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
      rhcollect(reduce.key, combine)
    }
  )
  job$combiner <- TRUE
  job$input <- rhfmt(FileInput , type = "sequence")
  job$output <- rhfmt(FileOutput, type = "sequence")
  job$mapred <- list(mapred.reduce.tasks = 100)
  job$mon.sec <- 10
  job$jobname <- FileOutput
  job$readback <- FALSE  

  job.mr <- do.call("rhwatch", job)
  
  return(FileOutput)

} 




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

a1950.STLfit <- function(input, reduce, sw, sd, tw, td, fcw=NULL, fcd=NULL) {

  tuning <- list(sw=sw, sd=sd, tw=tw, td=td, fcw=fcw, fcd=fcd)

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      value <- arrange(map.values[[r]], year, match(month, month.abb))
      value$station.id <- map.keys[[r]]
      value$date <- 1:nrow(value)
      Index <- which(is.na(value$resp))
      Resp <- value$resp
      value$flag <- 1
      value$flag[Index] <- 0
      Resp[Index] <- value$fitted[Index]
      if (is.null(par$fcw)) {
        
        fit <- stl2(
          x=Resp, t=value$date, n.p=12, s.window=par$sw, s.degree=par$sd, 
          t.window=par$tw, t.degree=par$td, inner=10, outer=0
        )$data
        value <- cbind(value, subset(fit, select = c(seasonal, trend, remainder)))

      } else {

        fit <- do.call("cbind", stl2(
          x=Resp, t=value$date, n.p=12, s.window=par$sw, s.degree=par$sd, t.window=par$tw, 
          t.degree=par$td, fc.window=c(par$tw,par$fcw), fc.degree=c(par$td,par$fcd), inner=10, outer=0
        )[c("data","fc")])
        value <- cbind(value, subset(fit, select = -c(data.raw, data.trend, data.remainder, data.weights, data.sub.labels)))
        names(value)[grep("fc.fc", names(value))] <- c("fc.first", "fc.second")

      }
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
      suppressMessages(library(stl2))
      library(plyr)
    }
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(
    file.path(
      rh.root, par$dataset, "a1950", "STL", paste("t",tuning$tw, "td", tuning$td, "_s", tuning$sw, 
      "sd", tuning$sd, "_f", tuning$fcw, "fd", tuning$fcd, sep="")
    ), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = reduce,  #cdh3,4
    mapreduce.job.reduces = reduce  #cdh5
  )
  job$readback <- FALSE
  job$jobname <- file.path(
    rh.root, par$dataset, "a1950", "STL", paste("t",tuning$tw, "td", tuning$td, "_s", 
    tuning$sw, "sd", tuning$sd, "_f", tuning$fcw, "fd", tuning$fcd, sep="")
  )
  
  job.mr <- do.call("rhwatch", job)

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r){
      if(map.keys[[r]] %in% sample.a1950$station.id) {
        rhcollect(1, map.values[[r]])
      }
    })
  })
  job$reduce <- expression(
    pre = {
      combine <- data.frame()
    },
    reduce = {
      combine <- rbind(combine, do.call("rbind", reduce.values))
    },
    post = {
      rhcollect(reduce.key, combine)
    }
  )
  job$shared <- file.path(rh.root, par$dataset, "a1950", "Rdata", "sample.a1950.RData")
  job$setup <- expression(
    map = {load("sample.a1950.RData")}
  )
  job$input <- rhfmt(
    file.path(
      rh.root, par$dataset, "a1950", "STL", paste("t",tuning$tw, "td", tuning$td, "_s", tuning$sw, 
      "sd", tuning$sd, "_f", tuning$fcw, "fd", tuning$fcd, sep="")
    ), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(
      rh.root, par$dataset, "a1950", "STL.plot", paste("t",tuning$tw, "td", tuning$td, "_s", tuning$sw, 
      "sd", tuning$sd, "_f", tuning$fcw, "fd", tuning$fcd, sep="")
    ), 
    type = "sequence"
  )

  job$mapred <- list(
    mapred.reduce.tasks = 1,  #cdh3,4
    mapreduce.job.reduces = 1  #cdh5
  )
  job$readback <- TRUE
  job$jobname <- file.path(
    rh.root, par$dataset, "a1950", "STL.plot", paste("t",tuning$tw, "td", tuning$td, "_s", 
    tuning$sw, "sd", tuning$sd, "_f", tuning$fcw, "fd", tuning$fcd, sep="")
  )
  
  job.mr <- do.call("rhwatch", job)

  return(job.mr)
  
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
bymonthSplit <- function(leaf = 100) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      v <- subset(map.values[[r]], !is.na(resp))
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
          if(nrow(value) == nrow(v)){
            rhcounter("Map","nrow",1)
          }
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
    file.path(file.path(rh.root, par$dataset, "shareRLib", "cppkdtree.so"))
  )
  job$parameters <- list(
    leaf = leaf,
    cppkdtree = cppkdtree
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bymonth"), 
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bymonthSplit"), 
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 100,  #cdh3,4
    mapreduce.job.reduces = 100  #cdh5
  )
  job$readback <- FALSE
  job$combiner <- TRUE
  job$jobname <- file.path(rh.root, par$dataset, "a1950", "bymonthSplit")
  job.mr <- do.call("rhwatch", job)

}

#############################################################################
##  Input is the a1950 bymonthSplit file, for each key-value pairs, 128    ##
##  flagged locations are predicted using the rest of locations. Then      ##
##  the keys are changed from c(year, month, rep) to c(year, month). Value ##
##  is the error of the prediction. In the Reduce, errors of each month    ##
##  is accumulated.                                                        ##
#############################################################################
newCrossValid <- function(Elev = TRUE, sp, Edeg, deg=2, fam="symmetric", surf="direct", error = "mse") {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      v <- map.values[[r]]
      orig <- subset(v, flag == 1)
      v$resp[v$flag == 1]<- NA
      if(Elev) {
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
        }
      } else {
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
      value <- merge(orig, lo.fit$pred, by= c("lon","lat"))
      if(er == "abs") {
        error <- with(value, abs(resp - fitted))
      } else if (er == "mse") {
        error <- with(value, (resp - fitted)^2)
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
    Elev = Elev,
    span = sp,
    degree = deg,
    family = fam,
    Edeg = Edeg,
    surf = surf,
    er = error
  )
  job$input <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bymonthSplit"),
    type = "sequence"
  )
  job$output <- rhfmt(
    file.path(rh.root, par$dataset, "a1950", "bymonth.fit.new", fam, surf, Edeg, paste("sp",sp, sep="")),
    type = "sequence"
  )
  job$mapred <- list(
    mapred.reduce.tasks = 20, #cdh3,4
    mapreduce.job.reduces = 20 #cdh5
  )
  job$readback <- FALSE
  job$jobname <- file.path(
    rh.root, par$dataset, "a1950", "bymonth.fit.new", fam, surf, Edeg, paste("sp",sp, sep="")
  )
  job$mon.sec <- 20
  job.mr <- do.call("rhwatch", job)

}

#########################################
##  Read in all sp files for a given fam, surf, and Edeg
##
#########################################
crossValidMerge <- function(fam, Edeg, surf, first = FALSE, span) {
  
  FileInput <- paste(
    file.path(rh.root, par$dataset, "a1950", "bymonth.fit.new", fam, surf, Edeg), "/sp", span, sep=""
  )
  FileOutput <- file.path(rh.root, par$dataset, "a1950", "bymonth.fit.new", fam, surf, Edeg, "MABSE")  

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
  job$input <- rhfmt(FileInput, type = "sequence")
  job$output <- rhfmt(FileOutput, type = "sequence")
  job$mapred <- list(
    mapred.reduce.tasks = 1,  #cdh3,4
    mapreduce.job.reduces = 1  #cdh5 
    #rhipe_reduce_buff_size = 10000
  )
  job$mon.sec <- 20
  job$jobname <- FileOutput  
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)

}
