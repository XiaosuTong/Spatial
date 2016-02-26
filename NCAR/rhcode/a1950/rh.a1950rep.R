repTime <- function(input, output, buffSize, Rep=9600){

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      Index <- which(is.na(map.values[[r]]$resp))
      map.values[[r]][Index, "resp"] <- map.values[[r]]$fitted[Index]
      value <- subset(arrange(map.values[[r]], year, match(month, month.abb)), select = -c(fitted, year))
      value <- rdply(Rep, value, .id=NULL)
      value$date <- 1:nrow(value)
      row.names(value) <- NULL
      attr(value, "loc") <- attributes(map.values[[r]])$loc
      rhcollect(map.keys[[r]], value)
    })
  })
  job$setup <- expression(
    map = {
      suppressMessages(library(plyr, lib.loc=lib.loc))
    }
  )
  job$parameters <- list(
    Rep = Rep
  )
  job$input <- rhfmt(input, type = "sequence")
  job$output <- rhfmt(output, type = "sequence")
  job$mapred <- list(
    mapred.reduce.tasks = 0,  #cdh3,4
    mapreduce.job.reduces = 0,  #cdh5
    rhipe_map_buff_size = buffSize
    #rhipe_reduce_buff_size = 10000
  )
  job$mon.sec <- 20
  job$jobname <- output  
  job$readback <- FALSE
  job.mr <- do.call("rhwatch", job)

}

swapTomonth <- function(input, output, io_sort, spill_percent, parallelcopies, reduce_merge_inmem, reduce_input_buffer) {

  job <- list()
  job$map <- expression({
    lapply(seq_along(map.values), function(r) {
      for(i in 1:nrow(map.values[[r]])) {
        key <- c(map.values[[r]]$date[i], as.character(map.values[[r]]$month[i]))
        value <- data.frame(
          station.id = map.keys[[r]],
          lon = as.numeric(attributes(map.values[[r]])$loc[1]),
          lat = as.numeric(attributes(map.values[[r]])$loc[2]),
          elev2 = as.numeric(attributes(map.values[[r]])$loc[3]),
          resp = map.values[[r]]$resp[i],
          stringsAsFactors = FALSE
        )
        rhcollect(key, value)
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
  job$setup <- expression(
    map = {library(plyr, lib.loc=lib.loc)}
  )
  job$mapred <- list(
    mapred.reduce.tasks = 100,  #cdh3,4
    mapreduce.job.reduces = 100,  #cdh5
    mapreduce.task.io.sort.mb = io_sort,
    mapreduce.map.sort.spill.percent = spill_percent,
    mapreduce.reduce.shuffle.parallelcopies = parallelcopies,
    mapreduce.reduce.merge.inmem.threshold = reduce_merge_inmem,
    mapreduce.reduce.input.buffer.percent = reduce_input_buffer,
    mapred.tasktimeout = 0,
    rhipe_reduce_buff_size = 10000
  )
  job$combiner <- TRUE
  job$input <- rhfmt(input, type="sequence")
  job$output <- rhfmt(output, type="sequence")
  job$mon.sec <- 10
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
        fml <- remainder ~ lon + lat + elev2
        dropSq <- FALSE
        condParam <- "elev2"
      } else if(argumt$Edeg == 1) {
        fml <- remainder ~ lon + lat + elev2
        dropSq <- "elev2"
        condParam <- "elev2"
      } else if (argumt$Edeg == 0) {
        fml <- remainder ~ lon + lat
        dropSq <- FALSE
        condParam <- FALSE
      }
      lo.fit <- spaloess( fml, 
        data    = v, 
        degree  = argumt$degree, 
        span    = argumt$span,
        para    = condParam,
        drop    = dropSq,
        family  = "symmetric",
        normalize = FALSE,
        distance = "Latlong",
        control = loess.control(surface = argumt$surf),
        napred = FALSE
      )
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
    mapred.reduce.tasks = 100,  #cdh3,4
    mapreduce.job.reduces = 100,  #cdh5
    mapred.tasktimeout = 0
  )
  job$input <- rhfmt(input, type="sequence")
  job$output <- rhfmt(output, type="sequence")
  job$mon.sec <- 20
  job$jobname <- output
  job$readback <- FALSE  

  job.mr <- do.call("rhwatch", job)  

}

swapTostation <- function(input, output, elevFlag=FALSE) {

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
      mapred.reduce.tasks = 72
    )
    job$mon.sec <- 10
    job$jobname <- output  
    job$readback <- FALSE  

    job.mr <- do.call("rhwatch", job)

}

a1950.STLfit <- function(input, output, reduce, tuning) {
  
  job <- list()
  job$map <- expression({
    lapply(seq_along(map.keys), function(r) {
      loc <- attributes(map.values[[r]])$loc
      value <- arrange(map.values[[r]], date)
      #value$station.id <- map.keys[[r]]
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
      library(dplyr, lib.loc=lib.loc)
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

}

a1950OptFit <- function(inSpan, inDegree, inEdeg, inSurf, sw, sd, tw, td, span, degree, Edeg, surf) {
  
  FileInput <- file.path(rh.root, par$dataset, "a1950", "bystation")
  FileOutput <- file.path(rh.root, par$dataset, "a1950Rep", "bystation")
  repTime(FileInput, FileOutput, Rep=8000, buffSize=7)

  FileInput <- FileOutput
  FileOutput <- file.path(rh.root, par$dataset, "a1950Rep", "bymonth")
  swapTomonth(FileInput, FileOutput, elevFlag=TRUE)

  paraImput <- list(span=inSpan, degree=inDegree, Edeg=inEdeg, surf=inSurf)
  FileInput <- FileOutput
  FileOutput <- file.path(rh.root, par$dataset, "a1950Rep", "spaimpute.bymonth")
  a1950.Spatialfit(FileInput, FileOutput, paraImput)

  FileInput <- FileOutput
  FileOutput <- file.path(rh.root, par$dataset, "a1950Rep", "spaimpute.bystation")
  swapTostation(FileInput, FileOutput, elevFlag=FALSE)
  
  paraSTL <- list(sw=sw, sd=sd, tw=tw, td=td, fcw=NULL, fcd=NULL)
  FileInput <- FileOutput
  FileOutput <- file.path(rh.root, par$dataset, "a1950Rep", "stlplus.bystation")
  a1950.STLfit(FileInput, FileOutput, paraSTL)

  FileInput <- FileOutput
  FileOutput <- file.path(rh.root, par$dataset, "a1950Rep", "stlplus.bymonth")
  swapTomonth(FileInput, FileOutput, elevFlag=FALSE)

  parafit <- list(span=span, degree=degree, Edeg=Edeg, surf=surf)
  FileInput <- FileOutput
  FileOutput <- file.path(rh.root, par$dataset, "a1950Rep", "spafit.bymonth")
  a1950.Spatialfit(FileInput, FileOutput, parafit)

  FileInput <- FileOutput
  FileOutput <- file.path(rh.root, par$dataset, "a1950Rep", "spafit.bystation")
  swapTostation(FileInput, FileOutput, elevFlag=FALSE)

}