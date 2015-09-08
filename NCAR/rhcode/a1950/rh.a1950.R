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
          lo.fit <- my.loess2( resp ~ lon + lat + elev2, 
            data    = v, 
            degree  = degree, 
            span    = span,
            para    = "elev2",
            family  = family,
            control = loess.control(surface = surf)
          )
        } else if(Edeg == 1) {
          v$elev2 <- log2(v$elev + 128)
          lo.fit <- my.loess2( resp ~ lon + lat + elev2, 
            data    = v, 
            degree  = degree, 
            span    = span,
            drop    = "elev2",
            para    = "elev2",
            family  = family,
            control = loess.control(surface = surf)
          )
        }
        fit <- my.predict.loess(
          object = lo.fit, 
          newdata = data.frame(
            lon = v$lon, 
            lat = v$lat,
            elev2 = v$elev2
          )
        ) 

      } else {

        lo.fit <- my.loess2( resp ~ lon + lat, 
          data    = v, 
          degree  = degree, 
          span    = span,
          family  = family
        )
        fit <- my.predict.loess(
          object = lo.fit, 
          newdata = data.frame(
            lon = v$lon, 
            lat = v$lat
          )
        )

      }
      v$fitted <- fit
      rhcollect(map.keys[[r]], v)
    })
  })
  job$setup <- expression(
    map = {
      system("chmod 777 myloess2.so")
      dyn.load("myloess2.so")
    }
  )
  job$shared <- c(
    file.path(file.path(rh.root, par$dataset, "shareRLib", "myloess2.so"))
  )
  job$parameters <- list(
    Elev = Elev,
    span = sp,
    degree = deg,
    family = fam,
    Edeg = Edeg,
    surf = surf,
    my.loess2 = my.loess2,
    my.simple2 = my.simple2,
    my.predict.loess = my.predict.loess,
    my.predLoess = my.predLoess
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

crossValid <- function(fam, Edeg, surf, first = FALSE) {
  
  FileInput <- file.path(rh.root, par$dataset, "a1950", "bymonth.fit", fam, surf, Edeg)
  FileOutput <- file.path(rh.root, par$dataset, "a1950", "bymonth.fit", fam, surf, Edeg, "MSE")  
  
  if(!first) {
    rhdel(FileOutput)
  }

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
