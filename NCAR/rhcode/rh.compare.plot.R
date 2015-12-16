fitcompare <- function(type = "a1950", mod1, mod2, summary=FALSE) {
  ## If summary is TRUE, we will calculate a summary statistic for each station
  ## and combine all station to one data.frame
  if(summary) {

    job <- list()
    job$map <- expression({
      lapply(seq_along(map.keys), function(r) {  

      })
    })
    job$reduce <- expression(
      pre = {  

      },
      reduce = {  

      },
      post = {  

      }
    )
    files <- c(
      file.path(rh.root, par$dataset, type, "STL", mod1),
      file.path(rh.root, par$dataset, type, "STL", mod2)
    )
    job$input <- rhfmt(files, type = "sequence")
    job$output <- rhfmt(
      file.path(rh.root, par$dataset, type, "STLcompare", paste(mod1, mod2, sep="VS")),
      type = "sequence"
    )
    job$mapred <- list(
      mapreduce.job.reduces = 1,  #cdh5
      mapred.reduce.tasks = 1
    )
    job$jobname <- file.path(rh.root, par$dataset, type, "STLcompare", paste(mod1, mod2, sep="VS"))
    job$combiner <- TRUE
    job$readback <- FALSE
    job$mon.sec <- 10
    job.mr <- do.call("rhwatch", job)

  } else { ## else if summary is FALSE, only keep 128 stations and all residual of each 128 stations to one data.frame

    job <- list()
    job$map <- expression({
      lapply(seq_along(map.keys), function(r) {  
      	if(map.keys[[r]] %in% sample.a1950$station.id) {
          value <- map.values[[r]][, "remainder", drop=FALSE]
          value$station.id <- map.keys[[r]]
          file <- Sys.getenv("mapred.input.file")
          span <- substr(tail(strsplit(file, "/")[[1]],3)[2], 3, 7)
          value$group <- 
          rhcollect(1, value)
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
    job$shared <- file.path(rh.datadir, par$dataset, type, "Rdata", "sample.a1950.RData")
    job$setup <- expression(
      map = {
      	library(plyr)
      	load("sample.a1950.RData")
      }
    )
    job$parameters <- list(
      mod1 = mod1,
      mod2 = mod2
    )
    files <- c(
      file.path(rh.root, par$dataset, type, "STL", mod1),
      file.path(rh.root, par$dataset, type, "STL", mod2)
    )
    job$input <- rhfmt(files, type = "sequence")
    job$output <- rhfmt(
      file.path(rh.root, par$dataset, type, "STLcompare", paste(mod1, mod2, sep="VS")),
      type = "sequence"
    )
    job$mapred <- list(
      mapreduce.job.reduces = 1,  #cdh5
      mapred.reduce.tasks = 1
    )
    job$jobname <- file.path(rh.root, par$dataset, type, "STLcompare", paste(mod1, mod2, sep="VS"))
    job$combiner <- TRUE
    job$readback <- FALSE
    job$mon.sec <- 10
    job.mr <- do.call("rhwatch", job)

  }
}