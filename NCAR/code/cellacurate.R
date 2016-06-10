cellacurate <- function(cell = 0.2, run = TRUE) {

  if (run) {
    mcontrol <- spacetime.control(
      vari="resp", time="date", n=576, n.p=12, stat_n=7738, mthbytime = 1, surf = "interpolate",
      s.window="periodic", t.window = 241, degree=2, span=0.015, Edeg=2, cell = cell
    )
    ccontrol <- mapreduce.control(
      libLoc= .libPaths(), reduceTask=95, io_sort=128, slow_starts = 0.5,
      map_jvm = "-Xmx200m", reduce_jvm = "-Xmx200m",
      map_memory = 1024, reduce_memory = 1024,
      reduce_input_buffer_percent=0.4, reduce_parallelcopies=10,
      reduce_merge_inmem=0, task_io_sort_factor=100,
      spill_percent=0.9, reduce_shuffle_input_buffer_percent = 0.8,
      reduce_shuffle_merge_percent = 0.4
    )
    drsstl(
      data = "/user/tongx/spatem/tmax.txt", output = paste0("/user/tongx/spatem/output", cell),
      stat_info = "/user/tongx/spatem/station_info.RData", model_control = mcontrol,
      cluster_control=ccontrol
    )
  }
  tmp <- rhread(paste0("/user/tongx/spatem/output", cell, "/output_bystat"))
  
  rst <- do.call("rbind", lapply(tmp, "[[", 2))
  rst$station.id <- rep(unlist(lapply(tmp, "[[", 1)), each=576)

  rst <- arrange(rst, station.id, date)

  rst$tmax <- tmax_all$tmax

  return(sum(with(rst, tmax - trend - seasonal - Rspa)^2, na.rm = TRUE)/nrow(rst))

}