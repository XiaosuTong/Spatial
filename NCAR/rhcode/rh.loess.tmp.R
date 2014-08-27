source("~/Rhipe/rhinitial.R")
#source("~/Projects/Spatial/NCAR/rhcode/rh.setup.R")
source("~/Projects/Spatial/NCAR/rhcode/my.loess02.R")

par <- list()
par$N <- 10
par$span <- 0.2
par$degree <- 2
job <- list()
job$map <- expression({
	lapply(seq_along(map.values), function(r) {
		set.seed(99)
		x = rnorm(100)
		set.seed(100)
		y = rnorm(100)
		set.seed(map.keys[[r]])
		z = rnorm(100)
		df = data.frame(
			x = x,
			y = y,
			z = z
		)
		lo.fit <- my.loess2(z ~ x + y,
			data = df,
			span = par$span,
			degree = par$degree
		)
		rhcollect(map.keys[[r]], lo.fit$kd$vert2)
	})
})
job$reduce <- expression(
	pre = {
	},
	reduce = {
	},
	post = {
		rhcollect(reduce.key, reduce.values)
	}	
)
job$setup <- expression(
	map = {
		dyn.load("myloess2.so")
	}
)
job$shared <- c(
	"/ln/tongx/myloess/myloess2.so",
)
job$parameters <- list(
	par = par, 
	my.loess = my.loess2, 
	my.simple = my.simple2
)
job$input <- c(par$N, 10) 
job$output <- rhfmt(
	"/ln/tongx/userhipe/loess", 
	type = "sequence"
)
job$mapred <- list(
	mapred.reduce.tasks = 1, 
	rhipe_reduce_buff_size = 10000
)
job$mon.sec <- 5
job$jobname <- "/ln/tongx/userhipe/loess"
job$readback <- FALSE
job.mr <- do.call("rhwatch", job)

