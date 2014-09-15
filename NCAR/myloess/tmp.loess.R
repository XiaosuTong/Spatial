source("my.loess02.R")

set.seed(10)
df <- data.frame(
	x = runif(100, -125, -65), 
	y = runif(100, 20, 50), 
	z = rnorm(100)
)

lo.fit <- my.loess2(z~x+y, df)