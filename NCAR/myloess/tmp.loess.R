set.seed(100)
df <- data.frame(x = rnorm(100), y = rnorm(100), z = rnorm(100), w = rnorm(100))
newx <- data.frame(x = runif(10), y = runif(10))
lo.fit1 <- my.loess1(
	z ~ x + y, 
	data = df, 
	span = 0.5,
	degree = 1, 
	normalize = FALSE, 
	control = loess.control(surface = "direct")
)
lo.fit2 <- my.loess2(
	z ~ x + y, 
	data = df, 
	span = 0.5,
	degree = 1, 
	normalize = FALSE
)

tmp <- data.frame(matrix(lo.fit2$kd$vval, byrow=TRUE, ncol=3))
tmp <- cbind(tmp, lo.fit2$kd$vert2)
tmp <- setNames(tmp, c("fitted", "b1", "b2", "x","y"))



#kd <- lo.fit2$kd$vert2
#value <-predict(lo.fit1, data.frame(x = kd$X1, y = kd$X2))
#tmp <- data.frame(matrix(lo.fit2$kd$vval, byrow=TRUE, ncol=3))
#names(tmp) <- c("b0", "b1", "b2")
#fit <- cbind(tmp, kd)
#fit$fitted <- with(fit, b1*X1 + b2*X2)
#lo.fit3 <- my.loess2(w ~ x + y, data = df, span = 0.5, normalize = FALSE)
#lo.fit13 <- my.loess1(w ~ x + y, data = df, span = 0.5, normalize = FALSE, control = loess.control(surface = "direct"))
#value3 <-predict(lo.fit13, data.frame(x = lo.fit3$kd$vert2$X1, y = lo.fit3$kd$vert2$X2))
