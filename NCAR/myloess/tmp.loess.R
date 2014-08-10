source("my.loess02.R")

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

library(lattice)
xyplot(X2 ~ X1,
	data = lo.fit2$kd$vert2,
	panel = function(x,y,...){
		panel.xyplot(x,y,...)
		for(i in seq(1,(length(x)),2)) {
			if( i ==1){
				panel.segments(x[i],y[i],x[i+2],y[i+2])
			}else if(i ==3){
				panel.segments(x[i-1],y[i-1],x[i+1],y[i+1])
			}
			panel.segments(x[i],y[i],x[i+1],y[i+1])
		}
	}
)

tmp <- lo.fit2$kd$vert2[with(lo.fit2$kd$vert2, order(X2,X1)),]

for(i in 1:ncol(tmp)){
	x <- tmp[i,1]
	y <- tmp[i,2]
	yu <- sort(tmp[which(tmp[, 1] == x & tmp[, 2] < y), 2], decreasing = TRUE)[1]
	yl <- sort(tmp[which(tmp[, 1] == x & tmp[, 2] > y), 2])[1]
	xu <- sort(tmp[which(tmp[, 2])])

}

#kd <- lo.fit2$kd$vert2
#value <-predict(lo.fit1, data.frame(x = kd$X1, y = kd$X2))
#tmp <- data.frame(matrix(lo.fit2$kd$vval, byrow=TRUE, ncol=3))
#names(tmp) <- c("b0", "b1", "b2")
#fit <- cbind(tmp, kd)
#fit$fitted <- with(fit, b1*X1 + b2*X2)
#lo.fit3 <- my.loess2(w ~ x + y, data = df, span = 0.5, normalize = FALSE)
#lo.fit13 <- my.loess1(w ~ x + y, data = df, span = 0.5, normalize = FALSE, control = loess.control(surface = "direct"))
#value3 <-predict(lo.fit13, data.frame(x = lo.fit3$kd$vert2$X1, y = lo.fit3$kd$vert2$X2))
