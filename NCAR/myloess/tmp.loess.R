source("my.loess02.R")

set.seed(10)
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
tmp <- lo.fit2$kd$vert2
fit <- setNames(data.frame(matrix(lo.fit2$kd$vval, byrow=TRUE, ncol=3)), c("fitted","dx","dy"))
fit$dxx <- rep(0, ncol(tmp))
fit$dyy <- rep(0, ncol(tmp))
for(i in 1:nrow(tmp)){
	x <- tmp[i,1]
	y <- tmp[i,2]
	yl <- sort(
		tmp[which(tmp[, 1] == x & tmp[, 2] < y), 2], 
		decreasing = TRUE
	)[1]
	yu <- sort(
		tmp[which(tmp[, 1] == x & tmp[, 2] > y), 2]
	)[1]
	xl <- sort(
		tmp[which(tmp[, 2] == y & tmp[, 1] < x), 1], 
		decreasing = TRUE
	)[1]
	xu <- sort(
		tmp[which(tmp[, 2] == y & tmp[, 1] > x), 1]
	)[1]
	N <- c(x, yu)
	S <- c(x, yl)
	E <- c(xu, y)
	W <- c(xl, y)
	F <- fit[i, 1]
	NF <- fit[which(tmp$X1 == N[1] & tmp$X2 == N[2]), 1]
	SF <- fit[which(tmp$X1 == S[1] & tmp$X2 == S[2]), 1]
	EF <- fit[which(tmp$X1 == E[1] & tmp$X2 == E[2]), 1]
	WF <- fit[which(tmp$X1 == W[1] & tmp$X2 == W[2]), 1]
	if(length(NF)!=0){
		dn <- (NF - F)/(yu - y)
	}else{
		dn <- 0
	}
	if(length(SF)!=0){
		ds <- (F - SF)/(y - yl)
	}else{
		ds <- 0
	}
	if(length(EF)!=0){
		de <- (EF - F)/(xu - x)
	}else{
		de <- 0
	}
	if(length(WF)!=0){
		dw <- (F - WF)/(x - xl)
	}else{
		dw <- 0
	}
	if((!is.na(xu)) & (!is.na(xl))){
		dx <- ((xu - x)/(xu - xl))*de + ((x - xl)/(xu - xl))*dw
	}else if(is.na(xu)){
		dx <- dw
	}else if(is.na(xl)){
		dx <- de  
	}
	if((!is.na(yu)) & (!is.na(yl))){
		dy <- ((yu - y)/(yu - yl))*dn + ((y - yl)/(yu - yl))*ds
	}else if(is.na(yu)){
		dy <- ds
	}else if(is.na(yl)){
		dy <- dn
	}
	fit[i, "dxx"] <- dx
	fit[i, "dyy"] <- dy
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
