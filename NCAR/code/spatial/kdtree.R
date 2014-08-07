#
kdtree <- function(df, alpha)
{
	v <- data.frame()
	li <- list(0)
	cell <- list()
	depth <- 1
	N <- dim(df)[1]
	D <- dim(df)[2]
	oname <- names(df)[(D-1):D]
	names(df)[(D-1):D] <- c("x","y")
	rx <- range(df$x)
	rx1 <- rx[1] - 0.005*diff(rx)
	rx2 <- rx[2] + 0.005*diff(rx)
	rx <- c(rx1,rx2)
	ry <- range(df$y)
	ry1 <- ry[1] - 0.005*diff(ry)
	ry2 <- ry[2] + 0.005*diff(ry)
	ry <- c(ry1, ry2)
	dif.rx <- diff(rx)
	dif.ry <- diff(ry)
	box.x <- rx
	box.y <- ry
	package <- list(df, box.x, box.y)
	if(depth == 1){
		v <- expand.grid(x=rx, y=ry)
	}
	while(length(li) != 0){
	  length(li) <- length(li) - 1
	  if(dim(package[[1]])[1] > max(alpha*N, 1)){
		rx <- range(package[[1]]$x)
		ry <- range(package[[1]]$y)
		dif.rx <- diff(rx)
		dif.ry <- diff(ry)
		if(dif.rx >= dif.ry){
			m <- median(package[[1]]$x)
			v <- rbind(v, expand.grid(x=m, y=package[[3]]))
			box.x <- c(package[[2]][1], m)
			box.y <- package[[3]]
			li <- append(li, list(list(package[[1]][with(package[[1]], x<=m),], box.x, box.y))) #left
			box.x <- c(m, package[[2]][2])
			box.y <- package[[3]]
			li <- append(li, list(list(package[[1]][with(package[[1]], x>m),], box.x, box.y))) #right
			depth <- depth + 1
		}else{
			m <- median(package[[1]]$y)
			v <- rbind(v, expand.grid(x=package[[2]], y=m))
			box.x <- package[[2]]
			box.y <- c(package[[3]][1], m)
			li <- append(li, list(list(package[[1]][with(package[[1]], y<=m),], box.x, box.y))) #down
			box.x <- package[[2]]
			box.y <- c(m, package[[3]][2])
			li <- append(li, list(list(package[[1]][with(package[[1]], y>m),], box.x, box.y))) #up
			depth <- depth + 1
		}
	  }else{
		names(package[[1]])[(D-1):D] <- oname
		package[[1]]$sub <- sample(1:nrow(package[[1]]))
		cell <- append(cell, list(package[[1]]))
	  }
	  if(length(li) != 0){
		package <- li[[length(li)]]
	  }
	  else{
		package <- list()
	  } 			
	}
	names(v) <- oname
	return(list(vertix=v, cells=cell))
}


