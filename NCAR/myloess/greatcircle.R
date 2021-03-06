rdist.earth <- function(x1, x2, miles = TRUE, R = NULL) {
	if (is.null(R)) {
		if (miles)
			R <- 3963.34
		else R <- 6378.388
	}
	coslat1 <- cos((x1[, 2] * pi)/180)
	sinlat1 <- sin((x1[, 2] * pi)/180)
	coslon1 <- cos((x1[, 1] * pi)/180)
	sinlon1 <- sin((x1[, 1] * pi)/180)
	if (missing(x2)) {
		pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*%
		t(cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1))
		return(R * acos(ifelse(abs(pp) > 1, 1 * sign(pp), pp)))
	}
	else {
		coslat2 <- cos((x2[, 2] * pi)/180)
		sinlat2 <- sin((x2[, 2] * pi)/180)
		coslon2 <- cos((x2[, 1] * pi)/180)
		sinlon2 <- sin((x2[, 1] * pi)/180)
		pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*%
		t(cbind(coslat2 * coslon2, coslat2 * sinlon2, sinlat2))
		return(R * acos(ifelse(abs(pp) > 1, 1 * sign(pp), pp)))
	}
}
