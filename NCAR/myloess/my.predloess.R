##############################################
##predict loess function
##############################################
my.predict.loess <- function (object, newdata = NULL, se = FALSE, na.action = na.pass, 
    ...) 
{
    if (!inherits(object, "loess")) 
        stop("first argument must be a \"loess\" object")
    if (is.null(newdata) && !se) 
        return(fitted(object))
    newx <- if (is.null(newdata)) 
        object$x
    else if (is.data.frame(newdata)) 
        as.matrix(model.frame(delete.response(terms(object)), 
            newdata, na.action = na.action))
    else as.matrix(newdata)
    res <- with(object, my.predLoess(y, x, newx, s, weights, pars$robust, 
        pars$span, pars$degree, pars$normalize, pars$parametric, 
        pars$drop.square, pars$surface, pars$cell, pars$family, 
        kd, divisor, se = se))
    if (!is.null(out.attrs <- attr(newdata, "out.attrs"))) {
        if (se) {
            res$fit <- array(res$fit, out.attrs$dim, out.attrs$dimnames)
            res$se.fit <- array(res$se.fit, out.attrs$dim, out.attrs$dimnames)
        }
        else res <- array(res, out.attrs$dim, out.attrs$dimnames)
    }
    if (se) 
        res$df <- object$one.delta^2/object$two.delta
    res
}

my.predLoess <- function (y, x, newx, s, weights, robust, span, degree, normalize, 
    parametric, drop.square, surface, cell, family, kd, divisor, 
    se = FALSE) 
{
    D <- NCOL(x)
    N <- NROW(x)
    M <- NROW(newx)
    x <- as.matrix(x)
    newx <- as.matrix(newx)
    newx <- newx/rep(divisor, rep(M, D))
    x <- x/rep(divisor, rep(N, D))
    sum.drop.sqr <- sum(drop.square)
    nonparametric <- sum(!parametric)
    order.parametric <- order(parametric)
    x <- x[, order.parametric, drop = FALSE]
    x.evaluate <- newx[, order.parametric, drop = FALSE]
    order.drop.sqr <- (2L - drop.square)[order.parametric]
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    if (surface == "direct") {
        nas <- rowSums(is.na(newx)) > 0L
        fit <- rep(NA_real_, length(nas))
        x.evaluate <- x.evaluate[!nas, , drop = FALSE]
        M <- nrow(x.evaluate)
        if (se) {
            se.fit <- fit
            z <- .C("loess_dfitse", y, x, as.double(x.evaluate), 
                as.double(weights * robust), as.double(robust), 
                as.integer(family == "gaussian"), as.double(span), 
                as.integer(degree), as.integer(nonparametric), 
                as.integer(order.drop.sqr), as.integer(sum.drop.sqr), 
                as.integer(D), as.integer(N), as.integer(M), 
                fit = double(M), L = double(N * M))[c("fit", 
                "L")]
            fit[!nas] <- z$fit
            ses <- (matrix(z$L^2, M, N)/rep(weights, rep(M, N))) %*% 
                rep(1, N)
            se.fit[!nas] <- drop(s * sqrt(ses))
        }
        else {
            fit[!nas] <- .C("loess_dfit", y, x, as.double(x.evaluate), 
                as.double(weights * robust), as.double(span), 
                as.integer(degree), as.integer(nonparametric), 
                as.integer(order.drop.sqr), as.integer(sum.drop.sqr), 
                as.integer(D), as.integer(N), as.integer(M), 
                fit = double(M))$fit
        }
    }
    else {
        inside <- matrix(FALSE, M, ncol = D)
        ranges <- apply(x, 2L, range)
        inside <- (x.evaluate <= rep(ranges[2L, ], rep(M, D))) & 
            (x.evaluate >= rep(ranges[1L, ], rep(M, D)))
        inside <- inside %*% rep(1, D) == D
        inside[is.na(inside)] <- FALSE
        #M1 is the number of new points want to be fitted
        M1 <- sum(inside)
        fit <- rep(NA_real_, M)
        if (any(inside)) 
            fit[inside] <- .C("loess_ifit", as.integer(kd$parameter), 
                as.integer(kd$a), as.double(kd$xi), as.double(kd$vert), 
                as.double(kd$vval), as.integer(M1), as.double(x.evaluate[inside, 
                  ]), fit = double(M1))$fit
        if (se) {
            se.fit <- rep(NA_real_, M)
            if (any(inside)) {
                L <- .C("loess_ise", y, x, as.double(x.evaluate[inside, 
                  ]), as.double(weights), as.double(span), as.integer(degree), 
                  as.integer(nonparametric), as.integer(order.drop.sqr), 
                  as.integer(sum.drop.sqr), as.double(span * 
                    cell), as.integer(D), as.integer(N), as.integer(M1), 
                  double(M1), L = double(N * M1))$L
                tmp <- (matrix(L^2, M1, N)/rep(weights, rep(M1, 
                  N))) %*% rep(1, N)
                se.fit[inside] <- drop(s * sqrt(tmp))
            }
        }
    }
    rn <- rownames(newx)
    if (se) {
        if (!is.null(rn)) 
            names(fit) <- names(se.fit) <- rn
        list(fit = fit, se.fit = drop(se.fit), residual.scale = s)
    }
    else {
        if (!is.null(rn)) 
            names(fit) <- rn
        fit
    }
}


