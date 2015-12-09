
source("aa_kdtree.R")

r = 100
set.seed(r)
p = 4
n = 2^28
x = runif(n*p)  ## unif
x = matrix(x,ncol=p)
nb = 2^18

########### in C
t2 <- system.time({res2 <- cppkdtree(x, nb)})


r = 100
set.seed(r)
p = 4
n = 2^30
x = runif(n*p)  ## unif
x = data.frame(matrix(x,ncol=p))
nb = 2^10

t1 <- system.time({
kdwhole <- my.loess2(
	X1 ~ X2 + X3 + X4, 
	data    = x,
	degree  = 1, 
	span    = 5/nb,
	normalize = FALSE,
	parametric = "X4"
)$kd$vert2
})
