
source("aa_kdtree.R")

r = 100
set.seed(r)
p = 3
n = 2^20
x = runif(n*p)  ## unif
x = matrix(x,ncol=p)
nb = 2^15

########### in C
t2 <- system.time({res2 <- cppkdtree(x, nb)})


