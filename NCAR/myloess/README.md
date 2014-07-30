### Fortran code ###
- lowesb: kd-tree construction and fitting
v is a long vector with length lv which contains the whole information about kd-tree and some other
things. Memory is assigned by using Calloc() function in C, `loess_workspace()`.
1. vval is vector with length nvmax = max(200, N) in v. It starts at v(iv(13)) in Fortran, which is 
v[iv[12]-1] in C. Length of vval is (d+1) $\times$ nvmax, but useful length is (d+1) $\times$ nv.
2. vert is vector only has max and min of kd-tree vertices for every dimension of predictors. All kd-tree
vertices can be found starting from v(iv(11)) in Fortran, which is v[iv[10]-1] in C, and from 
v(iv(11)+nvmax) in Fortran, which is v[iv[10]-1+nvmax]. Length of vertices is nv which is iv[5].
3. xi is vector of all node points from original predictors. Length of xi is nc which is iv[4]. xi
can be found in v from v(iv(12)) in Fortran, v[iv[11]-1] in C.
- lowese: interpolation based on kd-tree

```
lowesb -> ehg131 --> ehg126(built kd-tree)
                 |-> ehg124
                 |-> ehg139(fit at vertices, vval passed into as s(0:od, nv))
                       |---> ehg127(called for each vertex(nv), s(0:od) is passed into)
```

### R code: "kd" element of loess object ###
- kd$xi is the nodes from original data, if loess(z \~ x+y), xi can be either x-coordinate or 
y-coordinate of the cutting point.
- kd$a specifies which dimension the point in kd$xi comes from.
2 means the second dimension
1 means the first dimension
- kd$vert is the max and min vertex coordinate. For instance, d=2, vert is a vector of
(min(x), min(y), max(x), max(y)) 
- kd$vval is the fitted value for kd-tree vertices. Thera are (d+1) values for each vertex, each of
which is result from a dot product. The first value out of (d+1) is the fitted value at kd-tree
vertices.

