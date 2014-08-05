### Fortran code ###
- lowesb: kd-tree construction and fitting
v is a long vector with length lv which contains the whole information about kd-tree and some other
things. Memory is assigned by using Calloc() function in C, `loess_workspace()`.
```
lowesb -> ehg131 -> ehg126(built kd-tree)
                |-> ehg124(not sure what this function is for)
                |-> ehg139(fit at vertices, vval passed into as s(0:od, nv))
                       |-> ehg127(called for each vertex(nv), s(0:od) is passed into)
                       |-> ehg137(try to compare the cutting points xi with vertex)
                       |-> ehg128(interpolation function is called here based on vval2)
```	
  1. vval is vector with length nvmax = max(200, N) in v. It starts at v(iv(13)) in Fortran, which is 
v[iv[12]-1] in C. Length of vval is (d+1)\*nvmax, but useful length is (d+1)\*nv.
  2. vert is vector only has max and min of kd-tree vertices for every dimension of predictors. All kd-tree
vertices can be found starting from v(iv(11)) in Fortran, which is v[iv[10]-1] in C, and from 
v(iv(11)+nvmax) in Fortran, which is v[iv[10]-1+nvmax]. Length of vertices is nv which is iv[5].
  3. xi is vector of all node points from original predictors. Length of xi is nc which is iv[4]. xi
can be found in v from v(iv(12)) in Fortran, v[iv[11]-1] in C.
  4. In ehg127 function, "b" is the design matrix, b(nf, k). k is iv(29), which equal to 
(d+2)\*(d+1)/2, for example two predictors, degree is 2, then there are 6 terms in local 
regression fit. The maximum of k is 15, which means we only can have 4 predictors at most.
  5. In ehg127, for design matrix, a preliminary factorization X = QR into R and Q with Q'Q = I
followed by SVD of R allows the pseudo-inverse to be computed efficiently.
  6. Not sure what is vval2?

- lowese: interpolation based on kd-tree
```
lowese -> ehg133 -> ehg128(interpolation, delta is X for each newobs)
```
  1. For each vertex, the fitted value g(hat) and d derivatives of g(hat) which estimated by taking the
slopes of the locally linear or locally quadratic fit are saved and used to do the interpolation.
  2. Ech cell boundary consists of four segments that meet at vertices. On each segment, function value
g(hat) are interpolated using the unique *cubic polynomial* determined by the function and derivative 
data at the vertices, this cubic polynomial should be an univariate interpolation since there is only one
dimension at edges of cells; normal derivatives are *interpolated linearly* along the segment.
  3. Finally, blending functions interpolate across the cell by using *cubic polynomial* as well. Certain 
cross derivative terms are neeeded, as described by Barnhill(1977), but we have obtained acceptable 
results by setting those cross derivative to 0.
  4. Cubic spline/Cubic interpolation/Cubic Hermite spline:
On the unit interval (0,1), given a starting point p0 at t=0 and an ending point p1 at t=1 with starting 
tangent m0 at t=0 and ending tangent m1 at t=1, the polynomial can be defined by
```
P_t = (2t^3-3t^2+1)P_0 + (t^3-2t^2+t)M_0 + (-2t^3+3t^2)P_1 +(t^3-t^2)M_1 
```
In `ehg128`, the cubic interpolation on boundaries are done as following:
```
c Hermite basis
phi0=(1-h)**2*(1+2*h) --> P_0
phi1=h**2*(3-2*h)     --> P_1
psi0=h*(1-h)**2       --> M_0
psi1=h**2*(h-1)       --> M_1
```
where the `h` is a standardized value of a particular edge of cells. Cubic polynomial of edges of 
cells are done in lop 7.

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
vertices, and the rest of d values are the derivatives of locally linear or locally quadratic fit,
which will be used in interpolation

