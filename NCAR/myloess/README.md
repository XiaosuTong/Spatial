# Fortran code #
lowesb: kd-tree construction and fitting

lowese: interpolation based on kd-tree

# R code: "kd" element of loess object #
- kd$xi is the nodes from original data, if loess(z \~ x+y), xi can be either x-coordinate or 
y-coordinate of the cutting point.
- kd$a specifies which dimension the point in kd$xi comes from.
2 means the second dimension
1 means the first dimension
- kd$vert is the max and min vertex coordinate. For instance, d=2, vert is a vector of
(min(x), min(y), max(x), max(y)) 


