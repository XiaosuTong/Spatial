In myloess directory, there are two R shared library files: myloess1.so and myloess2.so

myloess1.so is the default one which is loading the original C source functions in loess.

myloess2.so is my shared library file which loading my own C source functions.

`myloess1.so and` and `myloess2.so` are coming from loessc1.c, loessf.f and loessc2.c
loessf.f respectively. `loessc1.c` is same as `loessc.c` which is from the source code
of loess. `loessc2.c` is my personal .c file. In "interpolate/1.approx" section, I included 
a new object named `vert2` in the output list, which is a data.frame including all kd-tree
vertices. At the same time I commented out the `F77_CALL(lowese)` which did the interpolation.

