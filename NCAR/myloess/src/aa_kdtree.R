
###########################
## cpp program to create kd-tree
###########################

## notes:
## 1. "aa_cppkdtree.cpp" is the cpp program
## 2. "aa_auto_compile" is the code to compile and clean up
## 3. "aa_kdtree.R" is a wrapper in R
## 4. run "aa_auto_compile", then dynamic load the .so file

dyn.load("cppkdtree.so")

cppkdtree <- function(data, nb)  ## data matrix, no.of.leaves
{
        D <- ncol(data)  ## number of columns in dm
        ND <- nrow(data)  ## number of rows in dm
        bucketSize <- ND/nb  ## number of rows in each leaf

        ## void getkdtree(double *data, int *D, int *ND, int *bucketSize, int* pidx, int *owner)
        res <- .C("getkdtree"
                , as.numeric(data)  ## data matrix as a vector, by column first
                , as.integer(D)  ## number of columns in dm
                , as.integer(ND)  ## number of rows in dm
                , as.integer(bucketSize)  ## number of rows in each leaf
                , idx = integer(ND)  ## pre-allocate index of rows
                , leaf = integer(ND)  ## pre-allocate leaf
        )  
	## output is a list of 6 elements:
	## 1. dm as a vector
	## 2. number of columns in dm
	## 3. number of rows in dm
	## 4. number of rows in each leaf
	## 5. index of rows
	## 6. index of leaves, can be matched up with index of rows to find leaves

#	return(tapply(res$idx, res$leaf, function(r) data[r,], simplify=FALSE))
#	return(data.frame(idx=res[[5]], leaf=res[[6]]))
	return(res)
	## return the index of rows and index of leaf
}

