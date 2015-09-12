#include <math.h>    // math routines
#include <R.h>       // R header
#include "ANN/ANN.h"     // ANN library header
#include "kd_tree.h"

//------------------------------------------------------------------------------------------------
//				 Near Neighbours Program
//------------------------------------------------------------------------------------------------

extern "C"
{
	void getkdtree(double *data, int *D, int *ND, int *bucketSize, int* pidx, int *owner)
	{
	const int d = *D;		// Number of Dimensions for points
	const int nd = *ND;		// Number of Data points
	const int bs = *bucketSize;
  vector<ANNkd_leaf> vec;
	
	ANNkd_tree *kdtree = NULL;	// Search structure

	ANNpointArray data_pts 	= annAllocPts(nd,d);		// Allocate data points
/*
	int *d_ptr = new int[d];
	
	// set up column offsets for data point matrix (to convert Row/Col major)
	for(int i = 0; i < d; i++)
	{
		d_ptr[i] = i*nd;
	}

	for(int i = 0; i < nd; i++) // now construct the points
	{
		for(int j = 0; j < d; j++)
		{
			data_pts[i][j]=data[ d_ptr[j]++ ];
		}
	}
*/
	
	for(int i = 0; i < nd; i++) // now construct the points
	{
		for(int j = 0; j < d; j++)
		{
			data_pts[i][j]=data[ j*nd + i ];
		}
	}

	kdtree = new ANNkd_tree( data_pts, nd, d, bs, ANN_KD_STD); // Std split rule

  ANNidxArray kdtree_pidx = kdtree->getIdx();
  for (int i = 0; i < nd; i++) {
    pidx[i] = kdtree_pidx[i] + 1;
  }

	vector<ANNkd_leaf>* v = kdtree->getLeaves();
  int o_idx = 0;
  int l_idx = 0;
  for (vector<ANNkd_leaf>::iterator it = v->begin(); it!=v->end(); ++it) {
    for (int i = 0; i < it->getCount(); i++)
      owner[o_idx++] = l_idx + 1;
    l_idx++;
  }

	// Do a little bit of memory management......
	annDeallocPts(data_pts);
//	delete [] d_ptr;
	delete kdtree;
	}
}

