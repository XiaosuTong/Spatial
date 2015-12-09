//----------------------------------------------------------------------
// File:			kd_tree.h
// Programmer:		Sunil Arya and David Mount
// Description:		Declarations for standard kd-tree routines
// Last modified:	05/03/05 (Version 1.1)
//----------------------------------------------------------------------
// Copyright (c) 1997-2005 University of Maryland and Sunil Arya and
// David Mount.  All Rights Reserved.
// 
// This software and related documentation is part of the Approximate
// Nearest Neighbor Library (ANN).  This software is provided under
// the provisions of the Lesser GNU Public License (LGPL).  See the
// file ../ReadMe.txt for further information.
// 
// The University of Maryland (U.M.) and the authors make no
// representations about the suitability or fitness of this software for
// any purpose.  It is provided "as is" without express or implied
// warranty.
//----------------------------------------------------------------------
// History:
//	Revision 0.1  03/04/98
//		Initial release
//	Revision 1.1  05/03/05
//		Added fixed radius kNN search
//----------------------------------------------------------------------

#ifndef ANN_kd_tree_H
#define ANN_kd_tree_H

#include <ANN/ANNx.h>					// all ANN includes
#include <iostream>
using namespace std;					// make std:: available

//----------------------------------------------------------------------
//	kd-tree:
//		The main search data structure supported by ANN is a kd-tree.
//		The main constructor is given a set of points and a choice of
//		splitting method to use in building the tree.
//
//		Construction:
//		-------------
//		The constructor is given the point array, number of points,
//		dimension, bucket size (default = 1), and the splitting rule
//		(default = ANN_KD_SUGGEST).  The point array is not copied, and
//		is assumed to be kept constant throughout the lifetime of the
//		search structure.  There is also a "load" constructor that
//		builds a tree from a file description that was created by the
//		Dump operation.
//
//		Search:
//		-------
//		There are two search methods:
//
//			Standard search (annkSearch()):
//				Searches nodes in tree-traversal order, always visiting
//				the closer child first.
//			Priority search (annkPriSearch()):
//				Searches nodes in order of increasing distance of the
//				associated cell from the query point.  For many
//				distributions the standard search seems to work just
//				fine, but priority search is safer for worst-case
//				performance.
//
//		Printing:
//		---------
//		There are two methods provided for printing the tree.  Print()
//		is used to produce a "human-readable" display of the tree, with
//		indenation, which is handy for debugging.  Dump() produces a
//		format that is suitable reading by another program.  There is a
//		"load" constructor, which constructs a tree which is assumed to
//		have been saved by the Dump() procedure.
//		
//		Performance and Structure Statistics:
//		-------------------------------------
//		The procedure getStats() collects statistics information on the
//		tree (its size, height, etc.)  See ANNperf.h for information on
//		the stats structure it returns.
//
//		Internal information:
//		---------------------
//		The data structure consists of three major chunks of storage.
//		The first (implicit) storage are the points themselves (pts),
//		which have been provided by the users as an argument to the
//		constructor, or are allocated dynamically if the tree is built
//		using the load constructor).  These should not be changed during
//		the lifetime of the search structure.  It is the user's
//		responsibility to delete these after the tree is destroyed.
//
//		The second is the tree itself (which is dynamically allocated in
//		the constructor) and is given as a pointer to its root node
//		(root).  These nodes are automatically deallocated when the tree
//		is deleted.  See the file src/kd_tree.h for further information
//		on the structure of the tree nodes.
//
//		Each leaf of the tree does not contain a pointer directly to a
//		point, but rather contains a pointer to a "bucket", which is an
//		array consisting of point indices.  The third major chunk of
//		storage is an array (pidx), which is a large array in which all
//		these bucket subarrays reside.  (The reason for storing them
//		separately is the buckets are typically small, but of varying
//		sizes.  This was done to avoid fragmentation.)  This array is
//		also deallocated when the tree is deleted.
//
//		In addition to this, the tree consists of a number of other
//		pieces of information which are used in searching and for
//		subsequent tree operations.  These consist of the following:
//
//		dim						Dimension of space
//		n_pts					Number of points currently in the tree
//		n_max					Maximum number of points that are allowed
//								in the tree
//		bkt_size				Maximum bucket size (no. of points per leaf)
//		bnd_box_lo				Bounding box low point
//		bnd_box_hi				Bounding box high point
//		splitRule				Splitting method used
//
//----------------------------------------------------------------------

//----------------------------------------------------------------------
// Some types and objects used by kd-tree functions
// See src/kd_tree.h and src/kd_tree.cpp for definitions
//----------------------------------------------------------------------
class ANNkdStats;				// stats on kd-tree
class ANNkd_node;
typedef ANNkd_node*	ANNkd_ptr;	// pointer to a kd-tree node
class ANNkd_leaf;

//----------------------------------------------------------------------
//	kd-splitting function:
//		kd_splitter is a pointer to a splitting routine for preprocessing.
//		Different splitting procedures result in different strategies
//		for building the tree.
//----------------------------------------------------------------------

typedef void (*ANNkd_splitter)(			// splitting routine for kd-trees
	ANNpointArray		pa,				// point array (unaltered)
	ANNidxArray			pidx,			// point indices (permuted on return)
	const ANNorthRect	&bnds,			// bounding rectangle for cell
	int					n,				// number of points
	int					dim,			// dimension of space
	int					&cut_dim,		// cutting dimension (returned)
	ANNcoord			&cut_val,		// cutting value (returned)
	int					&n_lo);			// num of points on low side (returned)

class DLL_API ANNkd_tree: public ANNpointSet {
protected:
	int				dim;				// dimension of space
	int				n_pts;				// number of points in tree
	int				bkt_size;			// bucket size
	ANNpointArray	pts;				// the points
	ANNidxArray		pidx;				// point indices (to pts array)
	ANNkd_ptr		root;				// root of kd-tree
	ANNpoint		bnd_box_lo;			// bounding box low point
	ANNpoint		bnd_box_hi;			// bounding box high point
  vector<ANNkd_leaf> leaves; // vector to store leaves

	void SkeletonTree(					// construct skeleton tree
		int				n,				// number of points
		int				dd,				// dimension
		int				bs,				// bucket size
		ANNpointArray pa = NULL,		// point array (optional)
		ANNidxArray pi = NULL);			// point indices (optional)

  ANNkd_ptr rkd_tree(				// recursive construction of kd-tree
    ANNpointArray		pa,				// point array
    ANNidxArray			pidx,			// point indices to store in subtree
    int					n,				// number of points
    int					dim,			// dimension of space
    int					bsp,			// bucket space
    ANNorthRect			&bnd_box,		// bounding box for current node
    ANNkd_splitter		splitter);		// splitting routine
public:
	ANNkd_tree(							// build skeleton tree
		int				n = 0,			// number of points
		int				dd = 0,			// dimension
		int				bs = 1);		// bucket size

	ANNkd_tree(							// build from point array
		ANNpointArray	pa,				// point array
		int				n,				// number of points
		int				dd,				// dimension
		int				bs = 1,			// bucket size
		ANNsplitRule	split = ANN_KD_SUGGEST);	// splitting method

	ANNkd_tree(							// build from dump file
		std::istream&	in);			// input stream for dump file

	~ANNkd_tree();						// tree destructor

	void annkSearch(					// approx k near neighbor search
		ANNpoint		q,				// query point
		int				k,				// number of near neighbors to return
		ANNidxArray		nn_idx,			// nearest neighbor array (modified)
		ANNdistArray	dd,				// dist to near neighbors (modified)
		double			eps=0.0);		// error bound

	void annkPriSearch( 				// priority k near neighbor search
		ANNpoint		q,				// query point
		int				k,				// number of near neighbors to return
		ANNidxArray		nn_idx,			// nearest neighbor array (modified)
		ANNdistArray	dd,				// dist to near neighbors (modified)
		double			eps=0.0);		// error bound

	int annkFRSearch(					// approx fixed-radius kNN search
		ANNpoint		q,				// the query point
		ANNdist			sqRad,			// squared radius of query ball
		int				k,				// number of neighbors to return
		ANNidxArray		nn_idx = NULL,	// nearest neighbor array (modified)
		ANNdistArray	dd = NULL,		// dist to near neighbors (modified)
		double			eps=0.0);		// error bound

	int theDim()						// return dimension of space
		{ return dim; }

	int nPoints()						// return number of points
		{ return n_pts; }

	ANNpointArray thePoints()			// return pointer to points
		{ return pts; }

	virtual void Print(					// print the tree (for debugging)
		ANNbool			with_pts,		// print points as well?
		std::ostream&	out);			// output stream

	virtual void Dump(					// dump entire tree
		ANNbool			with_pts,		// print points as well?
		std::ostream&	out);			// output stream
								
	virtual void getStats(				// compute tree statistics
		ANNkdStats&		st);			// the statistics (modified)

	vector<ANNkd_leaf> *getLeaves()
		{ return &leaves; }

	ANNidxArray getIdx()
    { return pidx; }
};

//----------------------------------------------------------------------
//	Box decomposition tree (bd-tree)
//		The bd-tree is inherited from a kd-tree.  The main difference
//		in the bd-tree and the kd-tree is a new type of internal node
//		called a shrinking node (in the kd-tree there is only one type
//		of internal node, a splitting node).  The shrinking node
//		makes it possible to generate balanced trees in which the
//		cells have bounded aspect ratio, by allowing the decomposition
//		to zoom in on regions of dense point concentration.  Although
//		this is a nice idea in theory, few point distributions are so
//		densely clustered that this is really needed.
//----------------------------------------------------------------------

class DLL_API ANNbd_tree: public ANNkd_tree {
public:
	ANNbd_tree(							// build skeleton tree
		int				n,				// number of points
		int				dd,				// dimension
		int				bs = 1)			// bucket size
		: ANNkd_tree(n, dd, bs) {}		// build base kd-tree

	ANNbd_tree(							// build from point array
		ANNpointArray	pa,				// point array
		int				n,				// number of points
		int				dd,				// dimension
		int				bs = 1,			// bucket size
		ANNsplitRule	split  = ANN_KD_SUGGEST,	// splitting rule
		ANNshrinkRule	shrink = ANN_BD_SUGGEST);	// shrinking rule

	ANNbd_tree(							// build from dump file
		std::istream&	in);			// input stream for dump file
};

//----------------------------------------------------------------------
//	Generic kd-tree node
//
//		Nodes in kd-trees are of two types, splitting nodes which contain
//		splitting information (a splitting hyperplane orthogonal to one
//		of the coordinate axes) and leaf nodes which contain point
//		information (an array of points stored in a bucket).  This is
//		handled by making a generic class kd_node, which is essentially an
//		empty shell, and then deriving the leaf and splitting nodes from
//		this.
//----------------------------------------------------------------------

class ANNkd_node{						// generic kd-tree node (empty shell)
public:
	virtual ~ANNkd_node() {}					// virtual distroyer

	virtual void ann_search(ANNdist) = 0;		// tree search
	virtual void ann_pri_search(ANNdist) = 0;	// priority search
	virtual void ann_FR_search(ANNdist) = 0;	// fixed-radius search

	virtual void getStats(						// get tree statistics
				int dim,						// dimension of space
				ANNkdStats &st,					// statistics
				ANNorthRect &bnd_box) = 0;		// bounding box
												// print node
	virtual void print(int level, ostream &out) = 0;
	virtual void dump(ostream &out) = 0;		// dump node

	friend class ANNkd_tree;					// allow kd-tree to access us
};

//----------------------------------------------------------------------
//	Leaf kd-tree node
//		Leaf nodes of the kd-tree store the set of points associated
//		with this bucket, stored as an array of point indices.  These
//		are indices in the array points, which resides with the
//		root of the kd-tree.  We also store the number of points
//		that reside in this bucket.
//----------------------------------------------------------------------

class ANNkd_leaf: public ANNkd_node		// leaf node for kd-tree
{
	int					n_pts;			// no. points in bucket
	ANNidxArray			bkt;			// bucket of points
public:
	ANNkd_leaf(							// constructor
		int				n,				// number of points
		ANNidxArray		b)				// bucket
  {
    n_pts		= n;			// number of points in bucket
    bkt			= b;			// the bucket
  }
  
	virtual void getStats(						// get tree statistics
				int dim,						// dimension of space
				ANNkdStats &st,					// statistics
				ANNorthRect &bnd_box);			// bounding box
	virtual void print(int level, ostream &out);// print node
	virtual void dump(ostream &out);			// dump node

	virtual void ann_search(ANNdist);			// standard search
	virtual void ann_pri_search(ANNdist);		// priority search
	virtual void ann_FR_search(ANNdist);		// fixed-radius search
  friend ostream& operator<<(ostream& stream, ANNkd_leaf leaf)
  {
    for (int i = 0; i < leaf.n_pts; i++)
      stream << leaf.bkt[i] << ",";
    //go back and erase the last comma from stream
    stream << '\b';
    stream << ' ';
    return stream;
  }
  int getCount()
  {
    return n_pts;
  }
};

//----------------------------------------------------------------------
//		KD_TRIVIAL is a special pointer to an empty leaf node. Since
//		some splitting rules generate many (more than 50%) trivial
//		leaves, we use this one shared node to save space.
//
//		The pointer is initialized to NULL, but whenever a kd-tree is
//		created, we allocate this node, if it has not already been
//		allocated. This node is *never* deallocated, so it produces
//		a small memory leak.
//----------------------------------------------------------------------

extern ANNkd_leaf *KD_TRIVIAL;					// trivial (empty) leaf node

//----------------------------------------------------------------------
//	kd-tree splitting node.
//		Splitting nodes contain a cutting dimension and a cutting value.
//		These indicate the axis-parellel plane which subdivide the
//		box for this node. The extent of the bounding box along the
//		cutting dimension is maintained (this is used to speed up point
//		to box distance calculations) [we do not store the entire bounding
//		box since this may be wasteful of space in high dimensions].
//		We also store pointers to the 2 children.
//----------------------------------------------------------------------

class ANNkd_split : public ANNkd_node	// splitting node of a kd-tree
{
	int					cut_dim;		// dim orthogonal to cutting plane
	ANNcoord			cut_val;		// location of cutting plane
	ANNcoord			cd_bnds[2];		// lower and upper bounds of
										// rectangle along cut_dim
	ANNkd_ptr			child[2];		// left and right children
public:
	ANNkd_split(						// constructor
		int cd,							// cutting dimension
		ANNcoord cv,					// cutting value
		ANNcoord lv, ANNcoord hv,				// low and high values
		ANNkd_ptr lc=NULL, ANNkd_ptr hc=NULL)	// children
		{
			cut_dim		= cd;					// cutting dimension
			cut_val		= cv;					// cutting value
			cd_bnds[ANN_LO] = lv;				// lower bound for rectangle
			cd_bnds[ANN_HI] = hv;				// upper bound for rectangle
			child[ANN_LO]	= lc;				// left child
			child[ANN_HI]	= hc;				// right child
		}

	~ANNkd_split()						// destructor
		{
			if (child[ANN_LO]!= NULL && child[ANN_LO]!= KD_TRIVIAL)
				delete child[ANN_LO];
			if (child[ANN_HI]!= NULL && child[ANN_HI]!= KD_TRIVIAL)
				delete child[ANN_HI];
		}

	virtual void getStats(						// get tree statistics
				int dim,						// dimension of space
				ANNkdStats &st,					// statistics
				ANNorthRect &bnd_box);			// bounding box
	virtual void print(int level, ostream &out);// print node
	virtual void dump(ostream &out);			// dump node

	virtual void ann_search(ANNdist);			// standard search
	virtual void ann_pri_search(ANNdist);		// priority search
	virtual void ann_FR_search(ANNdist);		// fixed-radius search
};

//----------------------------------------------------------------------
//		External entry points
//----------------------------------------------------------------------

ANNkd_ptr rkd_tree(				// recursive construction of kd-tree
	ANNpointArray		pa,				// point array (unaltered)
	ANNidxArray			pidx,			// point indices to store in subtree
	int					n,				// number of points
	int					dim,			// dimension of space
	int					bsp,			// bucket space
	ANNorthRect			&bnd_box,		// bounding box for current node
	ANNkd_splitter		splitter);		// splitting routine

#endif
