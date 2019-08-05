/* plot_phylo.c (2017-04-25) */

/* Copyright 2004-2017 Emmanuel Paradis */

/* Code from the R-package `ape'. */

#include <R.h>

// [[Rcpp::export]]
void node_depth(int *ntip, int *e1, int *e2,
		int *nedge, double *xx, int *method)
/* method == 1: the node depths are proportional to the number of tips
   method == 2: the node depths are evenly spaced */
{
    int i;

    /* First set the coordinates for all tips */
    for (i = 0; i < *ntip; i++) xx[i] = 1;

    /* Then compute recursively for the nodes; we assume `xx' has */
    /* been initialized with 0's which is true if it has been */
    /* created in R (the tree must be in pruningwise order) */
    if (*method == 1) {
	for (i = 0; i < *nedge; i++)
	    xx[e1[i] - 1] = xx[e1[i] - 1] + xx[e2[i] - 1];
    } else { /* *method == 2 */
	for (i = 0; i < *nedge; i++) {
	    /* if a value > 0 has already been assigned to the ancestor
	       node of this edge, check that the descendant node is not
	       at the same level or more */
	    if (xx[e1[i] - 1])
		if (xx[e1[i] - 1] >= xx[e2[i] - 1] + 1) continue;
	    xx[e1[i] - 1] = xx[e2[i] - 1] + 1;
	}
    }
}
