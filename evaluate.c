#include <math.h>
#include "evaluate.h"

void calc_distribution_param(struct point pts[], int n_pts, struct kdtree *tree, double *mean, double *std, double *skewness, double *kurtosis)
{
	struct point *current, *nearest;
	struct kdtree* tree_copy;
	int i;
	double d1, d2, d3, d4;
    if (tree == NULL) {
        tree_copy = build_kdtree(pts, n_pts);
    } else {
        tree_copy = copy_kdtree(tree);
    }

	d1 = 0.0; d2 = 0.0; d3 = 0.0; d4 = 0.0;
	for(i=0;i<n_pts;i++) {
		double dist;
		current = &pts[i];
		nearest = search_nearest(&pts[i], tree_copy);
		dist = sqrt((current->x-nearest->x)*(current->x-nearest->x)+(current->y-nearest->y)*(current->y-nearest->y));
		d1 += dist;
		d2 += dist * dist;
		d3 += dist * dist * dist;
		d4 += dist * dist * dist * dist;
	}
	d1 = d1 / n_pts; d2 = d2 / n_pts; d3 = d3 / n_pts; d4 = d4 / n_pts;

	*mean = d1;
	*std = sqrt(d2 - d1 * d1);
	*skewness = (d3 - 3 * d1 * d2 + 2 * d1 * d1 * d1) / (*std * *std * *std);
	*kurtosis = (d4 - 4 * d1 * d3 + 6 * d1 * d1 * d2 - 3 * d1 * d1 * d1 * d1) / (*std * *std * *std * *std); 
}
