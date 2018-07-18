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

void calc_distribution_param_kernel(struct point pts[], int n_pts, 
									double *mean, double *std, double *skewness, double *kurtosis,
									double h, int site)
{
	int i, j, k;
	//int data[site][site];
	int x_max = -1, y_max = -1;
	double x_site, y_site;
    double d1, d2, d3, d4;
	
	for(i=0;i<n_pts;i++) {
		if(pts[i].x > x_max) x_max = pts[i].x;
		if(pts[i].y > y_max) y_max = pts[i].y;
	}
	x_site = x_max / (double)(site);
	y_site = y_max / (double)(site);

    d1 = 0.0; d2 = 0.0; d3 = 0.0; d4 = 0.0;
	for(i=0;i<site;i++) {
		for(j=0;j<site;j++) {
			double px = i*y_site, py = j*x_site;
			double tmp = 0.0;
			for(k=0;k<n_pts;k++) {
				double dist = sqrt((px-pts[k].x) * (px-pts[k].x) + (py-pts[k].y) * (py-pts[k].y));
				tmp += 1000.0/sqrt(2.0*M_PI)/n_pts*pow(M_E, (0-1.0)*(dist/h)*(dist/h)/2.0);
			}
			tmp /= (n_pts*h);
			//printf("%lf\n", tmp);
			d1 += tmp;
			d2 += (tmp * tmp);
			d3 += (tmp * tmp * tmp);
			d4 += (tmp * tmp * tmp * tmp);
		}
	}
    d1 = d1 / (site * site);
	d2 = d2 / (site * site);
	d3 = d3 / (site * site);
	d4 = d4 / (site * site);
	//printf("%lf, %lf\n", d1, d2);

    *mean = d1;
    *std = sqrt(d2 - d1 * d1);
    *skewness = (d3 - 3 * d1 * d2 + 2 * d1 * d1 * d1) / (*std * *std * *std);
    *kurtosis = (d4 - 4 * d1 * d3 + 6 * d1 * d1 * d2 - 3 * d1 * d1 * d1 * d1) / (*std * *std * *std * *std); 
}
