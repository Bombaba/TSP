#pragma once

#include "point.h"
#include "kdtree.h"

void calc_distribution_param(struct point pts[], int n_pts, struct kdtree *tree, double *mean, double *std, double *skewness, double *kurtosis);
