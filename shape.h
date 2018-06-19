#pragma once

#include "point.h"
//#include "kdtree.h"

//#define REDUCE_THRE 5.0
//#define GRV_THRE 20.0
//#define ALPHA 0.7

void build_reduced_kd_tree(struct point pts[], int n_pts);

void shape_map(struct point pts[], int n_pts, struct point copy_pts[], double grv_thre, double alpha);
