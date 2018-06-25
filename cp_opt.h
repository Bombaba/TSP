#pragma once

#include "point.h"
#include "kdtree.h"

void extract_cp(struct point pts[], int n_pts, 
				struct point c_pts[], struct point p_pts[], int *c_num_ptr, 
				double mean);

void restore_cp(struct point c_pts[], int c_num, int c_tour[], struct point p_pts[],
				int n_pts, int tour[]);
