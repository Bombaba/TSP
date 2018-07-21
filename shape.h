#pragma once

#include "point.h"

/*
How to use "reduce map":
	int reduce_num;
	struct kdtree* tree2;
	int copy_tour[MAX_N];
	struct point reduced_pts[MAX_N];
	struct point reduce_memo[MAX_N];

	reduce_num = reduce_map(pts, n_pts, reduced_pts, reduce_memo, 5.0);
	tree2 = build_kdtree(reduced_pts, n_pts-reduce_num);

	build_tour_nn(reduced_pts, n_pts-reduce_num, start, tour, tree2);
	restore_reduced_tour(reduced_pts, reduce_memo, reduce_num, tour, copy_tour, n_pts);

	double length = tour_length(pts, n_pts-reduce_num, copy_tour);
	
	memcpy(best_tour, copy_tour, sizeof(int) * n_pts);
*/
int reduce_map(struct point pts[], int n_pts, struct point copy_pts[], 
				struct point reduced_pts[], const int cp[], int cn, int copy_cp[], double thre);

void restore_reduced_tour(struct point pts[], struct point reduced_pts[], int n_rpts, int tour[], int copy_tour[], int n_pts);


/*
How to use "shape map":
	struct point shaped_pts[MAX_n];
	shape_map(pts, n_pts, shaped_pts, 100.0, 0.6);

	build_tour_nn(shaped_pts, n_pts, start, tour, tree);

	double length = tour_length(pts, n_pts, tour);
*/
void shape_map(struct point pts[], int n_pts, struct point copy_pts[], 
			   const int cp[], int cn, double grv_thre, double alpha);
//void shape_map(struct point pts[], int n_pts, struct point copy_pts[], double grv_thre, double alpha);
