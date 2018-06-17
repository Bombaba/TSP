#pragma once

bool two_opt_fast(struct point pts[], int n_pts, double radius,
                  struct kdtree* tree, struct kdheap* heap);
