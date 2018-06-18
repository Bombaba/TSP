#pragma once

bool two_opt(struct point pts[], int n_pts, double radius,
             const struct kdtree* tree, struct kdheap* heap);
