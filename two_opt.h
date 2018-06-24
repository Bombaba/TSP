#pragma once

bool two_opt(struct point pts[], int n_pts,
             double mul_factor, double add_factor,
             const struct kdtree* tree, struct kdheap* heap);
