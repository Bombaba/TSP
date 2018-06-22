#pragma once

bool or_opt(int n, struct point pts[], int n_pts,
            double mul_factor, double add_factor,
            struct kdtree* tree, struct kdheap* heap);

