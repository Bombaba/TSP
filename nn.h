#pragma once

struct point* build_tour_nn(struct point pts[], int n_pts, int start,
                            const struct kdtree* tree);
