#pragma once

void build_tour_nn(struct point pts[], int n_pts, int ixstart,
                   int tour[], const struct kdtree* tree);

void build_tour_nn2(struct point pts[], int n_pts, int ixstart, int tour[],
                    const struct kdtree* tree, struct kdheap* heap);
