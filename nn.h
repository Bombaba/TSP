#pragma once

void build_tour_nn(struct point pts[], int n_pts, int ixstart,
                   int tour[], const struct kdtree* tree);
