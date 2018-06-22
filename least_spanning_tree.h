#pragma once

#include <stdlib.h>
#include "point.h"

struct edge {
    int pix;
    int qix;
    double sqrdist;
};

struct lst_edges {
    struct edge* content;
    size_t current;
    size_t size;
};

void build_lst_path(struct point pts[], int n_pts, int* tour);
