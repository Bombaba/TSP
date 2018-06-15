#pragma once

#include "point.h"

struct kdnode {
    struct point* point;
    struct kdnode* left;
    struct kdnode* right;
    struct kdnode* parent;
    int dim;
};

struct kdnode* build_kdtree(struct point pts[], int n_pts);
void free_kdtree(struct kdnode* root);

