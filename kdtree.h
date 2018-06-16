#pragma once

#include "point.h"

struct kdnode {
    struct point* point;
    struct kdnode* left;
    struct kdnode* right;
    struct kdnode* parent;
    int dim;
    bool valid;
};

struct kdtree {
    struct kdnode* head;
    struct kdnode* root;
    int max_index;
};

struct kdtree* build_kdtree(struct point pts[], int n_pts);
void free_kdtree(struct kdtree* tree);
bool remove_point_from_tree(int point_index, struct kdtree* tree);

