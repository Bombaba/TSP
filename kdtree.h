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

struct kdnearest {
    struct kdnode* node;
    double sqrdist;
};

struct kdtree* build_kdtree(struct point pts[], int n_pts);
void free_kdtree(struct kdtree* tree);
//bool remove_point_from_tree(struct point* p, struct kdtree* tree);
bool remove_point_from_tree(int pindex, struct kdtree* tree);
void print_tree(struct kdtree* tree);
int search_nearest(struct point* p, struct kdtree* tree);
