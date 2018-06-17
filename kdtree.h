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
    int size;
};

struct kdnearest {
    struct kdnode* node;
    double sqrdist;
};

struct kdtree* build_kdtree(struct point pts[], int n_pts);
void free_kdtree(struct kdtree* tree);
bool remove_point_from_tree(int pindex, struct kdtree* tree);
void print_kdtree(const struct kdtree* tree);
struct point* search_nearest(const struct point* p, const struct kdtree* tree);
struct kdtree* copy_kdtree(const struct kdtree* src);
