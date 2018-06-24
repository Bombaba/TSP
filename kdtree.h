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

struct kdnear {
    struct point* point;
    double sqrdist;
};

struct kdheap {
    struct kdnear* content;
    int length;
    int maxsize;
};

struct kdtree* build_kdtree(struct point pts[], int n_pts);
void free_kdtree(struct kdtree* tree);
bool remove_point_from_tree(int pindex, struct kdtree* tree);
//void add_point_toward_tree(
void print_kdtree(const struct kdtree* tree);
struct point* search_nearest(const struct point* p, const struct kdtree* tree);
struct kdtree* copy_kdtree(const struct kdtree* src);
struct kdheap* create_kdheap(struct kdtree* tree);
int search_nearby_points(const struct point* p, const struct kdtree* tree,
                         struct kdheap* heap, double maxdist, int maxsize);
struct point* kdh_pop(struct kdheap* heap);
struct kdnear* kdh_look(int ix, struct kdheap* heap);
void free_kdheap(struct kdheap* heap);
