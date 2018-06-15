#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "kdtree.h"
#include "point.h"

enum LR {
    LEFT,
    RIGHT,
    ROOT
};

static void pprint(struct kdnode* node, int depth, enum LR lr)
{
    for (int i = 0; i < depth; i++) {
        printf("| ");
    }

    switch (lr) {
        case LEFT:
            printf("L");
            break;
        case RIGHT:
            printf("R");
            break;
        default:
            break;
    }

    if (node == NULL) {
        printf("[NULL]\n");
        return;
    }

    switch(node->dim) {
        case 0:
            printf("X");
            break;
        case 1:
            printf("Y");
            break;
    }

    printf("[%d:(%d,%d)]\n", node->point->index, node->point->x, node->point->y);
    pprint(node->left, depth+1, LEFT);
    pprint(node->right, depth+1,RIGHT);
}


int key_compare_x(const void* p1, const void* p2)
{
    int i1 = (*(struct point**)p1)->x;
    int i2 = (*(struct point**)p2)->x;

    return (i1 < i2) ? -1 : (i1 > i2) ? 1 : 0;
}

int key_compare_y(const void* p1, const void* p2)
{
    int i1 = (*(struct point**)p1)->y;
    int i2 = (*(struct point**)p2)->y;

    return (i1 < i2) ? -1 : (i1 > i2) ? 1 : 0;
}

struct kdnode* build_subtree(struct kdnode nodes[], struct point* p_pts[], int n_pts, int dim)
{
    if (n_pts == 0) return NULL;

    if (n_pts == 1) {
        struct kdnode* nd = &nodes[p_pts[0]->index];
        nd->dim = dim;
        return nd;
    }

    if (dim == 0) {
        qsort(p_pts, n_pts, sizeof(struct point*), key_compare_x);
    } else {
        qsort(p_pts, n_pts, sizeof(struct point*), key_compare_y);
    }
    int mid = n_pts / 2;
    struct kdnode* nd = &nodes[p_pts[mid]->index];

    nd->dim = dim;
    nd->left = build_subtree(nodes, p_pts, mid, 1-dim);
    if (nd->left) nd->left->parent = nd;

    nd->right = build_subtree(nodes, p_pts+mid+1, n_pts-mid-1, 1-dim);
    if (nd->right) nd->right->parent = nd;

    return nd;
}

struct kdnode* build_kdtree(struct point pts[], int n_pts)
{
    if (n_pts == 0) return NULL;

    struct point** p_pts = (struct point**) malloc(sizeof(struct point*) * n_pts);
    struct kdnode* nodes = (struct kdnode*) malloc(sizeof(struct kdnode) * n_pts + 1);
    
    if (p_pts == NULL || nodes == NULL) return NULL;

    for (int i = 0; i < n_pts; i++) {
        p_pts[i] = pts + i;
        nodes[i].point = pts + i;
        nodes[i].left = NULL;
        nodes[i].right = NULL;
        nodes[i].parent = NULL;
    }
    struct kdnode* root = build_subtree(nodes, p_pts, n_pts, 0);

    pprint(root, 0, ROOT);

    return root;
}

#define kdhead(node) ((node) - (node)->point->index)
#define metric(p, q) (((p).x-(q).x) * ((p).x-(q).x) + ((p).y-(q).y) * ((p).y-(q).y))

void free_kdtree(struct kdnode* root)
{
    if (root == NULL) return;

    free(kdhead(root));
}

struct point* search_nearest(struct point* p, struct kdnode* root)
{
    double min_dist = DBL_MAX;
    return;
} 
