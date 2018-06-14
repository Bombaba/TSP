#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kdtree.h"
#include "point.h"

static void pprint(struct kdnode* node, int depth)
{
    if (node == NULL) return;

    for (int i = 0; i < depth; i++) {
        printf("  ");
    }
    printf("%3d\n", node->point->index);
    pprint(node->left, depth+1);
    pprint(node->right, depth+1);
}


int key_compare_x(const struct point** p1, const struct point** p2)
{
    int i1 = (*p1)->x;
    int i2 = (*p2)->x;

    return (i1 < i2) ? -1 : (i1 > i2) ? 1 : 0;
}

int key_compare_y(const struct point** p1, const struct point** p2)
{
    int i1 = (*p1)->y;
    int i2 = (*p2)->y;

    return (i1 < i2) ? -1 : (i1 > i2) ? 1 : 0;
}

void build_subtree(struct kdnode* node, struct point* p_pts[], int n_pts, int dim)
{
    if (node == NULL || n_pts == 0) return;

    if (n_pts == 1) {
        node->point = p_pts[0];
        return;
    }

    if (dim == 0) {
        qsort(p_pts, n_pts, sizeof(struct point*), key_compare_x);
    } else {
        qsort(p_pts, n_pts, sizeof(struct point*), key_compare_y);
    }

    int mid = n_pts / 2;
    node->point = p_pts[mid];

    build_subtree(node->left, p_pts, mid, 1-dim);
    build_subtree(node->right, p_pts+mid+1, n_pts-mid-1, 1-dim);
}

struct kdnode* build_kdtree(struct point pts[], int n_pts)
{
    if (n_pts == 0) return NULL;

    struct point** p_pts = (struct point**) malloc(sizeof(struct point*) * n_pts);
    for (int i = 0; i < n_pts; i++) {
        p_pts[i] = pts + i;
    }

    struct kdnode* root = (struct kdnode*) malloc(sizeof(struct kdnode) * n_pts + 1);
    for (int i = 0; i < n_pts; i++) {
        root[i].left = (2*i+1 < n_pts) ? &root[2*i+1] : NULL;
        root[i].right = (2*i+2 < n_pts) ? &root[2*i+2] : NULL;
        root[i].parent = &root[(i-1)/2];
    }
    build_subtree(root, p_pts, n_pts, 0);

    pprint(root, 0);

    return root;
}


void free_kdtree(struct kdnode* root)
{
    if (root) free(root);
}
