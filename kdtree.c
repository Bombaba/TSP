#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <float.h>
#include "kdtree.h"
#include "point.h"

#define metric(p, q) (((p).x-(q).x) * ((p).x-(q).x) + ((p).y-(q).y) * ((p).y-(q).y))
#define point_to_kdnode(p, tree) ((tree)->head + (p)->index)

struct kdnode* find_min(struct kdnode* node, int dim);

enum LR {
    LEFT,
    RIGHT,
    ROOT
};

static void print_tree(struct kdnode* node, int depth, enum LR lr)
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
        //printf("[NULL]\n");
        putchar('\n');
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
    print_tree(node->left, depth+1, LEFT);
    print_tree(node->right, depth+1,RIGHT);
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
        nd->valid = true;
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
    nd->valid = true;
    nd->left = build_subtree(nodes, p_pts, mid, 1-dim);
    if (nd->left) nd->left->parent = nd;

    nd->right = build_subtree(nodes, p_pts+mid+1, n_pts-mid-1, 1-dim);
    if (nd->right) nd->right->parent = nd;

    return nd;
}

struct kdtree* build_kdtree(struct point pts[], int n_pts)
{
    if (n_pts == 0) return NULL;

    struct point** p_pts = (struct point**) malloc(sizeof(struct point*) * n_pts);
    struct kdnode* head = (struct kdnode*) malloc(sizeof(struct kdnode) * n_pts + 1);
    
    if (p_pts == NULL || head == NULL) return NULL;

    for (int i = 0; i < n_pts; i++) {
        p_pts[i] = pts + i;
        head[i].point = pts + i;
        head[i].left = NULL;
        head[i].right = NULL;
        head[i].parent = NULL;
        head[i].dim = -1;
        head[i].valid = false;
    }
    struct kdnode* root = build_subtree(head, p_pts, n_pts, 0);
    free(p_pts);

    struct kdtree* tree = (struct kdtree*) malloc(sizeof(struct kdtree));
    tree->head = head;
    tree->root = root;
    tree->max_index = n_pts - 1;

    print_tree(tree->root, 0, ROOT);

    putchar('\n');
    remove_point_from_tree(6, tree);
    remove_point_from_tree(6, tree);
    remove_point_from_tree(100, tree);
    print_tree(tree->root, 0, ROOT);
    remove_point_from_tree(13, tree);
    remove_point_from_tree(4, tree);
    print_tree(tree->root, 0, ROOT);
    remove_point_from_tree(17, tree);
    remove_point_from_tree(3, tree);
    remove_point_from_tree(15, tree);
    remove_point_from_tree(1, tree);
    remove_point_from_tree(9, tree);
    print_tree(tree->root, 0, ROOT);

    return tree;
}


void free_kdtree(struct kdtree* tree)
{
    if (tree->head == NULL) return;

    free(tree->head);
}


struct kdnode* find_min(struct kdnode* node, int dim)
{
    if (node == NULL) return NULL;

    if (node->dim == dim) {
        if (node->left != NULL) {
            return find_min(node->left, dim);
        }
        return node;
    }

    struct kdnode* left = find_min(node->left, dim);
    struct kdnode* right = find_min(node->right, dim);

    struct kdnode* min = node;
    if (left != NULL && left->point->pos[dim] < min->point->pos[dim]) {
        min = left;
    }
    if (right != NULL && right->point->pos[dim] < min->point->pos[dim]) {
        min = right;
    }

    return min;
}

void remove_node(struct kdnode* node, struct kdtree* tree)
{
    //printf("*** Removed point #%d ***\n", node->point->index);

    node->valid = false;

    if (node->left == NULL && node->right == NULL) {
        if (node->parent == NULL) {
            tree->root = NULL;
        } else if (node->parent->left == node) {
            node->parent->left = NULL;
        } else {
            node->parent->right = NULL;
        }
        
        return;
    }

    if (node->right == NULL) {
        node->right = node->left;
        node->left = NULL;
    }

    struct kdnode* next = find_min(node->right, node->dim);
    remove_node(next, tree);
    next->valid = true;
    next->parent = node->parent;
    next->left = node->left;
    next->right = node->right;
    next->dim = node->dim;
    if (node->parent == NULL) {
        tree->root = next;
    } else if (node->parent->left == node) {
        node->parent->left = next;
    } else {
        node->parent->right = next;
    }
}

bool remove_point_from_tree(int point_index, struct kdtree* tree)
{
    if (point_index > tree->max_index) {
        fprintf(
            stderr,
            "Error: Trying to remove point#%d from a kdtree,\n"
            "       but the maximum index of the tree is #%d\n",
            point_index, tree->max_index
        );
        return false;
    }

    struct kdnode* node = tree->head + point_index;

    if (node->valid == false) {
        fprintf(
            stderr,
            "Warning: Trying to remove point#%d from a kdtree,\n"
            "         but point#%d was already removed from the tree.\n",
            point_index, point_index
        );
        return false;
    }

    remove_node(node, tree);
    return true;
}

struct point* search_nearest(struct point* p, struct kdtree* tree)
{
    double min_dist = DBL_MAX;
    return p;
} 
