#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <float.h>
#include "kdtree.h"
#include "point.h"


struct kdnode* find_min(struct kdnode* node, int dim);

void ptree(const struct kdnode* node, int depth)
{
    if (node == NULL) return;

    for (int i = 0; i < depth; i++) {
        printf("| ");
    }

    //if (node == NULL) {
    //    printf("[NULL]\n");
    //    return;
    //}

    switch(node->dim) {
        case 0:
            printf("X");
            break;
        case 1:
            printf("Y");
            break;
    }

    printf("[#%d:(%d,%d) -> #%d]\n",
           node->point->index, node->point->x, node->point->y,
           depth ? node->parent->point->index : -1);
    ptree(node->left, depth+1);
    ptree(node->right, depth+1);
}

void print_kdtree(const struct kdtree* tree)
{
    ptree(tree->root, 0);
}

struct kdtree* copy_kdtree(const struct kdtree* src)
{
    struct kdnode* shead = src->head;
    struct kdnode* dhead = (struct kdnode*) malloc(sizeof(struct kdnode) * src->size);

    for (int i = 0; i < src->size; i++) {
        dhead[i].point = shead[i].point;
        dhead[i].left = shead[i].left ? (dhead + (shead[i].left - shead)) : NULL;
        dhead[i].right = shead[i].right ? (dhead + (shead[i].right - shead)) : NULL;
        dhead[i].parent = shead[i].parent ? (dhead + (shead[i].parent - shead)) : NULL;
        dhead[i].dim = shead[i].dim;
        dhead[i].valid = shead[i].valid;
    }

    struct kdtree* dst = (struct kdtree*) malloc(sizeof(struct kdtree));
    dst->head = dhead;
    dst->root = dst->head + src->root->point->index;
    dst->size = src->size;
    return dst;
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
    tree->size = n_pts;

    return tree;
}


void free_kdtree(struct kdtree* tree)
{
    if (tree->head == NULL) return;

    free(tree->head);
    free(tree);
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
    if (node->left) node->left->parent = next;
    if (node->right) node->right->parent = next;

    if (node->parent == NULL) {
        tree->root = next;
    } else if (node->parent->left == node) {
        node->parent->left = next;
    } else {
        node->parent->right = next;
    }
}

//bool remove_point_from_tree(struct point* p, struct kdtree* tree)
bool remove_point_from_tree(int pindex, struct kdtree* tree)
{
    //int pindex = p->index;
    if (pindex >= tree->size) {
        fprintf(
            stderr,
            "Error: Trying to remove point#%d from the kdtree,\n"
            "       but the maximum index of the tree is #%d\n",
            pindex, tree->size-1
        );
        exit(EXIT_FAILURE);
    }

    struct kdnode* node = tree->head + pindex;

    if (node->valid == false) {
        fprintf(
            stderr,
            "Warning: Trying to remove point#%d from the kdtree,\n"
            "         but point#%d was already removed.\n",
            pindex, pindex
        );
        return false;
    }

    remove_node(node, tree);
    return true;
}

void nearest_point(const struct point* p, struct kdnode* node, struct kdnear* kdn)
{
    double sqrdist = metric(p, node->point);

    if (node->left == NULL && node->right == NULL) {
        if (p != node->point && sqrdist < kdn->sqrdist) {
            kdn->point = node->point;
            kdn->sqrdist = sqrdist;
        }
        return;
    }

    int dim = node->dim;

    struct kdnode* next;
    if (node->right == NULL) {
        next = node->left;
    } else if (node->left == NULL) {
        next = node->right;
    } else {
        if (p->pos[dim] < node->point->pos[dim]) {
            next = node->left;
        } else {
            next = node->right;
        }
    }
    nearest_point(p, next, kdn);

    if (p != node->point && sqrdist < kdn->sqrdist) {
        kdn->point = node->point;
        kdn->sqrdist = sqrdist;
    }

    double lin_sqrdist = p->pos[dim] - node->point->pos[dim];
    lin_sqrdist *= lin_sqrdist;

    if (lin_sqrdist < kdn->sqrdist) {
        if(next == node->left) {
            next = node->right;
        } else {
            next = node->left;
        }
        if (next != NULL) {
            nearest_point(p, next, kdn);
        }
    }
}

struct point* search_nearest(const struct point* p, const struct kdtree* tree)
{
    struct kdnear kdn;
    kdn.point = NULL;
    kdn.sqrdist = DBL_MAX;

    if (tree->root == NULL) return NULL;

    nearest_point(p, tree->root, &kdn);

    return kdn.point;
} 

struct kdheap* create_kdheap(struct kdtree* tree)
{
    struct kdheap* heap = (struct kdheap*) malloc(sizeof(struct kdheap));

    heap->content = (struct kdnear*) malloc(sizeof(struct kdnear) * tree->size);
    heap->length = 0;
    heap->maxsize = tree->size;

    return heap;
}

void free_kdheap(struct kdheap* heap)
{
    if (heap) {
        if (heap->content) free(heap->content);
        free(heap);
    }
}

void kdh_clear(struct kdheap* heap)
{
    heap->length = 0;
}

#define heap_up(i) ( ((i)-1) / 2 )
void kdh_bubbleup(int ix, struct kdheap* heap)
{
    struct kdnear* content = heap->content;
    struct kdnear this = content[ix];

    while (ix > 0) {
        int pix = heap_up(ix);
        if (this.sqrdist < content[pix].sqrdist) {
            struct kdnear parent = content[pix];
            content[pix] = this;
            content[ix] = parent;
            ix = pix;
        } else {
            break;
        }
    }
}

#define heap_left(i) ( 2*(i) + 1 )
#define heap_right(i) ( 2*(i) + 2 )
void kdh_sinkdown(int ix, struct kdheap* heap)
{
    int length = heap->length;
    struct kdnear* content = heap->content;
    if (ix >= length) return;

    struct kdnear this = heap->content[ix];
    int ix_swap = ix;
    while (ix < length) {
        int lix = heap_left(ix);
        int rix = heap_right(ix);

        if (lix < length && content[lix].sqrdist < this.sqrdist) {
            ix_swap = lix;
        }
        if (rix < length && content[rix].sqrdist < content[ix_swap].sqrdist) {
            ix_swap = rix;
        }
        if (ix_swap != ix) {
            struct kdnear shorter = content[ix_swap];
            content[ix_swap] = this;
            content[ix] = shorter;
            ix = ix_swap;
        } else {
            break;
        }
    }
}

void kdh_push(struct kdnear* kdn, struct kdheap* heap)
{
    heap->content[heap->length] = *kdn;
    kdh_bubbleup(heap->length, heap);
    if (heap->length < heap->maxsize) heap->length++;
}

struct point* kdh_pop(struct kdheap* heap)
{
    if (heap->length == 0) return NULL;

    struct kdnear* content = heap->content;
    struct kdnear top = content[0];
    heap->length--;
    if (heap->length > 0) {
        content[0] = content[heap->length];
        kdh_sinkdown(0, heap);
    }
    return top.point;
}

struct kdnear* kdh_look(int ix, struct kdheap* heap)
{
    if (ix > heap->length) return NULL;
    return heap->content + ix;
}

void nearby_points(const struct point* p, struct kdnode* node,
                   struct kdheap* heap, double maxsqrdist)
{
    double sqrdist = metric(p, node->point);

    if (node->left == NULL && node->right == NULL) {
        if (p != node->point && sqrdist < maxsqrdist) {
            struct kdnear this = { .point=node->point, .sqrdist=sqrdist };
            kdh_push(&this, heap);
        }
        return;
    }

    int dim = node->dim;

    struct kdnode* next;
    if (node->right == NULL) {
        next = node->left;
    } else if (node->left == NULL) {
        next = node->right;
    } else {
        if (p->pos[dim] < node->point->pos[dim]) {
            next = node->left;
        } else {
            next = node->right;
        }
    }
    nearby_points(p, next, heap, maxsqrdist);

    if (p != node->point && sqrdist < maxsqrdist) {
        struct kdnear this = { .point=node->point, .sqrdist=sqrdist };
        kdh_push(&this, heap);
    }

    double lin_sqrdist = p->pos[dim] - node->point->pos[dim];
    lin_sqrdist *= lin_sqrdist;

    if (lin_sqrdist < maxsqrdist) {
        if(next == node->left) {
            next = node->right;
        } else {
            next = node->left;
        }
        if (next != NULL) {
            nearby_points(p, next, heap, maxsqrdist);
        }
    }
}

int search_nearby_points(const struct point* p, const struct kdtree* tree,
                         struct kdheap* heap, double maxsqrdist, int maxsize)
{
    if (tree->root == NULL || heap == NULL) return -1;

    kdh_clear(heap);

    if (maxsize > 0 && maxsize <= tree->size) {
        heap->maxsize = maxsize;
    } else {
        heap->maxsize = tree->size;
    }

    if (maxsqrdist < 0) maxsqrdist = DBL_MAX;

    nearby_points(p, tree->root, heap, maxsqrdist);

    return heap->length;
} 
