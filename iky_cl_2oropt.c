/*
 * 13.Monney-Men
 * Implemented Farthest-Insertion and 2opt.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <assert.h>


#define MAX_N 10000   // 点の数の最大値
#define INF 100000000 // 無限大の定義
#define EPSILON 0.00000001 //ε 小さい正の値

struct point {
    int index;
    int original_index;
    int x;
    int y;
    int pos[2];
    struct point* next;
    struct point* prev;
};

static inline void swap(int* a, int* b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

static inline void copy_point(struct point* origin, struct point* copy)
{
	copy->index = origin->index;
	copy->original_index = origin->original_index;
	copy->x = origin->x;
	copy->y = origin->y;
	copy->pos[0] = origin->pos[0];
	copy->pos[1] = origin->pos[1];
	copy->next = origin->next;
	copy->prev = origin->prev;
}

static inline void insert(struct point* a, struct point* b, struct point* c)
{
    assert( a->next == b && b->prev == a);
    a->next = c;
    c->next = b;
    c->prev = a;
    b->prev = c;        
}

static inline double metric(const struct point* p, const struct point* q)
{
    return (p->x - q->x) * (p->x - q->x) + (p->y - q->y) * (p->y - q->y);

}

static inline double distp(struct point* p, struct point* q)
{
    return sqrt((p->x - q->x) * (p->x - q->x) + (p->y - q->y) * (p->y - q->y));
}

static inline double dist(struct point p, struct point q)
{
    return sqrt((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y));
}

static inline void build_list_from_tour(struct point pts[], int n_pts, int tour[])
{
    int i;
    for (i = 0; i < n_pts-1; i++) {
        pts[tour[i]].next = &pts[tour[i+1]];
        pts[tour[i+1]].prev = &pts[tour[i]];
    }
    pts[tour[n_pts-1]].next = &pts[tour[0]];
    pts[tour[0]].prev = &pts[tour[n_pts-1]];
}

static inline void build_tour_from_list(struct point* start, int n_pts, int tour[])
{
    struct point* p = start;
    int i = 0;
    do {
        tour[i] = p->index;
        i++;
        p = p->next;
    } while (p != start);
    assert(i == n_pts);
}

static inline int check_list_from_tour(struct point pts[], int n_pts, int tour[])
{
	int i;
	int flag[n_pts];
	struct point *current;
	for(i=0;i<n_pts;i++) flag[i] = 0;

	for(i=0;i<n_pts;i++) {
		if(tour[i] < 0 || tour[i] >= n_pts) {
			printf("%d\n", i);
			return 0;
		}
		flag[tour[i]] += 1;
	}

	for(i=n_pts-1;i>=0;i--)
		if(flag[i] != 1) {
			printf("tour %d\n", i);
			return 0;
		}

	for(i=0;i<n_pts;i++) flag[i] = 0;
	current = &pts[tour[0]];
	for(i=0;i<n_pts;i++) {
		flag[current->index] += 1;
		current = current->next;
	}

	for(i=0;i<n_pts;i++)
		if(flag[i] != 1) {
			printf("list %d\n", i);
			return 0;
		}

	if(current == &pts[tour[0]]) return 1;
	return 0;
}

static inline void shuffle(int array[], int n)
{
    if (n > 1) {
        int i;
        for (i = 0; i < n-1; i++) {
            int j = i + rand() / (RAND_MAX / (n - i) + 1);
            int temp = array[j];
            array[j] = array[i];
            array[i] = temp;
        }
    }
}


struct vec2 {
    double x;
    double y;
};

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
    int n_valid;
    int* pix_to_nix_map;
    int pix_max;
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

static inline double dot(struct vec2 v1, struct vec2 v2)
{
    return v1.x * v2.x + v1.y + v2.y;
}

static inline double L2norm(struct vec2 v)
{
    return sqrt(v.x * v.x + v.y * v.y);
}

static inline double get_pathlength(struct point pts[], int path[], int n_path)
{
    int i;
    double len = 0;
    for (i = 0; i < n_path-1; i++) {
        len += distp(pts + path[i], pts + path[i+1]);
    }
    return len;
}

void ptree(const struct kdnode* node, int depth)
{
    if (node == NULL) return;

    int i;
    for (i = 0; i < depth; i++) {
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
    assert(dhead);

    int i;
    for (i = 0; i < src->size; i++) {
        dhead[i].point = shead[i].point;
        dhead[i].left = shead[i].left ? (dhead + (shead[i].left - shead)) : NULL;
        dhead[i].right = shead[i].right ? (dhead + (shead[i].right - shead)) : NULL;
        dhead[i].parent = shead[i].parent ? (dhead + (shead[i].parent - shead)) : NULL;
        dhead[i].dim = shead[i].dim;
        dhead[i].valid = shead[i].valid;
    }

    int* dmap = (int *) malloc(sizeof(int) * (src->pix_max+1));
    memcpy(dmap, src->pix_to_nix_map, sizeof(int) * (src->pix_max+1));

    struct kdtree* dst = (struct kdtree*) malloc(sizeof(struct kdtree));
    assert(dst);

    dst->head = dhead;
    dst->root = dst->head + src->root->point->index;
    dst->size = src->size;
    dst->n_valid = src->n_valid;
    dst->pix_to_nix_map = dmap;
    dst->pix_max = src->pix_max;
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

struct kdnode* build_subtree(struct kdnode nodes[], int pton_map[],
                             struct point* p_pts[], int n_pts, int dim)
{
    if (n_pts == 0) return NULL;

    if (n_pts == 1) {
        struct kdnode* nd = &nodes[pton_map[p_pts[0]->index]];
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
    //for (int i = mid+1; i < n_pts; i++) {
    //    if (p_pts[i-1]->pos[dim] == p_pts[i]->pos[dim]) mid++;
    //}
    struct kdnode* nd = &nodes[pton_map[p_pts[mid]->index]];

    nd->dim = dim;
    nd->valid = true;
    nd->left = build_subtree(nodes, pton_map, p_pts, mid, 1-dim);
    if (nd->left) nd->left->parent = nd;

    nd->right = build_subtree(nodes, pton_map, p_pts+mid+1, n_pts-mid-1, 1-dim);
    if (nd->right) nd->right->parent = nd;

    return nd;
}

struct kdtree* build_kdtree(struct point pts[], int n_pts)
{
    int i;

    if (n_pts == 0) return NULL;

    struct point** p_pts = (struct point**) malloc(sizeof(struct point*) * n_pts);
    assert(p_pts);
    struct kdnode* head = (struct kdnode*) malloc(sizeof(struct kdnode) * n_pts);
    assert(head);

    int pix_max = -1;
    for (i = 0; i < n_pts; i++) {
        if (pts[i].index > pix_max) {
            pix_max = pts[i].index;
        }
    }
    assert(pix_max > 0);

    int* pix_to_nix_map = (int*) malloc(sizeof(int) * (pix_max+1));
    assert(pix_to_nix_map);
    for (i = 0; i < pix_max+1; i++) {
        pix_to_nix_map[i] = -1; 
    }

    for (i = 0; i < n_pts; i++) {
        p_pts[i] = pts + i;
        head[i].point = pts + i;
        head[i].left = NULL;
        head[i].right = NULL;
        head[i].parent = NULL;
        head[i].dim = -1;
        head[i].valid = false;

        pix_to_nix_map[pts[i].index] = i;
    }
    struct kdnode* root = build_subtree(head, pix_to_nix_map, p_pts, n_pts, 0);
    free(p_pts);

    struct kdtree* tree = (struct kdtree*) malloc(sizeof(struct kdtree));
    assert(tree);
    tree->head = head;
    tree->root = root;
    tree->size = n_pts;
    tree->n_valid = n_pts;
    tree->pix_to_nix_map = pix_to_nix_map;
    tree->pix_max = pix_max;

    return tree;
}

struct kdtree* build_kdtree_from_indices(struct point pts[], int n_pts,
                                         int pointixs[], int n_indices)
{
    int i, j;

    assert(n_indices > 0);
    assert(n_pts >= n_indices);

    struct point** p_pts = (struct point**) malloc(sizeof(struct point*) * n_indices);
    assert(p_pts);
    struct kdnode* head = (struct kdnode*) malloc(sizeof(struct kdnode) * n_indices);
    assert(head);

    int pix_max = -1;
    for (i = 0; i < n_indices; i++) {
        if (pointixs[i] > pix_max) {
            pix_max = pointixs[i];
        }
    }
    assert(pix_max > 0);

    int* pix_to_nix_map = (int*) malloc(sizeof(int) * (pix_max+1));
    assert(pix_to_nix_map);
    for (i = 0; i < pix_max+1; i++) {
        pix_to_nix_map[i] = -1; 
    }

    for (i = 0; i < n_indices; i++) {
        struct point* p;
        for (j = 0; j < n_pts; j++) {
            if (pts[j].index == pointixs[i]) {
                p = pts + j;
                break;
            }
        }
        assert(j < n_pts);

        p_pts[i] = p;
        head[i].point = p;
        head[i].left = NULL;
        head[i].right = NULL;
        head[i].parent = NULL;
        head[i].dim = -1;
        head[i].valid = false;

        pix_to_nix_map[pointixs[i]] = i;
    }
    struct kdnode* root = build_subtree(head, pix_to_nix_map, p_pts, n_indices, 0);
    free(p_pts);

    struct kdtree* tree = (struct kdtree*) malloc(sizeof(struct kdtree));
    assert(tree);
    tree->head = head;
    tree->root = root;
    tree->size = n_indices;
    tree->n_valid = n_indices;
    tree->pix_to_nix_map = pix_to_nix_map;
    tree->pix_max = pix_max;

    return tree;
}



void free_kdtree(struct kdtree* tree)
{
    if (tree->head == NULL) return;

    free(tree->head);
    free(tree->pix_to_nix_map);
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
    if (pindex > tree->pix_max) {
        fprintf(
            stderr,
            "Error: Trying to remove point#%d from the kdtree,\n"
            "       but the maximum index of the tree is #%d\n",
            pindex, tree->pix_max
        );
        exit(EXIT_FAILURE);
    }

    int nodeix = tree->pix_to_nix_map[pindex];

    if (nodeix < 0) {
        fprintf(
            stderr,
            "Error: Trying to remove point#%d from the kdtree,\n"
            "       but point#%d is not included in the tree.\n",
            pindex, pindex
        );
        exit(EXIT_FAILURE);
    }

    struct kdnode* node = tree->head + nodeix;

    if (node->valid == false) {
        fprintf(
            stderr,
            "Warning: Trying to remove point#%d from the kdtree,\n"
            "         but point#%d was already removed.\n",
            pindex, pindex
        );
        exit(EXIT_FAILURE);
    }

    remove_node(node, tree);
    tree->n_valid--;
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
                         struct kdheap* heap, double maxdist, int maxsize)
{
    if (tree->root == NULL || heap == NULL) return -1;

    kdh_clear(heap);

    if (maxsize > 0 && maxsize <= tree->size) {
        heap->maxsize = maxsize;
    } else {
        heap->maxsize = tree->size;
    }

    if (maxdist < 0) maxdist = DBL_MAX;

    nearby_points(p, tree->root, heap, maxdist * maxdist);

    return heap->length;
} 

void build_clusters(struct point pts[], int n_pts,
                    int prec[], int n_prec,
                    int out_clusters[], int out_n_clusters[])
{
    int i, j;

    int nearby_prec[n_pts];
    int n_nearby_points[n_pts];
    bool is_in_prec[n_pts];
    for (i = 0; i < n_pts; i++) {
        nearby_prec[i] = -1;
        n_nearby_points[i] = -1;
        is_in_prec[i] = false;
    }
    for (i = 0; i < n_prec; i++) {
        n_nearby_points[prec[i]] = 0;
        is_in_prec[prec[i]] = true;
    }

    struct kdtree* kdprec = build_kdtree_from_indices(pts, n_pts, prec, n_prec);

    for (i = 0; i < n_pts; i++) {
        if (is_in_prec[pts[i].index]) {
            //printf("#%d is in prec.\n", pts[i].index);
        } else {
            struct point* nearest = search_nearest(pts + i, kdprec);
            nearby_prec[i] = nearest->index;
            n_nearby_points[nearest->index]++;
            //printf("#%d -> #%d\n", pts[i].index, nearest->index);
        }
    }

    int n_clusters_filled = 0;

    for (i = 0; i < n_prec; i++) {
        int n_cluster = n_nearby_points[prec[i]] + 1;

        if (n_cluster == 1) {
            out_clusters[n_clusters_filled] = prec[i];
            n_clusters_filled++;
            out_n_clusters[i] = 1;
            //printf("#%d : %d\n", prec[i], prec[i]);
            continue;
        }
            
        struct point *s = pts + prec[(i-1 + n_prec) % n_prec];
        struct point *t = pts + prec[i];
        struct point *u = pts + prec[(i+1) % n_prec];

        struct vec2 v1 = {.x = (t->x - s->x), .y = (t->y - s->y) };
        double norm_v1 = L2norm(v1);
        v1.x /= norm_v1;
        v1.y /= norm_v1;

        struct vec2 v2 = {.x = (u->x - t->x), .y = (u->y - t->y) };
        double norm_v2 = L2norm(v2);
        v2.x /= norm_v2;
        v2.y /= norm_v2;

        struct vec2 v3 = {.x = (v1.x + v2.x), .y = (v1.y + v2.y)};
        if (L2norm(v3) == 0.0) {
            v3.x = -v2.y;
            v3.y = v2.x;
        }
        //printf("||v3|| = %lf\n", L2norm(v3));
        //assert(L2norm(v3) > 0.0);

        //printf("%4d(%3d): [", prec[i], n_nearby_points[prec[i]]);
        //for (j = 0; j < n_pts; j++) {
        //    if (nearby_prec[j] == prec[i]) printf("%d, ",j); 
        //}
        //printf("]\n");

        int cluster[n_pts];
        int n_left = 0;
        int n_right = 0;
        double dist_left = DBL_MAX;
        double dist_right = DBL_MAX;
        for (j = 0; j < n_pts; j++) {
            if (nearby_prec[j] != t->index) continue;

            struct point *p = pts + j;
            struct vec2 v_tp = {.x = (p->x - t->x), .y = (p->y - t->y)};
            double direction = dot(v_tp, v3);

            if (direction <= 0.0) {
                // Left
                n_left++;
                double dist_sp = metric(s, p);
                if (dist_sp < dist_left) {
                    dist_left = dist_sp;
                    cluster[n_left-1] = cluster[0];
                    cluster[0] = p->index;
                } else {
                    cluster[n_left-1] = p->index;
                }
            } else {
                // Right
                n_right++;
                double dist_up = metric(u, p);
                if (dist_up < dist_right) {
                    dist_right = dist_up;
                    cluster[n_cluster-n_right] = cluster[n_cluster-1];
                    cluster[n_cluster-1] = p->index;
                } else {
                    cluster[n_cluster-n_right] = p->index;
                }
            }
        }
        cluster[n_left] = t->index;
        assert(n_left + n_right == n_cluster-1);

        memcpy(out_clusters + n_clusters_filled, cluster, sizeof(int) * n_cluster);
        n_clusters_filled += n_cluster;

        out_n_clusters[i] = n_cluster;
    }
    assert(n_clusters_filled == n_pts);

    free_kdtree(kdprec);
}

void build_tour_cl(struct point pts[], int n_pts,
                   int prec[], int n_prec,
                   int clusters[], int n_clusters[],
                   int tour[], int seed)
{
    int i;

    if (seed != 0) srand((unsigned) seed);

    int ix_cluster_begin = 0;
    if (seed != 0) {
        for (i = 0; i < n_prec; i++) {
            shuffle(&clusters[ix_cluster_begin+1], n_clusters[i]-2);
            ix_cluster_begin += n_clusters[i];
        }
    }

    memcpy(tour, clusters, sizeof(int) * n_pts);
}



bool two_opt_prec(struct point pts[], int n_pts,
                  int prec[], int n_prec, int tour[])
{
    int i, j;
    bool is_in_prec[n_pts];
    for (i = 0; i < n_pts; i++) is_in_prec[i] = false;
    for (i = 0; i < n_prec; i++) is_in_prec[prec[i]] = true;

    bool success = false;
    for (i = 0; i < n_pts-3; i++) {
        int co_prec = 0;
        int a_ix = tour[i];
        int b_ix = tour[i + 1];
        if (is_in_prec[b_ix]) co_prec++;
        double dist_ab = dist(pts[a_ix], pts[b_ix]);

        for (j = i+2; j < n_pts-1; j++){
            int c_ix = tour[j];
            int d_ix = tour[j + 1];
            if (is_in_prec[c_ix]) {
                co_prec++;
                if (co_prec >= 2) break;
            }
            double delta = (dist_ab + dist(pts[c_ix], pts[d_ix]))
                           - (dist(pts[a_ix], pts[c_ix]) + dist(pts[b_ix], pts[d_ix]));

            if (delta > 0) {
                success = true;
                int g = i + 1;
                int h = j;
                while (g < h) {
                    swap(tour + g, tour + h);
                    g++;
                    h--;
                }
                i--;
                break;
            }
        }
    }

    return success;
}


bool or_opt_prec(struct point pts[], int n_pts,
                 int prec[], int n_prec, int tour[])
{
    int i;

    build_list_from_tour(pts, n_pts, tour);

    bool is_in_prec[n_pts];
    for (i = 0; i < n_pts; i++) is_in_prec[i] = false;
    for (i = 0; i < n_prec; i++) is_in_prec[prec[i]] = true;

    int n_skip = 0;
    bool success = false;

    struct point* y = pts;
    while (n_skip <= n_pts) {
        struct point* x = y->prev;
        struct point* z = y->next;

        double dist_xyz = distp(x, y) + distp(y, z);
        double dist_xz = distp(x, z);

        struct point* a = z;
        struct point* b = a->next;

        while (b != x) {
            if (is_in_prec[y->index] && is_in_prec[a->index]) {
                break;
            }

            double dist_ab = distp(a, b);
            double dist_ayb = distp(a, y) + distp(y, b);
            if ( (dist_xyz + dist_ab) > (dist_xz + dist_ayb) ) {
                success = true;
                x->next = z;
                z->prev = x;
                insert(a, b, y);
                n_skip = 0;
                break;
            }
            a = b;
            b = a->next;
        }
        n_skip++;
        y = z;
    }
    
    build_tour_from_list(pts, n_pts, tour);

    return success;
}

bool or_opt_prec2(struct point pts[], int n_pts,
                  int prec[], int n_prec, int tour[], int length)
{
    int i;

    build_list_from_tour(pts, n_pts, tour);

    bool is_in_prec[n_pts];
    for (i = 0; i < n_pts; i++) is_in_prec[i] = false;
    for (i = 0; i < n_prec; i++) is_in_prec[prec[i]] = true;

    int n_skip = 0;
    bool success = false;

    bool have_prec = false;
    struct point* y1 = pts;
        
    while (n_skip <= n_pts) {
        have_prec = false;

        struct point* y2 = y1;
        if (is_in_prec[y2->index]) have_prec = true;

        for (i = 0; i < length; i++) {
            y2 = y2->next;
            if (is_in_prec[y2->index]) have_prec = true;
        }

        struct point* x = y1->prev;
        struct point* z = y2->next;

        double dist_xyz = distp(x, y1) + distp(y2, z);
        double dist_xz = distp(x, z);

        struct point* a = z;
        struct point* b = a->next;

        while (b != x) {
            if (have_prec && is_in_prec[a->index]) {
                break;
            }

            double dist_ab = distp(a, b);
            double dist_ayb = distp(a, y1) + distp(y2, b);
            if ( (dist_xyz + dist_ab) > (dist_xz + dist_ayb)+0.00000001 ) {
                success = true;
                x->next = z;
                z->prev = x;
                a->next = y1;
                y1->prev = a;
                y2->next = b;
                b->prev = y2;

                n_skip = -1;
                break;
            }
            a = b;
            b = a->next;
        }
        n_skip++;

        if (n_skip) {
            y1 = y1->next;;
        } else {
            y1 = z;
        }
    }
    
    build_tour_from_list(pts, n_pts, tour);

    return success;
}



int num = 0;
char tourFileName[20];


double tour_length(struct point p[MAX_N], int n, int tour[MAX_N]) {
    int i;
    double sum=0.0;
    for(i=0;i<n;i++) sum+=dist(p[tour[i]],p[tour[(i+1)%n]]);
    return sum;// 総距離が関数の戻り値
}

void print_prec(int prec[], int n_prec)
{
    printf("prec: ");
    int i;
    for (i = 0; i < n_prec; i++) {
        printf("%4d ", prec[i]);
    }
    printf("\n");
}

void print_tour(struct point pts[], int n_pts, int tour[])
{
    printf("\n*** tour ***\n");
    int i;
    for (i = 0; i < n_pts; i++)
    {
        int n = tour[i];
        printf("%4d  %5d %5d\n", pts[n].index, pts[n].x, pts[n].y);
    }
    printf("***********\n\n");
}

int check_tour(int tour[], int n_pts) {
    int i;
    int flag[n_pts];
    for(i=0;i<n_pts;i++) flag[i] = 0;

    for(i=0;i<n_pts;i++) {
            if(tour[i] < 0 || tour[i] >= n_pts)
                    return 0;
            flag[tour[i]] += 1;
    }

    for(i=0;i<n_pts;i++) {
            if(flag[i] != 1) return 0;
    }

    return 1;
}

void read_tsp_data(char *filename, struct point p[MAX_N],int *np, int prec[MAX_N], int *mp) {
  FILE *fp;
  char buff[500];
  int i;

  if ((fp=fopen(filename,"rt")) == NULL) {// 指定ファイルを読み込み用に開く
    fprintf(stderr,"Error: File %s open failed.\n",filename);
    exit(EXIT_FAILURE);
  }   

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // PRECEDENCE_CONSTRAINTS:で始まる行に出会う
	&&(strncmp("PRECEDENCE_CONSTRAINTS:",buff,23)!=0)) ; // まで読み飛ばす. 
  if(strncmp("PRECEDENCE_CONSTRAINTS:",buff,23)==0)  {
    sscanf(buff+24,"%d",mp);
    for(i=0;i<*mp;i++) fscanf(fp,"%d ",&prec[i]);
  } else {
    fprintf(stderr,"Error: There is no precedence constraint in file %s.\n",filename);
    exit(EXIT_FAILURE);
  }

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // DIMENSION で始まる行に出会う
	&&(strncmp("DIMENSION",buff,9)!=0)) ; // まで読み飛ばす. 
      sscanf(buff,"DIMENSION : %d",np);           // 点の数 *np を読み込む

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // NODE_COORD_SECTIONで始まる
          &&(strncmp("NODE_COORD_SECTION",buff,18)!=0)) ; // 行に出会うまで, // 読み飛ばす. 
  for(i=0;i<*np;i++) {                       // i=0 から i=(*np)-1まで
      if(fgets(buff,sizeof(buff),fp)!=NULL) {
          sscanf(buff,"%*d %d %d", &(p[i].x), &(p[i].y)); // i番目の点の座標を読み込む
          p[i].index = i;
          p[i].pos[0] = p[i].x;
          p[i].pos[1] = p[i].y;
          p[i].next = NULL;
          p[i].prev = NULL;
      }
  }                                 

  fclose(fp);
}

void write_tour_data(char *filename, int n, int tour[MAX_N])
{
    FILE *fp; 
    int i;

    // 構築した巡回路をfilenameという名前のファイルに書き出すためにopen
    if((fp=fopen(filename,"wt"))==NULL){ 
            fprintf(stderr,"Error: File %s open failed.\n",filename);
            exit(EXIT_FAILURE);
    }
    fprintf(fp,"%d\n",n);
    for(i=0;i<n; i++){
            fprintf(fp,"%d ",tour[i]);
    }
    fprintf(fp,"\n");
    fclose(fp);
}

void print_points(struct point pts[], int n_pts)
{
    int i;
    for (i = 0; i < n_pts; i++)
    {
        printf("%3d: %5d %5d\n", pts[i].index, pts[i].x, pts[i].y);
    }
}

void save_tour_if_shortest(struct point pts[], int n_pts, int tour[], int best_tour[], double* min_length)
{
    double length = tour_length(pts, n_pts, tour);
    if (length < *min_length) {
        *min_length = length;
        memcpy(best_tour, tour, sizeof(int) * n_pts);
        sprintf(tourFileName, "tour%08d.dat", ++num);
        write_tour_data(tourFileName, n_pts, tour);
        printf("\n%s: %lf\n", tourFileName, *min_length);
        fflush(stdout);
    }
}


int main(int argc, char *argv[])
{
    int i;
    int n_pts;
    int n_prec;
    struct point pts[MAX_N];
    int prec[MAX_N];   // 順序制約を表現する配列

    if(argc != 2) {
        fprintf(stderr,"Usage: %s <tsp_filename>\n",argv[0]);
        return EXIT_FAILURE;
    }

    read_tsp_data(argv[1], pts, &n_pts, prec, &n_prec);
    //print_prec(prec, n_prec);

    int tour[MAX_N];
    int best_tour[MAX_N];
    double min_length = DBL_MAX;

    int clusters[MAX_N];
    int n_clusters[MAX_N];

    build_clusters(pts, n_pts, prec, n_prec, clusters, n_clusters);

    int seed;
    for (seed = 0; seed < 1000000; seed++) {
        printf("-");
        fflush(stdout);

        build_tour_cl(pts, n_pts, prec, n_prec, clusters, n_clusters, tour, seed);
        save_tour_if_shortest(pts, n_pts, tour, best_tour, &min_length);

        bool success;
        do {
            success = false;
            success |= two_opt_prec(pts, n_pts, prec, n_prec, tour);
            for (i = 0; i < 6; i++) {
                success |= or_opt_prec2(pts, n_pts, prec, n_prec, tour, i);
            }
            save_tour_if_shortest(pts, n_pts, tour, best_tour, &min_length);

        } while(success);
    }

    return EXIT_SUCCESS;
}


