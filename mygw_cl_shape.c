/*
 * 13.Monney-Men
 *
 * 制約付き都市を元に都市のクラスタリングを行い，クラスタごとに
 * 部分経路を作成してつなぎ合わせる構築法を実装．
 *
 * また，都市の分布を変形してクラスタリングしやすくするという手法
 * を上記の構築法の前処理として行なっている．分布の変形において，
 * より密な地域に年が少し移動するような処理を行なっている．
 * 
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

    memcpy(tour, clusters, sizeof(int) * n_pts);

    int ix_cluster_begin = 0;
    if (seed != 0) {
        for (i = 0; i < n_prec; i++) {
            shuffle(&tour[ix_cluster_begin+1], n_clusters[i]-2);
            ix_cluster_begin += n_clusters[i];
        }
    }
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

void calc_distribution_param(struct point pts[], int n_pts, struct kdtree *tree, double *mean, double *std, double *skewness, double *kurtosis)
{
    struct point *current, *nearest;
    struct kdtree* tree_copy;
    int i;
    double d1, d2, d3, d4;
    if (tree == NULL) {
        tree_copy = build_kdtree(pts, n_pts);
    } else {
        tree_copy = copy_kdtree(tree);
    }

    d1 = 0.0; d2 = 0.0; d3 = 0.0; d4 = 0.0;
    for(i=0;i<n_pts;i++) {
            double dist;
            current = &pts[i];
            nearest = search_nearest(&pts[i], tree_copy);
            dist = sqrt((current->x-nearest->x)*(current->x-nearest->x)+(current->y-nearest->y)*(current->y-nearest->y));
            d1 += dist;
            d2 += dist * dist;
            d3 += dist * dist * dist;
            d4 += dist * dist * dist * dist;
    }
    d1 = d1 / n_pts; d2 = d2 / n_pts; d3 = d3 / n_pts; d4 = d4 / n_pts;

    *mean = d1;
    *std = sqrt(d2 - d1 * d1);
    *skewness = (d3 - 3 * d1 * d2 + 2 * d1 * d1 * d1) / (*std * *std * *std);
    *kurtosis = (d4 - 4 * d1 * d3 + 6 * d1 * d1 * d2 - 3 * d1 * d1 * d1 * d1) / (*std * *std * *std * *std); 
}

double optimize_thre(struct point pts[], int n_pts, const struct kdtree* tree, double rate)
{
	int i, sum;
	int distribution[1000];
	struct point *current, *nearest;
	struct kdtree* tree_copy;
	//double tmp = (n_pts / 1000) * 0.05 + rate;
	double tmp = rate;
    if (tree == NULL) {
        tree_copy = build_kdtree(pts, n_pts);
    } else {
        tree_copy = copy_kdtree(tree);
    }

	for(i=0;i<1000;i++) distribution[i] = 0;

	for(i=0;i<n_pts;i++) {
		double dist;
		current = &pts[i];
		nearest = search_nearest(&pts[i], tree_copy);
		dist = sqrt((current->x-nearest->x)*(current->x-nearest->x)+(current->y-nearest->y)*(current->y-nearest->y));
		distribution[(int)dist]++;
	}

	sum = distribution[0];
	for(i=1;i<1000;i++) {
		sum += distribution[i];
		printf("%d ", distribution[i]);
		if(sum >= 500) {
			double thre = ((500.0 / n_pts) - (sum - distribution[i])) / distribution[i] + i;
			printf("\nthre : %lf\n", thre);
			return thre;
		} else if((double)sum / n_pts > tmp) {
			double thre = (tmp - (sum - distribution[i])) / distribution[i] + i;
			printf("\nthre : %lf\n", thre);
			return thre;
		}
	}
	return 1000.0;
}

int reduce_map(struct point pts[], int n_pts, struct point copy_pts[], 
				struct point reduced_pts[], const int cp[], int cn, int copy_cp[], double thre)
{
	int i, j, num, index;
	struct point list[n_pts*2];
	struct point start, goal;
	struct point *current, *current2;
	struct kdtree *tree = NULL;
	int weight[n_pts];
	start.index = -1; goal.index = -2;

	for(i=0;i<n_pts;i++) weight[i] = 1;

	if(thre == 0.0) {
		for(i=0;i<n_pts;i++)
			copy_point(&pts[i], &copy_pts[i]);
		return 0;
	}

	if(thre < 0)
		thre = optimize_thre(pts, n_pts, tree, 0.1);

	copy_point(&pts[0], &list[0]);
	list[0].prev = &start; start.next = &list[0];
	for(i=1;i<n_pts;i++) {
		copy_point(&pts[i], &list[i]);
		list[i].prev = &list[i-1];
		list[i-1].next = &list[i];
	}
	list[n_pts-1].next = &goal; goal.prev = &list[n_pts-1];

	num = 0;
	current = start.next;
	do {
		current2 = current->next;
		//printf("c");
		do {
			//printf("a");
			struct point a = *current, b = *current2;
			int success = 0;
			
			for(i=0;i<cn;i++) {
				if(a.original_index == cp[i]) success += 1;
				if(b.original_index == cp[i]) success += 2;
			}
				//printf("%d %d\t", a.index, b.index);
			if(success < 3 && sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)) < thre) {
				//printf("%d + %d\n", current->original_index, current2->original_index);
				// nの操作
				struct point *n = &list[n_pts + num];
				copy_point(current, n);
				//n->x = (current->x+current2->x) / 2.0;
				n->x = (weight[current->original_index]*current->x+weight[current2->original_index]*current2->x) / (double)(weight[current->original_index]+weight[current2->original_index]);
				//n->y = (current->y+current2->y) / 2.0;
				n->y = (weight[current->original_index]*current->y+weight[current2->original_index]*current2->y) / (double)(weight[current->original_index]+weight[current2->original_index]);
				n->pos[0] = n->x;
				n->pos[1] = n->y;
				weight[current->original_index] += weight[current2->original_index];
				weight[current2->original_index] += weight[current->original_index];
				
				if(success < 2) {
					// reduced_ptrsに記録
					copy_point(current, &reduced_pts[num*2]);
					copy_point(current2, &reduced_pts[num*2+1]);

				} else {
					// reduced_ptrsに記録
					//printf("back connect\n");
					copy_point(current2, &reduced_pts[num*2]);
					copy_point(current, &reduced_pts[num*2+1]);

					n->index = current2->index;
					n->original_index = current2->original_index;
				}
				current->prev->next = n;
				current->next->prev = n;

				current2->prev->next = current2->next;
				current2->next->prev = current2->prev;
				//printf("%d %d", current2->prev->index, current2->prev->next->index);

				//printf("%d + %d\nrest: %d, reduced: %d\n", current->original_index, current2->original_index, reduced_pts[num*2].index, reduced_pts[num*2+1].index);
				current = n;
				current2 = current2->next;
				num++;
			} else {
				current2 = current2->next;
			}
			//printf("\t%d\n", current2->index);
		} while(current2 != &goal);
		current = current->next;
		//printf("b %d ", current->index);
		fflush(stdout);
	} while(current != &goal && current->next != &goal);
	//printf("\n");

	current = start.next;
	index = 0;
	do {
		copy_pts[index].index = index;
		//copy_pts[index].original_index = current->index;
		copy_pts[index].original_index = current->original_index;
		copy_pts[index].x = current->x;
		copy_pts[index].y = current->y;
		copy_pts[index].pos[0] = current->x;
		copy_pts[index].pos[1] = current->y;
		copy_pts[index].prev = NULL;
		copy_pts[index].next = NULL;
		index++;
		current = current->next;
	} while(current != &goal);

	for(i=0;i<cn;i++) {
		for(j=0;j<n_pts-num;j++) {
			if(cp[i] == copy_pts[j].original_index) break;
		}
		copy_cp[i] = copy_pts[j].index;
		//printf("(%d, %d) ", cp[i], copy_cp[i]);
	}
	//printf("\n");

	printf("%d cities has reduced!\n", num);
	return num;
}

void restore_reduced_tour(struct point pts[], struct point reduced_pts[], int n_rpts, int tour[], int copy_tour[], int n_pts)
{
	int i, index = 0;
	struct point list[n_pts*2];
	struct point start, goal;
	struct point *current;
	start.index = -1;
	goal.index = -2;

	if(n_rpts == 0) {
		for(i=0;i<n_pts;i++)
			copy_tour[i] = tour[i];
		return;
	}


	// リストの作成
	copy_point(&pts[tour[0]], &list[0]);
	list[0].index = list[0].original_index;
	start.next = &list[0];
	list[0].prev = &start;
	index++;
	for(i=1;i<n_pts-n_rpts;i++) {
		copy_point(&pts[tour[i]], &list[index]);
		list[index].index = list[index].original_index;
		list[index].prev = &list[index-1];
		list[index-1].next = &list[index];
		index++;
	}
	list[index-1].next = &goal;
	goal.prev = &list[index-1];

	index = 0;
	current = start.next;

	// リストにreduced_ptsの内容を追加
	index = n_pts-n_rpts;
	for(i=n_rpts-1;i>=0;i--) {
		struct point a = reduced_pts[i*2];
		struct point b = reduced_pts[i*2+1];
		copy_point(&a, &list[index]);
		copy_point(&b, &list[index+1]);
		index += 2;
		current = start.next;
		do {
			if(current->index == list[index-2].index) break;
			current = current->next;
		} while(current != &goal);
		if((current->x-list[index-2].x)*(current->x-list[index-2].x)+(current->y-list[index-2].y)*(current->y-list[index-2].y) < 
				(current->x-list[index-1].x)*(current->x-list[index-1].x)+(current->y-list[index-1].y)*(current->y-list[index-1].y)) {
			list[index-2].prev = current->prev;
			list[index-2].next = &list[index-1];
			list[index-1].prev = &list[index-2];
			list[index-1].next = current->next;
			current->prev->next = &list[index-2];
			current->next->prev = &list[index-1];
		} else {
			list[index-1].prev = current->prev;
			list[index-1].next = &list[index-2];
			list[index-2].prev = &list[index-1];
			list[index-2].next = current->next;
			current->prev->next = &list[index-1];
			current->next->prev = &list[index-2];
		}
	}

	current = start.next;
	i = 0;

	do {
		copy_tour[i++] = current->index;
		current = current->next;
		//printf("%d ", copy_tour[i-1]);
	} while(current != &goal);
	//printf("\n");
}

// 都市の場所情報を整形する
// 元となるpoint型配列ptsを参考に，copy_ptsの中に整形後のデータを保存する
// grv_thre : gravitation threshold <- どの距離の都市まで引力の影響を受けるか
// alpha : 移動量にかかる係数
void shape_map(struct point pts[], int n_pts, struct point copy_pts[], 
			   const int cp[], int cn, double grv_thre, double alpha)
{
	int i, j, n, s_x, s_y, s_m;
	int mass[n_pts];
	struct kdtree *tree = NULL;

	for(i=0;i<n_pts;i++)
		copy_point(&pts[i], &copy_pts[i]);

	if(grv_thre == 0.0 || alpha == 0.0)
		return;

	if(grv_thre < 0) {
		//grv_thre = optimize_thre(pts, n_pts, tree, 0.6);
		double mean, hoge;
		calc_distribution_param(pts, n_pts, tree, &mean, &hoge, &hoge, &hoge);
		grv_thre = (int)mean + 1;
	}

	for(i=0;i<n_pts;i++) {
		mass[i] = 0;
		for(j=0;j<n_pts;j++)
			if(i!=j && sqrt((pts[i].x-pts[j].x)*(pts[i].x-pts[j].x)+(pts[i].y-pts[j].y)*(pts[i].y-pts[j].y)) < grv_thre)
				mass[i] += 1;
	}

	for(i=0;i<n_pts;i++) {
		n = 0;
		s_x = 0;
		s_y = 0;
		s_m = 0;
		for(j=0;j<n_pts;j++) {
			if(i==j) continue;

			double dist = sqrt((pts[i].x-pts[j].x)*(pts[i].x-pts[j].x)+(pts[i].y-pts[j].y)*(pts[i].y-pts[j].y));
			if(dist < grv_thre) {
				s_x += pts[j].x * mass[j];
				s_y += pts[j].y * mass[j];
				s_m += mass[j];
				n++;
			}
		}
		if(n == 0) {
			//printf("n==0\n");
			continue;
		}
		//printf("n=%d\n", n);

		//copy_pts[i].x += alpha * (s_x / (double)s_m - pts[i].x);
		//copy_pts[i].y += alpha * (s_y / (double)s_m - pts[i].y);
		copy_pts[i].x += alpha * (s_x / (double)s_m - pts[i].x) / sqrt(n);
		copy_pts[i].y += alpha * (s_y / (double)s_m - pts[i].y) / sqrt(n);
		copy_pts[i].pos[0] = copy_pts[i].x;
		copy_pts[i].pos[1] = copy_pts[i].y;
	}
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

int save_tour_if_shortest(struct point pts[], int n_pts, int tour[], int best_tour[], double* min_length)
{
	double length = tour_length(pts, n_pts, tour);
	if (length < *min_length) {
		*min_length = length;
		memcpy(best_tour, tour, sizeof(int) * n_pts);
		sprintf(tourFileName, "tour%08d.dat", ++num);
		write_tour_data(tourFileName, n_pts, tour);
		printf("\n%s: %lf\n", tourFileName, *min_length);
		fflush(stdout);
		return 1;
	}
	return 0;
}


int main(int argc, char *argv[])
{
    int i, j;
    int n_pts;
    int n_prec;
    struct point pts[MAX_N];
	struct point shaped_pts[MAX_N];
    int prec[MAX_N];   // 順序制約を表現する配列
	double elite[100];
	int n_elite = 20;
	double ave_array[100];
	double ave, ave_max = -1.0;
	int idx = -1;
	int cnt = 0;

    if(argc != 2) {
        fprintf(stderr,"Usage: %s <tsp_filename>\n",argv[0]);
        return EXIT_FAILURE;
    }

    read_tsp_data(argv[1], pts, &n_pts, prec, &n_prec);
    //print_prec(prec, n_prec);

    int tour[MAX_N];
    int best_tour[MAX_N];
	double min_length = DBL_MAX;

	//for(j=0;j<n_pts;j++) copy_point(&pts[j], &shaped_pts[j]);

	double k;
	for(k=0.02;k<=1.2;k+=0.02) {
		printf("alpha=%lf\n", k);
		int clusters[MAX_N];
		int n_clusters[MAX_N];

		shape_map(pts, n_pts, shaped_pts, prec, n_prec, -1, k);
		build_clusters(shaped_pts, n_pts, prec, n_prec, clusters, n_clusters);

		int seed;
		double len = 0.0;
		for (seed = 0; seed < 100; seed++) {
			printf("-");
			fflush(stdout);

			build_tour_cl(shaped_pts, n_pts, prec, n_prec, clusters, n_clusters, tour, seed);
			int r = save_tour_if_shortest(pts, n_pts, tour, best_tour, &min_length);

			bool success;
			do {
				success = false;
				success |= two_opt_prec(pts, n_pts, prec, n_prec, tour);
				if(r > 0) {
				for (i = 0; i < 6; i++) {
					success |= or_opt_prec2(pts, n_pts, prec, n_prec, tour, i);
				}
				save_tour_if_shortest(pts, n_pts, tour, best_tour, &min_length);
				}

			} while(success);
			len += tour_length(pts, n_pts, tour);
		}

		printf("len = %lf\n", len);
		if(cnt < n_elite) {
			elite[cnt] = k;
			ave_array[cnt] = len;
			if(ave_max < len) {
				ave_max = len;
				idx = cnt;
			}
			cnt++;
		} else if(ave_max > len) {
			elite[idx] = k;
			ave_array[idx] = len;
			int z;
			ave_max = -1;
			for(z=0;z<n_elite;z++) {
				if(ave_max < ave_array[z]) {
					idx = z;
					ave_max = ave_array[z];
				}
			}
		}
	}

	for(i=0;i<n_elite;i++) printf("%f ", (float)(elite[i])); printf("\n");

	int x;
	for(x=0;x<1;x++) {
		//double k;
		int z;

		//for(k=0.1;k<=2.0;k+=0.1) {
		for(z=0;z<=20;z++) {
			k = elite[z];
			printf("alpha=%lf\n", k);
			int clusters[MAX_N];
			int n_clusters[MAX_N];

			shape_map(pts, n_pts, shaped_pts, prec, n_prec, -1, k);
			build_clusters(shaped_pts, n_pts, prec, n_prec, clusters, n_clusters);

			int seed;
			for (seed = 0; seed < 10000; seed++) {
				printf("-");
				fflush(stdout);

				build_tour_cl(shaped_pts, n_pts, prec, n_prec, clusters, n_clusters, tour, seed);
				save_tour_if_shortest(pts, n_pts, tour, best_tour, &min_length);

				bool success;
				do {
					success = false;
					success |= two_opt_prec(pts, n_pts, prec, n_prec, tour);
					int r;
					for (i = 0; i < 6; i++) {
						success |= or_opt_prec2(pts, n_pts, prec, n_prec, tour, i);
						if(i == 2) r = save_tour_if_shortest(pts, n_pts, tour, best_tour, &min_length);
						if(r == 0) break;
					}
					save_tour_if_shortest(pts, n_pts, tour, best_tour, &min_length);

				} while(success);
			}
		}
	}

    return EXIT_SUCCESS;
}
