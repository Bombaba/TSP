#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>

struct point {
    int index;
    int original_index;
    int x;
    int y;
    int pos[2];
    struct point* next;
    struct point* prev;
};

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

static inline void shuffle(int *array, int n)
{
    if (n > 1) {
        int i;
        for (i = 0; i < n-1; i++) {
            int j = i + rand() / (RAND_MAX / (n - 1) + 1);
            int temp = array[j];
            array[j] = array[i];
            array[i] = temp;
        }
    }
}


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

struct edge {
    int pix;
    int qix;
    double sqrdist;
};

struct mst_edges {
    struct edge* content;
    size_t current;
    size_t size;
};

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

    struct kdtree* dst = (struct kdtree*) malloc(sizeof(struct kdtree));
    assert(dst);

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
    //for (int i = mid+1; i < n_pts; i++) {
    //    if (p_pts[i-1]->pos[dim] == p_pts[i]->pos[dim]) mid++;
    //}
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
    assert(p_pts);
    struct kdnode* head = (struct kdnode*) malloc(sizeof(struct kdnode) * n_pts);
    assert(head);
    
    if (p_pts == NULL || head == NULL) return NULL;

    int i;
    for (i = 0; i < n_pts; i++) {
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
    assert(tree);
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
    }
    assert(node->valid);

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

void build_tour_nn(struct point pts[], int n_pts, int ixstart,
                   int tour[], const struct kdtree* tree)
{
    struct kdtree* tree_copy;
    if (tree == NULL) {
        tree_copy = build_kdtree(pts, n_pts);
    } else {
        tree_copy = copy_kdtree(tree);
    }

    struct point* start = pts + ixstart;
    struct point* current = start;

    int ix = 0;
    tour[ix++] = current->index;
    remove_point_from_tree(current->index, tree_copy);

    struct point* next;
    while((next = search_nearest(current, tree_copy)) != NULL) {
        current = next;
        tour[ix++] = current->index;
        remove_point_from_tree(current->index, tree_copy);
    }

    free_kdtree(tree_copy);
}

struct vector2 {
    double x;
    double y;
};

static inline double dot(struct vector2 v1, struct vector2 v2)
{
    return v1.x * v2.x + v1.y + v2.y;
}

static inline struct vector2 get_direction(struct point* from, struct point* to)
{
    double norm = sqrt(metric(from, to));
    assert(norm > 0.0);

    struct vector2 vec = {
        .x = (to->x - from->x) / norm,
        .y = (to->y - from->y) / norm,
    };
    return vec;
}


void build_tour_nn2(double factor,
                    struct point pts[], int n_pts, int ixstart, int tour[],
                    const struct kdtree* tree, struct kdheap* heap)
{
    struct kdtree* tree_copy;
    if (tree == NULL) {
        tree_copy = build_kdtree(pts, n_pts);
    } else {
        tree_copy = copy_kdtree(tree);
    }

    int ix = 0;

    struct point* start = pts + ixstart;
    tour[ix++] = start->index;
    remove_point_from_tree(start->index, tree_copy);

    struct point* next = search_nearest(start, tree_copy);
    tour[ix++] = next->index;
    remove_point_from_tree(next->index, tree_copy);

    struct vector2 vec1 = get_direction(start, next);
    struct vector2 vec2;

    struct point* current = next;
    for (; ix < n_pts; ix++) {
        struct point* nearest = search_nearest(current, tree_copy);

        search_nearby_points(
            current, tree_copy, heap,
            sqrt(metric(current, nearest))*factor, false);

        double score_max = -DBL_MAX;
        int j;
        for (j = 0; j < heap->length; j++) {
            struct point* p = kdh_look(j, heap)->point;

            vec2 = get_direction(current, p);
            double score = dot(vec1, vec2);
            if (score > score_max) {
                score_max = score;
                next = p;
            }
        }
        tour[ix] = next->index;
        remove_point_from_tree(next->index, tree_copy);
        current = next;
        vec1 = vec2;
    }

    free_kdtree(tree_copy);
}

int* union_init(int size)
{
    int* table = malloc(sizeof(int) * size);
    int i;
    for (i = 0; i < size; i++) {
        table[i] = -1;
    }
    return table;
}

int union_find(int x, int* table)
{
    while(table[x] >= 0) x = table[x];
    return x;
}

void union_union(int x, int y, int* table)
{
    int r1 = union_find(x, table);
    int r2 = union_find(y, table);

    if (r1 != r2) {
        if(table[r1] < table[r2]) {
            table[r2] = r1;
        } else if (table[r1] > table[r2]) {
            table[r1] = r2;
        } else {
            // assert(table[r1] == table[r2]);
            table[r1] -= 1;
            table[r2] = r1;
        }
    }
}

static int key_compare(const void* e1, const void* e2)
{
    double d1 = ((struct edge*)e1)->sqrdist;
    double d2 = ((struct edge*)e2)->sqrdist;

    return (d1 < d2) ? -1 : (d1 > d2) ? 1 : 0;
}

static inline struct edge pop_edges(struct mst_edges* edges)
{
    assert(edges->current < edges->size);

    return edges->content[edges->current++];
}

struct mst_edges* collect_edges(struct point pts[], int n_pts)
{
    size_t size = (n_pts * (n_pts - 1)) / 2;

    struct mst_edges* edges = (struct mst_edges*) malloc(sizeof(struct mst_edges));
    edges->content = (struct edge*) malloc(sizeof(struct edge) * size);
    edges->current = 0;
    edges->size = size;

    size_t ix = 0;
    int i, j;
    for (i = 0; i < n_pts-1; i++) {
        for (j = i+1; j < n_pts; j++) {
            struct edge e;
            e.pix = i;
            e.qix = j;
            e.sqrdist = metric(pts+i, pts+j);
            edges->content[ix] = e;
            ix++;
        }
    }
    assert(ix == size);
    qsort(edges->content, size, sizeof(struct edge), key_compare);

    return edges;
}

void free_edges(struct mst_edges* edges)
{
    if (edges) {
        if (edges->content) free(edges->content);
        free(edges);
    }
}


int* select_edges(struct point pts[], const int n_pts)
{

    struct mst_edges* all_edges = collect_edges(pts, n_pts);

    int* selected = (int*) malloc(sizeof(int) * n_pts * 2);
    int* count = (int *) malloc(sizeof(int) * n_pts);
    int* utable = union_init(n_pts);

    int i = 0;
    for (i = 0; i < n_pts; i++) count[i] = 0;

    i = 0;
    while (i < n_pts) {
        struct edge e = pop_edges(all_edges);
        int pix = e.pix;
        int qix = e.qix;

        if (count[pix] < 2 && count[qix] < 2
             && (union_find(pix, utable) != union_find(qix, utable) || i == n_pts-1)) {
            union_union(pix, qix, utable);
            selected[2*pix+count[pix]] = qix;
            selected[2*qix+count[qix]] = pix;
            count[pix]++;
            count[qix]++;
            i++;
        }
    }

    free(count);
    free(utable);
    free_edges(all_edges);

    return selected;
}

void build_mst_path(struct point pts[], int n_pts, int* tour)
{
    int* selected_edges = select_edges(pts, n_pts);
    tour[0] = 0;
    tour[1] = selected_edges[2*0+0];
    tour[n_pts-1] = selected_edges[2*0+1];
    int i;
    for (i = 1; i < n_pts-1; i++) {
        int ix = tour[i];
        int pix = selected_edges[2*ix+0];
        int qix = selected_edges[2*ix+1];
        if (tour[i-1] == pix) {
            tour[i+1] = qix;
        } else {
            // assert(tour[i-1] == qix);
            tour[i+1] = pix;
        }
    }
    free(selected_edges);
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
		if(sum >= 100) {
			double thre = ((100 / n_pts) - (sum - distribution[i])) / distribution[i] + i;
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

int reduce_map(struct point pts[], int n_pts, struct point copy_pts[], struct point reduced_pts[], double thre)
{
	int i, num, index;
	struct point list[n_pts*2];
	struct point start, goal;
	struct point *current, *current2;
	struct kdtree *tree = NULL;
	start.index = -1; goal.index = -2;

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
				//printf("%d %d\t", a.index, b.index);
			if(sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)) < thre) {
				// reduced_ptrsに記録
				copy_point(current, &reduced_pts[num*2]);
				copy_point(current2, &reduced_pts[num*2+1]);

				// nの操作
				struct point *n = &list[n_pts + num++];
				copy_point(current, n);
				n->x = (current->x+current2->x) / 2.0;
				n->y = (current->y+current2->y) / 2.0;
				n->pos[0] = n->x;
				n->pos[1] = n->y;

				current->prev->next = n;
				current->next->prev = n;

				current2->prev->next = current2->next;
				current2->next->prev = current2->prev;
				//printf("%d %d", current2->prev->index, current2->prev->next->index);


				current = n;
				current2 = current2->next;
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
		copy_pts[index].original_index = current->index;
		copy_pts[index].x = current->x;
		copy_pts[index].y = current->y;
		copy_pts[index].pos[0] = current->x;
		copy_pts[index].pos[1] = current->y;
		copy_pts[index].prev = NULL;
		copy_pts[index].next = NULL;
		index++;
		current = current->next;
	} while(current != &goal);

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
void shape_map(struct point pts[], int n_pts, struct point copy_pts[], double grv_thre, double alpha)
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

bool two_opt(struct point pts[], int n_pts,
             double mul_factor, double add_factor,
             const struct kdtree* tree, struct kdheap* heap)
{
    bool success = false;
    int i;
    for (i = 0; i < n_pts-2; i++) {
        struct point* pa = pts + i;
        struct point* pb = pa->next;
        double dist_ab = sqrt(metric(pa, pb));
        search_nearby_points(pa, tree, heap, dist_ab*mul_factor+add_factor, false);

        double max_delta = 0;
        struct point* pc = NULL;
        struct point* pd = NULL;
        int j;
        for (j = 0; j < heap->length; j++) {
            struct point* p1 = kdh_look(j, heap)->point;
            struct point* p2 = p1->next;

            double delta = (dist_ab + sqrt(metric(p1, p2)))
                - (sqrt(metric(pa, p1)) + sqrt(metric(pb, p2)));
            if (delta > max_delta) {
                success = true;
                max_delta = delta;
                pc = p1;
                pd = p2;
                //break;
            }
        }
        if (pc) {
            struct point* p = pb;
            while (p != pd) {
                struct point* temp = p->next;
                p->next = p->prev;
                p->prev = temp;
                p = temp;
            }
            pa->next = pc;
            pc->prev = pa;
            pb->next = pd;
            pd->prev = pb;

            i--;
        }
    }

    return success;
}

bool or_opt(int n, struct point pts[], int n_pts,
            double mul_factor, double add_factor,
            struct kdtree* tree, struct kdheap* heap)
{
    bool success = false;
    int i;
    for (i = 0; i < n_pts; i++) {
        struct point* pa1 = pts + i;
        struct point* pb = pa1->prev;
        double dist_ba = sqrt(metric(pb, pa1));

        struct point* pa2 = pa1;
        int j;
        for (j = 0; j < n-1; j++) pa2 = pa2->next;
        struct point* pc = pa2->next;
        double dist_ac = sqrt(metric(pa2, pc));

        double longer = dist_ba > dist_ac ? dist_ba : dist_ac;
        search_nearby_points(pa1, tree, heap, longer*mul_factor+add_factor, false);

        double max_delta = 0;
        struct point* px = NULL;
        struct point* py = NULL;
        double dist_bac = dist_ba + dist_ac;
        double dist_bc = sqrt(metric(pb, pc));
        for (j = 0; j < heap->length; j++) {
            struct point* p1 = kdh_look(j, heap)->point;
            struct point* p2 = p1->prev;

            struct point* p = pb;
            bool skip = false;
            while (p != pc->next) {
                if (p1 == p || p2 == p) {
                    skip = true;
                    break;
                }
                p = p->next;
            }
            if (skip) continue;

            double dist_orig = dist_bac + sqrt(metric(p1, p2));
            double delta = dist_orig
                - (sqrt(metric(pa1, p1)) + sqrt(metric(pa2, p2)) + dist_bc);
            if (delta > max_delta) {
                success = true;
                max_delta = delta;
                px = p1;
                py = p2;
                //break;
            }
        }
        if (px) {
            struct point* p = pa1;
            while (p != pc) {
                struct point* temp = p->next;
                p->next = p->prev;
                p->prev = temp;
                p = temp;
            }
            pa1->next = px;
            px->prev = pa1;
            py->next = pa2;
            pa2->prev = py;
            pb->next = pc;
            pc->prev = pb;

            //i--;
        }
    }

    return success;
}

#define MAX_N 10000   // 点の数の最大値
#define INF 100000000 // 無限大の定義
#define EPSILON 0.00000001 //ε 小さい正の値

int num = 0;
char tourFileName[20];


static inline double dist(struct point p, struct point q) { // pとq の間の距離を計算 
	return sqrt((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y));
}

double tour_length(struct point p[MAX_N], int n, int tour[MAX_N]) {
	int i;
	double sum=0.0;
	for(i=0;i<n;i++) sum+=dist(p[tour[i]],p[tour[(i+1)%n]]);
	return sum;// 総距離が関数の戻り値
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

void read_tsp_data(char *filename, struct point p[MAX_N],int *np) {
	FILE *fp;
	char buff[100];
	int i;

	if ((fp=fopen(filename,"rt")) == NULL) {// 指定ファイルを読み込み用に開く
		fprintf(stderr,"Error: File %s open failed.\n",filename);
		exit(0);
	}   

	while((fgets(buff,sizeof(buff),fp)!=NULL)   // DIMENSION で始まる行に出会う
			&&(strncmp("DIMENSION",buff,9)!=0)) ; // まで読み飛ばす. 
	sscanf(buff,"DIMENSION : %d",np);           // 点の数 *np を読み込む

	while((fgets(buff,sizeof(buff),fp)!=NULL)   // NODE_COORD_SECTIONで始まる
			&&(strncmp("NODE_COORD_SECTION",buff,18)!=0)) ; // 行に出会うまで, 
	// 読み飛ばす. 
	for(i=0;i<*np;i++) {                       // i=0 から i=(*np)-1まで
		if(fgets(buff,sizeof(buff),fp)!=NULL) 
			sscanf(buff,"%*d %d %d", &(p[i].x), &(p[i].y)); // i番目の点の座標を読み込む
		p[i].index = i;
		p[i].pos[0] = p[i].x;
		p[i].pos[1] = p[i].y;
		p[i].next = NULL;
		p[i].prev = NULL;
	}                                 

	fclose(fp);
}

void write_tour_data(char *filename, int n, int tour[MAX_N]){
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

bool calc_two_opt(double mul_factor, double add_factor,
		struct point pts[], int n_pts, int tour[],
		struct kdtree* tree, struct kdheap* heap)
{
	bool success = false;
	build_list_from_tour(pts, n_pts, tour);
	if(!check_list_from_tour(pts, n_pts, tour)) {int i; for(i=0;i<n_pts;i++) printf("%d ", tour[i]); printf("\n");}
	//if(!check_list_from_tour(pts, n_pts, tour)) {struct point* cur = &pts[tour[0]]; do{ printf("%d ", cur->index); cur = cur->next; } while(cur != &pts[tour[0]]); printf("\n"); }
	while (two_opt(pts, n_pts, mul_factor, add_factor, tree, heap)) {
		success = true;
	}
	if (success) {
		struct point* list_tour = pts;
		int i;
		for (i = 0; i < n_pts; i++) {
			tour[i] = list_tour->index;
			list_tour = list_tour->next;
		}
	}

	return success;
}

bool calc_or_opt(int len, double mul_factor, double add_factor,
		struct point pts[], int n_pts, int tour[],
		struct kdtree* tree, struct kdheap* heap)
{
	bool success = false;
	build_list_from_tour(pts, n_pts, tour);
	while (or_opt(len, pts, n_pts, mul_factor, add_factor, tree, heap)) {
		success = true;
	}
	if (success) {
		struct point* list_tour = pts;
		int i;
		for (i = 0; i < n_pts; i++) {
			tour[i] = list_tour->index;
			list_tour = list_tour->next;
		}
	}
	return success;
}

int main(int argc, char *argv[])
{
    int n_pts;
    struct point pts[MAX_N];

    if(argc != 2) {
            fprintf(stderr,"Usage: %s <tsp_filename>\n",argv[0]);
            return EXIT_FAILURE;
    }

    read_tsp_data(argv[1], pts, &n_pts);

    int tour[MAX_N];
    //int best_tour[MAX_N];
    double min_length = DBL_MAX;


    struct kdtree* tree = build_kdtree(pts, n_pts);
    struct kdheap* heap = create_kdheap(tree);

    int opt_max = n_pts < 100 ? n_pts/10 : 10;
    int start;

    double afactor = 0.0;
    double mfactor = 1.5;


    printf("\n############# Minimum Spanning Tree ############\n");
    int mst_tour[MAX_N];
    build_mst_path(pts, n_pts, mst_tour);

    for (afactor = 0.0; afactor < 5.05; afactor += 5.0) {
    for (mfactor = 1.5; mfactor > 1.05; mfactor -= 0.1) {
        memcpy(tour, mst_tour, sizeof(int) * n_pts);
        bool success;
        int count = 0;
        do {
            success = false;
            success |= calc_two_opt(mfactor, afactor*count, pts, n_pts, tour, tree, heap);
            int i;
            for (i = 1; i <= opt_max; i++) {
                success |= calc_or_opt(i, mfactor, afactor*count, pts, n_pts, tour, tree, heap);
            }
            count++;
        } while (success);

        double length = tour_length(pts, n_pts, tour);
        if (length < min_length) {
            min_length = length;
            //memcpy(best_tour, tour, sizeof(int) * n_pts);
            sprintf(tourFileName, "tour%08d.dat", ++num);
            write_tour_data(tourFileName, n_pts, tour);
            printf("\n%s: %lf\n", tourFileName, min_length);
        }
    }}

    
    for (afactor = 0.0; afactor < 11.0; afactor += 5.0) {
    for (mfactor = 1.5; mfactor > 1.0; mfactor -= 0.1) {

        printf("\n##############################################\n"
               "     ADD FACTOR = %0.1lf, MUL FACTOR = %0.1lf\n"
               "##############################################\n", afactor, mfactor);

        printf("\n############# Reduce Map ############\n");
        int reduce_num;
        struct kdtree* tree2;
        int copy_tour[MAX_N];
        struct point reduced_pts[MAX_N];
        struct point reduce_memo[MAX_N];

        reduce_num = reduce_map(pts, n_pts, reduced_pts, reduce_memo, -1);
        tree2 = build_kdtree(reduced_pts, n_pts-reduce_num);
        printf("reduced: %d\n", reduce_num);

        for (start = 0; start < n_pts-reduce_num; start++) {
            build_tour_nn(reduced_pts, n_pts-reduce_num, start, copy_tour, tree2);
            //for (int i; i < n_pts-reduce_num; i++) printf("%d ", copy_tour[i]);
            //printf("\n");
            restore_reduced_tour(reduced_pts, reduce_memo, reduce_num, copy_tour, tour, n_pts);
            //for (int i; i < n_pts; i++) printf("%d ", copy_tour[i]);
            //printf("\n");

            putchar('-');
            fflush(stdout);

            bool success;
            int count = 0;
            do {
                success = false;
                success |= calc_two_opt(1, 0, pts, n_pts, tour, tree, heap);
                int i;
                for (i = 1; i <= opt_max; i++) {
                    success |= calc_or_opt(i, mfactor, afactor*count, pts, n_pts, tour, tree, heap);
                }
                count++;
            } while (success);


            double length = tour_length(pts, n_pts, tour);
            if (length < min_length) {
                min_length = length;
                //memcpy(best_tour, tour, sizeof(int) * n_pts);
                sprintf(tourFileName, "tour%08d.dat", ++num);
                write_tour_data(tourFileName, n_pts, tour);
                printf("\n%s: %lf\n", tourFileName, min_length);

            }
        }

        printf("\n############# Shape Map ############\n");
        struct point shaped_pts[MAX_N];
        shape_map(pts, n_pts, shaped_pts, -1, 0.6);

        struct kdtree* tree_shaped = build_kdtree(shaped_pts, n_pts);

        for (start = 0; start < n_pts; start++) {
            putchar('-');
            fflush(stdout);

            build_tour_nn(shaped_pts, n_pts, start, tour, tree_shaped);

            bool success;
            int count = 0;
            do {
                success = false;
                success |= calc_two_opt(mfactor, afactor*count, pts, n_pts, tour, tree, heap);
                int i;
                for (i = 1; i <= opt_max; i++) {
                    success |= calc_or_opt(i, mfactor, afactor*count, pts, n_pts, tour, tree, heap);
                }
                count++;
            } while (success);

            double length = tour_length(pts, n_pts, tour);
            if (length < min_length) {
                min_length = length;
                //memcpy(best_tour, tour, sizeof(int) * n_pts);
                sprintf(tourFileName, "tour%08d.dat", ++num);
                write_tour_data(tourFileName, n_pts, tour);
                printf("\n%s: %lf\n", tourFileName, min_length);
            }
        }

        free_kdtree(tree_shaped);


        printf("\n############ NN 2 ############\n");

        for (start = 0; start < n_pts; start++) {
            putchar('-');
            fflush(stdout);

            build_tour_nn2(1.1, pts, n_pts, start, tour, tree, heap);

            bool success;
            int count = 0;
            do {
                success = false;
                success |= calc_two_opt(mfactor, afactor*count, pts, n_pts, tour, tree, heap);
                int i;
                for (i = 1; i <= opt_max; i++) {
                    success |= calc_or_opt(i, mfactor, afactor*count, pts, n_pts, tour, tree, heap);
                }
                count++;
            } while (success);

            double length = tour_length(pts, n_pts, tour);
            if (length < min_length) {
                min_length = length;
                //memcpy(best_tour, tour, sizeof(int) * n_pts);
                sprintf(tourFileName, "tour%08d.dat", ++num);
                write_tour_data(tourFileName, n_pts, tour);
                printf("\n%s: %lf\n", tourFileName, min_length);

            }
        }

        //printf("\n############ NN 1 ###########\n");

        //for (start = 0; start < n_pts; start++) {
        //    putchar('-');
        //    fflush(stdout);

        //    build_tour_nn(pts, n_pts, start, tour, tree);

        //    bool success;
        //    int count = 0;
        //    do {
        //        success = false;
        //        success |= calc_two_opt(mfactor, afactor*count, pts, n_pts, tour, tree, heap);
        //        int i;
        //        for (i = 1; i <= opt_max; i++) {
        //            success |= calc_or_opt(i, mfactor, afactor*count, pts, n_pts, tour, tree, heap);
        //        }
        //        count++;
        //    } while (success);

        //    double length = tour_length(pts, n_pts, tour);
        //    if (length < min_length) {
        //        min_length = length;
        //        //memcpy(best_tour, tour, sizeof(int) * n_pts);
        //        sprintf(tourFileName, "tour%08d.dat", ++num);
        //        write_tour_data(tourFileName, n_pts, tour);
        //        printf("\n%s: %lf\n", tourFileName, min_length);

        //    }
        //}
    }}

    free_kdtree(tree);
    free_kdheap(heap);
    printf("\n");

    return EXIT_SUCCESS;
}

