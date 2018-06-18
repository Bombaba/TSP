#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <math.h>

struct point {
    int index;
    int x;
    int y;
    int pos[2];
    struct point* next;
    struct point* prev;
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
    for (int i = 0; i < n_pts-1; i++) {
        pts[tour[i]].next = &pts[tour[i+1]];
        pts[tour[i+1]].prev = &pts[tour[i]];
    }
    pts[tour[n_pts-1]].next = &pts[tour[0]];
    pts[tour[0]].prev = &pts[tour[n_pts-1]];
}

static inline void shuffle(int *array, int n)
{
    if (n > 1) {
        for (int i = 0; i < n-1; i++) {
            int j = i + rand() / (RAND_MAX / (n - 1) + 1);
            int temp = array[j];
            array[j] = array[i];
            array[i] = temp;
        }
    }
}

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
    assert(dhead);

    for (int i = 0; i < src->size; i++) {
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


void build_tour_nn2(struct point pts[], int n_pts, int ixstart, int tour[],
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
            sqrt(metric(current, nearest))*1.1, false);

        double score_max = -DBL_MAX;
        for (int j = 0; j < heap->length; j++) {
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

bool two_opt(struct point pts[], int n_pts, double radius,
             const struct kdtree* tree, struct kdheap* heap)
{
    bool success = false;
    for (int i = 0; i < n_pts-2; i++) {
        struct point* pa = pts + i;
        struct point* pb = pa->next;
        double dist_ab = sqrt(metric(pa, pb));
        search_nearby_points(pa, tree, heap, dist_ab+radius, false);

        double max_delta = 0;
        struct point* pc = NULL;
        struct point* pd = NULL;
        for (int j = 0; j < heap->length; j++) {
            struct point* p1 = kdh_look(j, heap)->point;
            struct point* p2 = p1->next;

            double delta = (dist_ab + sqrt(metric(p1, p2)))
                - (sqrt(metric(pa, p1)) + sqrt(metric(pb, p2)));
            if (delta > max_delta) {
                success = true;
                max_delta = delta;
                pc = p1;
                pd = p2;
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

bool or_opt(int n, struct point pts[], int n_pts, double radius,
            struct kdtree* tree, struct kdheap* heap)
{
    bool success = false;
    for (int i = 0; i < n_pts; i++) {
        struct point* pa1 = pts + i;
        struct point* pb = pa1->prev;
        double dist_ba = sqrt(metric(pb, pa1));

        struct point* pa2 = pa1;
        for (int j = 0; j < n-1; j++) pa2 = pa2->next;
        struct point* pc = pa2->next;
        double dist_ac = sqrt(metric(pa2, pc));

        double longer = dist_ba > dist_ac ? dist_ba : dist_ac;
        search_nearby_points(pa1, tree, heap, longer+radius, false);

        double max_delta = 0;
        struct point* px = NULL;
        struct point* py = NULL;
        double dist_bac = dist_ba + dist_ac;
        double dist_bc = sqrt(metric(pb, pc));
        for (int j = 0; j < heap->length; j++) {
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


double dist(struct point p, struct point q) { // pとq の間の距離を計算 
  return sqrt((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y));
}

double tour_length(struct point p[MAX_N], int n, int tour[MAX_N]) {
  int i;
  double sum=0.0;
  for(i=0;i<n;i++) sum+=dist(p[tour[i]],p[tour[(i+1)%n]]);
  return sum;// 総距離が関数の戻り値
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
  for (int i = 0; i < n_pts; i++)
  {
    printf("%3d: %5d %5d\n", pts[i].index, pts[i].x, pts[i].y);
  }
}

bool calc_two_opt(double radius, struct point pts[], int n_pts, int tour[],
                  struct kdtree* tree, struct kdheap* heap)
{
    bool success = false;
    build_list_from_tour(pts, n_pts, tour);
    while (two_opt(pts, n_pts, radius, tree, heap)) {
        success = true;
    }
    if (success) {
        struct point* list_tour = pts;
        for (int i = 0; i < n_pts; i++) {
            tour[i] = list_tour->index;
            list_tour = list_tour->next;
        }
    }

    return success;
}

bool calc_or_opt(int len, double radius, struct point pts[], int n_pts, int tour[],
                 struct kdtree* tree, struct kdheap* heap)
{
    bool success = false;
    build_list_from_tour(pts, n_pts, tour);
    while (or_opt(len, pts, n_pts, radius, tree, heap)) {
        success = true;
    }
    if (success) {
        struct point* list_tour = pts;
        for (int i = 0; i < n_pts; i++) {
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
    //print_points(pts, n_pts);
    //putchar('\n');

    int tour[MAX_N];
    double min_length = DBL_MAX;
    struct kdtree* tree = build_kdtree(pts, n_pts);
    struct kdheap* heap = create_kdheap(tree);
    //int best_tour[MAX_N];
    //print_kdtree(tree);

    int opt_max = n_pts < 100 ? n_pts/10 : 10;

    printf("NN 1\n");
    for (int start = 0; start < n_pts; start++) {
        putchar('-');
        fflush(stdout);

        build_tour_nn(pts, n_pts, start, tour, tree);

        bool success;
        int count = 0;
        do {
            success = false;
            success |= calc_two_opt(5*count, pts, n_pts, tour, tree, heap);
            for (int i = 1; i <= opt_max; i++) {
                success |= calc_or_opt(i, 5*count, pts, n_pts, tour, tree, heap);
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

    printf("\nNN 2\n");
    for (int start = 0; start < n_pts; start++) {
        putchar('-');
        fflush(stdout);

        build_tour_nn2(pts, n_pts, start, tour, tree, heap);

        bool success;
        int count = 0;
        do {
            success = false;
            success |= calc_two_opt(5*count, pts, n_pts, tour, tree, heap);
            for (int i = 1; i <= opt_max; i++) {
                success |= calc_or_opt(i, 5*count, pts, n_pts, tour, tree, heap);
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

    free_kdtree(tree);
    free_kdheap(heap);
    printf("\n");

    return EXIT_SUCCESS;
}

