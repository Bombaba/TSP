#include <stdlib.h>
#include <assert.h>
#include "least_spanning_tree.h"

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

static inline struct edge pop_edges(struct lst_edges* edges)
{
    assert(edges->current < edges->size);

    return edges->content[edges->current++];
}

struct lst_edges* collect_edges(struct point pts[], int n_pts)
{
    size_t size = (n_pts * (n_pts - 1)) / 2;

    struct lst_edges* edges = (struct lst_edges*) malloc(sizeof(struct lst_edges));
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

void free_edges(struct lst_edges* edges)
{
    if (edges) {
        if (edges->content) free(edges->content);
        free(edges);
    }
}


int* select_edges(struct point pts[], const int n_pts)
{

    struct lst_edges* all_edges = collect_edges(pts, n_pts);

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

void build_lst_path(struct point pts[], int n_pts, int* tour)
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
