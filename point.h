#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

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
    for (int i = 0; i < n_pts-1; i++) {
        pts[tour[i]].next = &pts[tour[i+1]];
        pts[tour[i+1]].prev = &pts[tour[i]];
    }
    pts[tour[n_pts-1]].next = &pts[tour[0]];
    pts[tour[0]].prev = &pts[tour[n_pts-1]];
}
