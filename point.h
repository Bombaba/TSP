#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

struct point {
    int index;
    int x;
    int y;
    int pos[2];
    struct point* next;
    struct point* prev;
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

