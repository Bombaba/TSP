#pragma once

#include <stdbool.h>

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
    a->next = c; c->next = b; c->prev = a; b->prev = c;        
}

static inline double metric(const struct point* p, const struct point* q)
{
    return (p->x - q->x) * (p->x - q->x) + (p->y - q->y) * (p->y - q->y);

}
