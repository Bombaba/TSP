#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
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

