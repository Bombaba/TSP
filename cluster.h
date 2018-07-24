#include "point.h"
#pragma once

struct vec2 {
    double x;
    double y;
};

static inline double dot(struct vec2 v1, struct vec2 v2)
{
    return v1.x * v2.x + v1.y + v2.y;
}

static inline double L2norm(struct vec2 v)
{
    return sqrt(v.x * v.x + v.y * v.y);
}

void build_clusters(struct point pts[], int n_pts,
                    int prec[], int n_prec,
                    int out_clusters[], int out_n_clusters[]);

void build_tour_cl(struct point pts[], int n_pts,
                   int prec[], int n_prec,
                   int clusters[], int n_clusters[],
                   int tour[], int seed);
