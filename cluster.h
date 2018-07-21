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

static inline double get_pathlength(struct point pts[], int path[], int n_path)
{
    int i;
    double len = 0;
    for (i = 0; i < n_path-1; i++) {
        len += distp(pts + path[i], pts + path[i+1]);
    }
    return len;
}

void build_clusters(struct point pts[], int n_pts,
                    int prec[], int n_prec,
                    int out_clusters[], int out_n_clusters[]);

void build_tour_cl(struct point pts[], int n_pts,
                   int prec[], int n_prec,
                   int clusters[], int n_clusters[],
                   int tour[], int seed);
