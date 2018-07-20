#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include <float.h>

#include "point.h"
#include "kdtree.h"

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

void build_tour_cl(struct point pts[], int n_pts, int prec[], int n_prec, int tour[])
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

    int n_tour_filled = 0;

    for (i = 0; i < n_prec; i++) {
        int n_cluster = n_nearby_points[prec[i]] + 1;

        if (n_cluster == 1) {
            tour[n_tour_filled] = prec[i];
            n_tour_filled++;
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

        //printf("#%d : ", t->index);
        //for (j = 0; j < n_cluster; j++) {
        //    printf("%d, ", cluster[j]);
        //}
        //printf("\n");

        memcpy(tour + n_tour_filled, cluster, sizeof(int) * n_cluster);
        n_tour_filled += n_cluster;
    }
    assert(n_tour_filled == n_pts);

    free_kdtree(kdprec);
}

