#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <stdbool.h>
#include "point.h"

void build_tour_fi_prec(struct point pts[], int n_pts,
                        int prec[], int n_prec, int tour[])
{
    int i;

    // Number of points in tour
    int n_in_tour = 0;
    // Index array of the points in tour
    bool in_tour[n_pts];
    // There is no point in tour yet.
    for (i = 0; i < n_pts; i++) in_tour[i] = false;

    // Build tour from the points of prec.
    struct point* p = &pts[prec[n_prec-1]];
    p->next = &pts[prec[0]];
    p->next->prev = p;
    in_tour[prec[0]] = true;
    for (i = 1; i < n_prec; i++) {
        p = p->next;
        p->next = &pts[prec[i]];
        p->next->prev = p;
        in_tour[prec[i]] = true;
    }
    // The number of points in tour is the same as the number of points of prec.
    n_in_tour = n_prec;

    // Insert other points into the tour.
    while (n_in_tour < n_pts) {
        double max_dist = 0;
        struct point* a;
        struct point* r;

        // Find the farthest point from the tour.
        for (i = 0; i < n_pts; i++) {
            // Continue if `pts[i]` is alread in the tour.
            if (in_tour[i]) continue;

            // `r_tmp` is a point not in the tour.
            struct point* r_tmp = &pts[i];

            double min_dist = DBL_MAX;
            struct point* a_tmp;
            struct point* p_in_tour = &pts[prec[0]];
            // Find the nearest point in tour `a_tmp` from `r_tmp`.
            do {
                double d = metric(&pts[i], p_in_tour);
                if (d < min_dist) {
                    min_dist = d;
                    a_tmp = p_in_tour;
                }
                p_in_tour = p_in_tour->next;
            } while (p_in_tour != &pts[prec[0]]);

            // Remember `a_tmp` and `r_tmp` if the distance between
            // them is the current farthest.
            if (min_dist > max_dist) {
                max_dist = min_dist;
                a = a_tmp;
                r = r_tmp;
            }
        }
        struct point* b = a->prev;
        struct point* c = a->next;

       // Insert the fartest point from the tour `r` into either edge next to `a`.
        if (metric(r, b) - metric(b, a) < metric(r, c) - metric(a, c)) {
            // b->a ==> b->r->a
            insert(b, a, r);
        } else {
            // a->c ==> a->r->c
            insert(a, c, r);
        }
        in_tour[r->index] = true;
        n_in_tour++;
    }

    build_tour_from_list(&pts[prec[0]], n_pts, tour);
}
