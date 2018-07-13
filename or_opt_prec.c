#include <stdbool.h>
#include <float.h>
#include <math.h>
#include "point.h"
#include "or_opt_prec.h"

bool or_opt_prec(struct point pts[], int n_pts,
                 int prec[], int n_prec, int tour[])
{
    int i;

    build_list_from_tour(pts, n_pts, tour);

    bool is_in_prec[n_pts];
    for (i = 0; i < n_pts; i++) is_in_prec[i] = false;
    for (i = 0; i < n_prec; i++) is_in_prec[prec[i]] = true;

    int n_skip = 0;
    bool success = false;

    struct point* y = pts;
    while (n_skip <= n_pts) {
        struct point* x = y->prev;
        struct point* z = y->next;

        double dist_xyz = distp(x, y) + distp(y, z);
        double dist_xz = distp(x, z);

        struct point* a = z;
        struct point* b = a->next;

        while (b != x) {
            if (is_in_prec[y->index] && is_in_prec[a->index]) {
                break;
            }

            double dist_ab = distp(a, b);
            double dist_ayb = distp(a, y) + distp(y, b);
            if ( (dist_xyz + dist_ab) > (dist_xz + dist_ayb) ) {
                success = true;
                x->next = z;
                z->prev = x;
                insert(a, b, y);
                n_skip = 0;
                break;
            }
            a = b;
            b = a->next;
        }
        n_skip++;
        y = z;
    }
    
    build_tour_from_list(pts, n_pts, tour);

    return success;
}

bool or_opt_prec2(struct point pts[], int n_pts,
                 int prec[], int n_prec, int tour[], int length)
{
    int i;

    build_list_from_tour(pts, n_pts, tour);

    bool is_in_prec[n_pts];
    for (i = 0; i < n_pts; i++) is_in_prec[i] = false;
    for (i = 0; i < n_prec; i++) is_in_prec[prec[i]] = true;

    int n_skip = 0;
    bool success = false;

    bool have_prec = false;
    struct point* y1 = pts;
        
    while (n_skip <= n_pts) {
        have_prec = false;

        struct point* y2 = y1;
        if (is_in_prec[y2->index]) have_prec = true;

        for (i = 0; i < length; i++) {
            y2 = y2->next;
            if (is_in_prec[y2->index]) have_prec = true;
        }

        struct point* x = y1->prev;
        struct point* z = y2->next;

        double dist_xyz = distp(x, y1) + distp(y2, z);
        double dist_xz = distp(x, z);

        struct point* a = z;
        struct point* b = a->next;

        while (b != x) {
            if (have_prec && is_in_prec[a->index]) {
                break;
            }

            double dist_ab = distp(a, b);
            double dist_ayb = distp(a, y1) + distp(y2, b);
            if ( (dist_xyz + dist_ab) > (dist_xz + dist_ayb) ) {
                success = true;
                x->next = z;
                z->prev = x;
                a->next = y1;
                y1->prev = a;
                y2->next = b;
                b->prev = y2;

                n_skip = -1;
                break;
            }
            a = b;
            b = a->next;
        }
        n_skip++;

        if (n_skip) {
            y1 = y1->next;;
        } else {
            y1 = z;
        }
    }
    
    build_tour_from_list(pts, n_pts, tour);

    return success;
}

