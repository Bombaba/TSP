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

    int n_done = 0;
    bool success = false;

    struct point* y = pts;
    while (n_done <= n_pts) {
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
                n_done = 0;
                break;
            }
            a = b;
            b = a->next;
        }
        n_done++;
        y = z;
    }
    
    build_tour_from_list(pts, n_pts, tour);

    return success;
}

