#include <stdbool.h>
#include <float.h>
#include <math.h>
#include "point.h"
#include "two_opt_prec.h"


bool two_opt_prec(struct point pts[], int n_pts,
                  int prec[], int n_prec, int tour[])
{
    int i, j;
    bool is_in_prec[n_pts];
    for (i = 0; i < n_pts; i++) is_in_prec[i] = false;
    for (i = 0; i < n_prec; i++) is_in_prec[prec[i]] = true;

    bool success = false;
    for (i = 0; i < n_pts-3; i++) {
        int co_prec = 0;
        int a_ix = tour[i];
        int b_ix = tour[i + 1];
        if (is_in_prec[b_ix]) co_prec++;
        double dist_ab = dist(pts[a_ix], pts[b_ix]);

        for (j = i+2; j < n_pts-1; j++){
            int c_ix = tour[j];
            int d_ix = tour[j + 1];
            if (is_in_prec[c_ix]) {
                co_prec++;
                if (co_prec >= 2) break;
            }
            double delta = (dist_ab + dist(pts[c_ix], pts[d_ix]))
                           - (dist(pts[a_ix], pts[c_ix]) + dist(pts[b_ix], pts[d_ix]));

            if (delta > 0) {
                success = true;
                int g = i + 1;
                int h = j;
                while (g < h) {
                    swap(tour + g, tour + h);
                    g++;
                    h--;
                }
                i--;
                break;
            }
        }
    }

    return success;
}

