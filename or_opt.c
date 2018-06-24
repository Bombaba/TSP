#include <stdbool.h>
#include <math.h>
#include "point.h"
#include "kdtree.h"

bool or_opt(int n, struct point pts[], int n_pts,
            double mul_factor, double add_factor,
            struct kdtree* tree, struct kdheap* heap)
{
    bool success = false;
    int i;
    for (i = 0; i < n_pts; i++) {
        struct point* pa1 = pts + i;
        struct point* pb = pa1->prev;
        double dist_ba = sqrt(metric(pb, pa1));

        struct point* pa2 = pa1;
        int j;
        for (j = 0; j < n-1; j++) pa2 = pa2->next;
        struct point* pc = pa2->next;
        double dist_ac = sqrt(metric(pa2, pc));

        double longer = dist_ba > dist_ac ? dist_ba : dist_ac;
        search_nearby_points(pa1, tree, heap, longer*mul_factor+add_factor, false);

        double max_delta = 0;
        struct point* px = NULL;
        struct point* py = NULL;
        double dist_bac = dist_ba + dist_ac;
        double dist_bc = sqrt(metric(pb, pc));
        for (j = 0; j < heap->length; j++) {
            struct point* p1 = kdh_look(j, heap)->point;
            struct point* p2 = p1->prev;

            struct point* p = pb;
            bool skip = false;
            while (p != pc->next) {
                if (p1 == p || p2 == p) {
                    skip = true;
                    break;
                }
                p = p->next;
            }
            if (skip) continue;

            double dist_orig = dist_bac + sqrt(metric(p1, p2));
            double delta = dist_orig
                - (sqrt(metric(pa1, p1)) + sqrt(metric(pa2, p2)) + dist_bc);
            if (delta > max_delta) {
                success = true;
                max_delta = delta;
                px = p1;
                py = p2;
                //break;
            }
        }
        if (px) {
            struct point* p = pa1;
            while (p != pc) {
                struct point* temp = p->next;
                p->next = p->prev;
                p->prev = temp;
                p = temp;
            }
            pa1->next = px;
            px->prev = pa1;
            py->next = pa2;
            pa2->prev = py;
            pb->next = pc;
            pc->prev = pb;

            //i--;
        }
    }

    return success;
}
