#include <stdbool.h>
#include <float.h>
#include <math.h>
#include "kdtree.h"
#include "point.h"
#include "two_opt.h"

bool two_opt(struct point pts[], int n_pts,
             double mul_factor, double add_factor,
             const struct kdtree* tree, struct kdheap* heap)
{
    bool success = false;
    int i;
    for (i = 0; i < n_pts-2; i++) {
        struct point* pa = pts + i;
        struct point* pb = pa->next;
        double dist_ab = sqrt(metric(pa, pb));
        search_nearby_points(pa, tree, heap, dist_ab*mul_factor+add_factor, false);

        double max_delta = 0;
        struct point* pc = NULL;
        struct point* pd = NULL;
        int j;
        for (j = 0; j < heap->length; j++) {
            struct point* p1 = kdh_look(j, heap)->point;
            struct point* p2 = p1->next;

            double delta = (dist_ab + sqrt(metric(p1, p2)))
                - (sqrt(metric(pa, p1)) + sqrt(metric(pb, p2)));
            if (delta > max_delta) {
                success = true;
                max_delta = delta;
                pc = p1;
                pd = p2;
                //break;
            }
        }
        if (pc) {
            struct point* p = pb;
            while (p != pd) {
                struct point* temp = p->next;
                p->next = p->prev;
                p->prev = temp;
                p = temp;
            }
            pa->next = pc;
            pc->prev = pa;
            pb->next = pd;
            pd->prev = pb;

            i--;
        }
    }

    return success;
}

