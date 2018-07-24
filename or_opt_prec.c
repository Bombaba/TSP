#include <stdbool.h>
#include <float.h>
#include <math.h>
#include "point.h"
#include "or_opt_prec.h"


bool or_opt_prec(struct point pts[], int n_pts,
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

        while (b != y1) {
            if (have_prec && is_in_prec[a->index]) {
                break;
            }

            double dist_ab = distp(a, b);
            double dist_ayb = distp(a, y1) + distp(y2, b);
            if ( (dist_xyz + dist_ab) > (dist_xz + dist_ayb)+0.00000001 ) {
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

        struct point* p;
        struct point* q;

        double delta_max = 0;
        while (b != y1) {
            if (have_prec && is_in_prec[a->index]) {
                break;
            }

            double dist_ab = distp(a, b);
            double dist_ayb = distp(a, y1) + distp(y2, b);

            double delta = (dist_xyz + dist_ab) - (dist_xz + dist_ayb);

            if (delta > delta_max) {
                success = true;
                delta_max = delta;
                p = a;
                q = b;
            }

            a = b;
            b = a->next;
        }

        if (delta_max > 0) {
            x->next = z;
            z->prev = x;
            p->next = y1;
            y1->prev = p;
            y2->next = q;
            q->prev = y2;
            y1 = z;
            n_skip = 0;
        } else {
            y1 = y1->next;;
            n_skip++;
        }
    }
    
    build_tour_from_list(pts, n_pts, tour);

    return success;
}

bool or_opt_prec3(struct point pts[], int n_pts,
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

        struct point* p;
        struct point* q;

        double delta_max = 0;
        while (b != y1) {
            if (have_prec && is_in_prec[a->index]) {
                break;
            }

            double dist_ab = distp(a, b);
            double dist_ayb = distp(a, y1) + distp(y2, b);

            double delta = (dist_xyz + dist_ab) - (dist_xz + dist_ayb);

            if (delta > delta_max) {
                success = true;
                delta_max = delta;
                p = a;
                q = b;
            }

            a = b;
            b = a->next;
        }

        if (delta_max > 0) {
            x->next = z;
            z->prev = x;
            p->next = y1;
            y1->prev = p;
            y2->next = q;
            q->prev = y2;
            //y1 = z;
            n_skip = 0;
        } else {
            y1 = y1->next;;
            n_skip++;
        }
    }
    
    build_tour_from_list(pts, n_pts, tour);

    return success;
}


bool or_opt_prec4(struct point pts[], int n_pts,
                  int prec[], int n_prec, int tour[],
                  int length1, int length2)
{
    int i;

    build_list_from_tour(pts, n_pts, tour);

    bool is_in_prec[n_pts];
    for (i = 0; i < n_pts; i++) is_in_prec[i] = false;
    for (i = 0; i < n_prec; i++) is_in_prec[prec[i]] = true;

    int n_skip = 0;
    bool success = false;

    struct point* y1 = pts;
        
    while (n_skip <= n_pts) {
        bool have_prec = false;

        if (is_in_prec[y1->index]) {
            have_prec = true;
        }

        struct point* y2 = y1;
        for (i = 0; i < length1; i++) {
            y2 = y2->next;
            if (is_in_prec[y2->index]) {
                have_prec = true;
            }
        }

        struct point* x = y1->prev;
        struct point* z = y2->next;

        double dist_xyz = distp(x, y1) + distp(y2, z);

        struct point* p;
        struct point* q1;
        struct point* q2;
        struct point* r;

        double delta_max = 0;

        bool skip = false;

        struct point* a = z;
        if (is_in_prec[a->index]) {
            if (have_prec) {
                skip = true;
            } else {
                have_prec = true;
            }
        }

        struct point* b1 = a->next;
        if (have_prec && is_in_prec[b1->index]) {
            skip = true;
        }

        struct point* b2 = b1;
        for (i = 0; i < length2; i++) {
            b2 = b2->next;
            if (have_prec && is_in_prec[b2->index]) {
                b2 = b2->prev;
                //skip = true;
                break;
            }
        }
        struct point* c = b2->next;

        while (!skip && c != x) {
            double dist_abc = distp(a, b1) + distp(b2, c);

            double dist_xbz = distp(x, b1) + distp(b2, z);
            double dist_ayc = distp(a, y1) + distp(y2, c);

            double delta = (dist_xyz + dist_abc) - (dist_xbz + dist_ayc);

            if (delta > delta_max) {
                success = true;
                delta_max = delta;
                p = a;
                q1 = b1;
                q2 = b2;
                r = c;
            }

            a = a->next;
            if (is_in_prec[a->index]) {
                if (have_prec) {
                    skip = true;
                } else {
                    have_prec = true;
                }
            }

            b1 = a->next;
            if (have_prec && is_in_prec[b1->index]) {
                skip = true;
            }

            b2 = b1;
            for (i = 0; i < length2; i++) {
                b2 = b2->next;
                if (have_prec && is_in_prec[b2->index]) {
                    b2 = b2->prev;
                    //skip = true;
                    break;
                }
            }
            c = b2->next;
        }

        if (delta_max > 0) {
            x->next = q1;
            q1->prev = x;
            q2->next = z;
            z->prev = q2;

            p->next = y1;
            y1->prev = p;
            y2->next = r;
            r->prev = y2;

            y1 = q1;
            n_skip = 0;
        } else {
            y1 = y1->next;
            n_skip++;
        }
    }
    
    build_tour_from_list(pts, n_pts, tour);

    return success;
}

bool or_opt_prec5(struct point pts[], int n_pts,
                  int prec[], int n_prec, int tour[], int length,
                  const bool is_in_prec[])
{
    int i;

    build_list_from_tour(pts, n_pts, tour);

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

        struct point* p;
        struct point* q;

        double delta_max = 0;
        while (b != y1) {
            if (have_prec && is_in_prec[a->index]) {
                break;
            }

            double dist_ab = distp(a, b);
            double dist_ayb = distp(a, y1) + distp(y2, b);

            double delta = (dist_xyz + dist_ab) - (dist_xz + dist_ayb);

            if (delta > delta_max) {
                success = true;
                delta_max = delta;
                p = a;
                q = b;
            }

            a = b;
            b = a->next;
        }

        if (delta_max > 0) {
            x->next = z;
            z->prev = x;
            p->next = y1;
            y1->prev = p;
            y2->next = q;
            q->prev = y2;
            //y1 = z;
            n_skip = 0;
        } else {
            y1 = y1->next;;
            n_skip++;
        }
    }
    
    build_tour_from_list(pts, n_pts, tour);

    return success;
}

