#include <stdlib.h>
#include <float.h>
#include "point.h"
#include "kdtree.h"

struct point* build_tour_ni(struct point pts[], int n_pts, int start,
                            const struct kdtree* tree)
{
    struct kdtree* tree_copy;
    if (tree == NULL) {
        tree_copy = build_kdtree(pts, n_pts);
    } else {
        tree_copy = copy_kdtree(tree);
    }
    struct point* tour = pts + start;
    struct point* current = tour;
    remove_point_from_tree(current->index, tree_copy);

    struct point* nearest = search_nearest(tour, tree_copy);
    current->next = current->prev = nearest;
    nearest->next = nearest->prev = current;
    remove_point_from_tree(nearest->index, tree_copy);

    for (int n_list = 2; n_list < n_pts; n_list++) {
        current = tour;
        double min_dist = DBL_MAX;
        struct point* a;
        struct point* r;

        for (int i = 0; i < n_list; i++) {
            nearest = search_nearest(current, tree_copy);
            double d = metric(current, nearest);
            if (d < min_dist) {
                min_dist = d;
                a = current;
                r = nearest;
            }
            current = current->next;
        }

        struct point* b = a->prev;
        struct point* c = a->next;

        if (metric(r, b) - metric(b, a) < metric(r, c) - metric(a, c)) {
            insert(b, a, r);
        } else {
            insert(a, c, r);
        }
        remove_point_from_tree(r->index, tree_copy);
    }

    free_kdtree(tree_copy);

    return tour;
}
