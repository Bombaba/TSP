#include <stdlib.h>
#include "point.h"
#include "kdtree.h"

void build_tour_nn(struct point pts[], int n_pts, int ixstart,
                   int tour[], const struct kdtree* tree)
{
    struct kdtree* tree_copy;
    if (tree == NULL) {
        tree_copy = build_kdtree(pts, n_pts);
    } else {
        tree_copy = copy_kdtree(tree);
    }

    struct point* start = pts + ixstart;
    struct point* current = start;
    struct point* next;

    int ix = 0;
    tour[ix++] = current->index;
    remove_point_from_tree(current->index, tree_copy);
    while((next = search_nearest(current, tree_copy)) != NULL) {
        current = next;
        tour[ix++] = current->index;
        remove_point_from_tree(current->index, tree_copy);
    }

    free_kdtree(tree_copy);
}
