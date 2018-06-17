#include <stdlib.h>
#include "point.h"
#include "kdtree.h"

struct point* build_tour_nn(struct point pts[], int n_pts, int start,
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
    struct point* next;

    remove_point_from_tree(current->index, tree_copy);
    while((next = search_nearest(current, tree_copy)) != NULL) {
        current->next = next;
        next->prev = current;
        current = next;
        remove_point_from_tree(current->index, tree_copy);
    }
    current->next = tour;
    tour->prev = current;

    free_kdtree(tree_copy);
    return tour;
}
