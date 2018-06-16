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
    int next_ix;

    remove_point_from_tree(start, tree_copy);
    while((next_ix = search_nearest(current, tree_copy)) != -1) {
        current->next = pts + next_ix;
        current = current->next;
        remove_point_from_tree(next_ix, tree_copy);
    }
    current->next = tour;

    free_kdtree(tree_copy);
    return tour;
}
