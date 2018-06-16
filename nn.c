#include "point.h"
#include "kdtree.h"

struct point* build_tour_nn(struct point pts[], int n_pts, int start)
{
    struct kdtree* tree = build_kdtree(pts, n_pts);

    struct point* tour = pts + start;
    struct point* current = tour;
    int next_ix;

    remove_point_from_tree(start, tree);
    while((next_ix = search_nearest(current, tree)) != -1) {
        current->next = pts + next_ix;
        current = current->next;
        remove_point_from_tree(next_ix, tree);
    }
    current->next = tour;

    //free_kdtree(tree);
    return tour;
}
