#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
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

    int ix = 0;
    tour[ix++] = current->index;
    remove_point_from_tree(current->index, tree_copy);

    struct point* next;
    while((next = search_nearest(current, tree_copy)) != NULL) {
        current = next;
        tour[ix++] = current->index;
        remove_point_from_tree(current->index, tree_copy);
    }

    free_kdtree(tree_copy);
}

struct vector2 {
    double x;
    double y;
};

static inline double dot(struct vector2 v1, struct vector2 v2)
{
    return v1.x * v2.x + v1.y + v2.y;
}

static inline struct vector2 get_direction(struct point* from, struct point* to)
{
    double norm = sqrt(metric(from, to));
    assert(norm > 0.0);

    struct vector2 vec = {
        .x = (to->x - from->x) / norm,
        .y = (to->y - from->y) / norm,
    };
    return vec;
}


void build_tour_nn2(struct point pts[], int n_pts, int ixstart, int tour[],
                    const struct kdtree* tree, struct kdheap* heap)
{
    struct kdtree* tree_copy;
    if (tree == NULL) {
        tree_copy = build_kdtree(pts, n_pts);
    } else {
        tree_copy = copy_kdtree(tree);
    }

    int ix = 0;

    struct point* start = pts + ixstart;
    tour[ix++] = start->index;
    remove_point_from_tree(start->index, tree_copy);

    struct point* next = search_nearest(start, tree_copy);
    tour[ix++] = next->index;
    remove_point_from_tree(next->index, tree_copy);

    struct vector2 vec1 = get_direction(start, next);
    struct vector2 vec2;

    struct point* current = next;
    for (; ix < n_pts; ix++) {
        struct point* nearest = search_nearest(current, tree_copy);

        search_nearby_points(
            current, tree_copy, heap,
            sqrt(metric(current, nearest))*1.1, false);

        double score_max = -DBL_MAX;
        for (int j = 0; j < heap->length; j++) {
            struct point* p = kdh_look(j, heap)->point;

            vec2 = get_direction(current, p);
            double score = dot(vec1, vec2);
            if (score > score_max) {
                score_max = score;
                next = p;
            }
        }
        tour[ix] = next->index;
        remove_point_from_tree(next->index, tree_copy);
        current = next;
        vec1 = vec2;
    }

    free_kdtree(tree_copy);
}
