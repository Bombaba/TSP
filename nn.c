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

void build_tour_nn_prec(int start, struct point pts[], int n_pts,
                        int prec[], int n_prec,
                        int tour[], const struct kdtree* tree,
                        struct kdheap* heap)
{
    int i;

    struct kdtree* tree_copy;
    if (tree == NULL) {
        tree_copy = build_kdtree(pts, n_pts);
    } else {
        tree_copy = copy_kdtree(tree);
    }

    int prec_list[n_pts];
    for (i = 0; i < n_pts; i++) prec_list[i] = -1;
    for (i = 0; i < n_prec; i++) {
        prec_list[prec[i]] = prec[(i+1)%n_prec];
    }

    int current_ix = 0;

    tour[current_ix] = start;
    int next_prec = prec_list[start];
    remove_point_from_tree(start, tree_copy);
    struct point* current = pts + start;
    current_ix++;

    while (tree_copy->n_valid) {
        //printf("%d (next %d): ", tree_copy->n_valid, next_prec);
        struct point* nearest;
        struct point* temp = search_nearest(current, tree_copy);
        if (next_prec == -1) {
            //printf("Case 1\n");
            nearest = temp;
            next_prec = prec_list[nearest->index];
            //next_prec = prec_list[nearest->original_index];
        } else if (prec_list[temp->index] == -1) {
        //} else if (prec_list[temp->original_index] == -1) {
            //printf("Case 2\n");
            nearest = temp;
        } else if (next_prec == temp->index) {
        //} else if (next_prec == temp->original_index) {
            //printf("Case 3\n");
            nearest = temp;
            next_prec = prec_list[next_prec];
        } else {
            search_nearby_points(current, tree_copy, heap, -1, -1);
            //printf("Case 4 : ");
            while((temp = kdh_pop(heap)) != NULL) {
                //printf("%d ", temp->original_index);
                if (next_prec == temp->index) {
                //if (next_prec == temp->original_index) {
                    next_prec = prec_list[next_prec];
                    break;
                } else if (prec_list[temp->index] == -1) {
                //} else if (prec_list[temp->original_index] == -1) {
                    break;
                }
            }
            //printf("\n");
            assert(temp != NULL);
            nearest = temp;
        }

        tour[current_ix] = nearest->index;
        remove_point_from_tree(tour[current_ix], tree_copy);
        current = nearest;
        current_ix++;
    }
    assert(current_ix == n_pts);

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


void build_tour_nn2(double factor,
                    struct point pts[], int n_pts, int ixstart, int tour[],
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
            sqrt(metric(current, nearest))*factor, false);

        double score_max = -DBL_MAX;
        int j;
        for (j = 0; j < heap->length; j++) {
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
