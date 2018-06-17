#include <stdio.h>
#include "kdtree.h"
#include "point.h"
#include "test.h"

void remove_point_kdtree(struct kdtree* tree)
{
    print_kdtree(tree);
    while(tree->root) {
        int ix;
        printf("Point to remove: ");
        scanf("%d", &ix);
        remove_point_from_tree(ix, tree);
        print_kdtree(tree);
    }
}

void search_circle(struct kdtree* tree)
{
    print_kdtree(tree);
    struct kdheap* heap = create_kdheap(tree);
    while (true) {
        int ix1;
        printf("center point: ");
        scanf("%d", &ix1);
        if (ix1 < 0) break;
        int ix2;
        printf("end point: ");
        scanf("%d", &ix2);
        if (ix2 < 0) break;
        //double radius;
        //printf("radius: ");
        //scanf("%lf", &radius);

        struct point* center = tree->head[ix1].point;
        struct point* end = tree->head[ix2].point;

        printf("Radius Squared: %lf\n", metric(center, end));

        search_nearby_points(center, tree, heap, metric(center, end), -1);

        while (heap->length) {
            struct point* p = kdh_pop(heap);
            printf("#%d : %lf\n", p->index, metric(center, p));
        }
    }

    free_kdheap(heap);
}


