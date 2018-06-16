#include <stdio.h>
#include "test.h"
#include "kdtree.h"
#include "point.h"

void check_kdtree(struct point pts[], int n_pts)
{
  struct kdtree* tree = build_kdtree(pts, n_pts);
  print_kdtree(tree);
  while(tree->root) {
      int ix;
      printf("Point to remove: ");
      scanf("%d", &ix);
      remove_point_from_tree(ix, tree);
      print_kdtree(tree);
  }

  free_kdtree(tree);
}
