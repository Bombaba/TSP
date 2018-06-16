#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "kdtree.h"
#include "point.h"

#define MAX_N 10000   // 点の数の最大値
#define INF 100000000 // 無限大の定義
#define EPSILON 0.00000001 //ε 小さい正の値


void read_tsp_data(char *filename, struct point p[MAX_N],int *np) {
  FILE *fp;
  char buff[100];
  int i;

  if ((fp=fopen(filename,"rt")) == NULL) {// 指定ファイルを読み込み用に開く
    fprintf(stderr,"Error: File %s open failed.\n",filename);
    exit(0);
  }   

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // DIMENSION で始まる行に出会う
	&&(strncmp("DIMENSION",buff,9)!=0)) ; // まで読み飛ばす. 
  sscanf(buff,"DIMENSION : %d",np);           // 点の数 *np を読み込む

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // NODE_COORD_SECTIONで始まる
	&&(strncmp("NODE_COORD_SECTION",buff,18)!=0)) ; // 行に出会うまで, 
                                                        // 読み飛ばす. 
  for(i=0;i<*np;i++) {                       // i=0 から i=(*np)-1まで
    if(fgets(buff,sizeof(buff),fp)!=NULL) 
      sscanf(buff,"%d %d %d", &(p[i].index), &(p[i].x), &(p[i].y)); // i番目の点の座標を読み込む
      p[i].pos[0] = p[i].x;
      p[i].pos[1] = p[i].y;
      p[i].visited = false;
  }                                 

  fclose(fp);
}

void print_points(struct point pts[], int n_pts)
{
    for (int i = 0; i < n_pts; i++)
    {
        printf("%3d: %5d %5d\n", pts[i].index, pts[i].x, pts[i].y);
    }
}

int main(int argc, char *argv[])
{
  int n_pts;
  struct point pts[MAX_N];

  if(argc != 2) {
    fprintf(stderr,"Usage: %s <tsp_filename>\n",argv[0]);
    return EXIT_FAILURE;
  }

  read_tsp_data(argv[1], pts, &n_pts);
  print_points(pts, n_pts);

  struct kdtree* tree = build_kdtree(pts, n_pts); 
  free_kdtree(tree);
  tree = NULL;

  return EXIT_SUCCESS;
}

