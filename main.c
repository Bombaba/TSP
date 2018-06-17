#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "kdtree.h"
#include "point.h"
#include "nn.h"
#include "ni.h"
#include "test.h"

#define MAX_N 10000   // 点の数の最大値
#define INF 100000000 // 無限大の定義
#define EPSILON 0.00000001 //ε 小さい正の値


double dist(struct point p, struct point q) { // pとq の間の距離を計算 
  return sqrt((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y));
}

double tour_length(struct point p[MAX_N], int n, int tour[MAX_N]) {
  int i;
  double sum=0.0;
  for(i=0;i<n;i++) sum+=dist(p[tour[i]],p[tour[(i+1)%n]]);
  return sum;// 総距離が関数の戻り値
}

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
      p[i].next = NULL;
  }                                 

  fclose(fp);
}

void write_tour_data(char *filename, int n, int tour[MAX_N]){
  FILE *fp; 
  int i;
 
 // 構築した巡回路をfilenameという名前のファイルに書き出すためにopen
  if((fp=fopen(filename,"wt"))==NULL){ 
    fprintf(stderr,"Error: File %s open failed.\n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(fp,"%d\n",n);
  for(i=0;i<n; i++){
   fprintf(fp,"%d ",tour[i]);
  }
  fprintf(fp,"\n");
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
  //print_points(pts, n_pts);
  //putchar('\n');

  //check_kdtree(pts, n_pts);

  int tour[MAX_N];
  int best_tour[MAX_N];
  double min_length = DBL_MAX;
  struct kdtree* tree = build_kdtree(pts, n_pts);

  int num = 0;
  char tourFileName[20];

  printf("\n********** Nearest Neighbor **********\n");
  for (int start = 0; start < n_pts; start++) {
      struct point* list_tour = build_tour_nn(pts, n_pts, start, tree);
      for (int i = 0; i < n_pts; i++) {
          tour[i] = list_tour->index;
          list_tour = list_tour->next;
      }
      double length = tour_length(pts, n_pts, tour);
      printf("length: %lf", length);
      if (length < min_length) {
          printf(" [*]");
          min_length = length;
          memcpy(best_tour, tour, sizeof(int) * n_pts);
          sprintf(tourFileName, "tour%08d.dat", ++num);
          write_tour_data(tourFileName, n_pts, best_tour);
      }
      putchar('\n');
  }

  printf("\n********** Nearest Insertion **********\n");
  for (int start = 0; start < n_pts; start++) {
      struct point* list_tour = build_tour_ni(pts, n_pts, start, tree);
      for (int i = 0; i < n_pts; i++) {
          tour[i] = list_tour->index;
          list_tour = list_tour->next;
      }
      double length = tour_length(pts, n_pts, tour);
      printf("length: %lf", length);
      if (length < min_length) {
          printf(" [*]");
          min_length = length;
          memcpy(best_tour, tour, sizeof(int) * n_pts);
          sprintf(tourFileName, "tour%08d.dat", ++num);
          write_tour_data(tourFileName, n_pts, best_tour);
      }
      putchar('\n');
  }


  return EXIT_SUCCESS;
}

