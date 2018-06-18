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
#include "two_opt.h"
#include "or_opt.h"

#define MAX_N 10000   // 点の数の最大値
#define INF 100000000 // 無限大の定義
#define EPSILON 0.00000001 //ε 小さい正の値

int num = 0;
char tourFileName[20];


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
      sscanf(buff,"%*d %d %d", &(p[i].x), &(p[i].y)); // i番目の点の座標を読み込む
      p[i].index = i;
      p[i].pos[0] = p[i].x;
      p[i].pos[1] = p[i].y;
      p[i].next = NULL;
      p[i].prev = NULL;
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

bool calc_two_opt(double radius, struct point pts[], int n_pts, int tour[],
                  struct kdtree* tree, struct kdheap* heap)
{
    bool success = false;
    build_list_from_tour(pts, n_pts, tour);
    while (two_opt(pts, n_pts, radius, tree, heap)) {
        success = true;
    }
    if (success) {
        struct point* list_tour = pts;
        for (int i = 0; i < n_pts; i++) {
            tour[i] = list_tour->index;
            list_tour = list_tour->next;
        }
    }

    return success;
}

bool calc_or_opt(int len, double radius, struct point pts[], int n_pts, int tour[],
                 struct kdtree* tree, struct kdheap* heap)
{
    bool success = false;
    build_list_from_tour(pts, n_pts, tour);
    while (or_opt(len, pts, n_pts, radius, tree, heap)) {
        success = true;
    }
    if (success) {
        struct point* list_tour = pts;
        for (int i = 0; i < n_pts; i++) {
            tour[i] = list_tour->index;
            list_tour = list_tour->next;
        }
    }
    return success;
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

    int tour[MAX_N];
    double min_length = DBL_MAX;
    struct kdtree* tree = build_kdtree(pts, n_pts);
    struct kdheap* heap = create_kdheap(tree);
    //int best_tour[MAX_N];
    //print_kdtree(tree);

    int opt_max = n_pts < 100 ? n_pts/10 : 10;

    if (n_pts < 100) {
        opt_max = n_pts / 10;
    }

    for (int start = 0; start < n_pts; start++) {
        putchar('-');
        fflush(stdout);

        build_tour_nn(pts, n_pts, start, tour, tree);

        bool success;
        int count = 0;
        do {
            success = false;
            success |= calc_two_opt(5*count, pts, n_pts, tour, tree, heap);
            for (int i = 1; i <= opt_max; i++) {
                success |= calc_or_opt(i, 5*count, pts, n_pts, tour, tree, heap);
            }
            count++;
        } while (success);

        double length = tour_length(pts, n_pts, tour);
        if (length < min_length) {
            min_length = length;
            //memcpy(best_tour, tour, sizeof(int) * n_pts);
            sprintf(tourFileName, "tour%08d.dat", ++num);
            write_tour_data(tourFileName, n_pts, tour);
            printf("\n%s: %lf\n", tourFileName, min_length);

        }
    }

    free_kdtree(tree);
    free_kdheap(heap);
    printf("\n");

    return EXIT_SUCCESS;
}

