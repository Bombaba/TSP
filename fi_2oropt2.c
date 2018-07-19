/*
 * 13.Monney-Men
 * Implemented Farthest-Insertion and 2opt.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <assert.h>

#include "point.h"
#include "kdtree.h"
#include "fi_prec.h"
#include "two_opt_prec.h"
#include "or_opt_prec.h"

#define MAX_N 10000   // 点の数の最大値
#define INF 100000000 // 無限大の定義
#define EPSILON 0.00000001 //ε 小さい正の値


int num = 0;
char tourFileName[20];


double tour_length(struct point p[MAX_N], int n, int tour[MAX_N]) {
    int i;
    double sum=0.0;
    for(i=0;i<n;i++) sum+=dist(p[tour[i]],p[tour[(i+1)%n]]);
    return sum;// 総距離が関数の戻り値
}

void print_prec(int prec[], int n_prec)
{
    printf("prec: ");
    int i;
    for (i = 0; i < n_prec; i++) {
        printf("%4d ", prec[i]);
    }
    printf("\n");
}

void print_tour(struct point pts[], int n_pts, int tour[])
{
    printf("\n*** tour ***\n");
    int i;
    for (i = 0; i < n_pts; i++)
    {
        int n = tour[i];
        printf("%4d  %5d %5d\n", pts[n].index, pts[n].x, pts[n].y);
    }
    printf("***********\n\n");
}

int check_tour(int tour[], int n_pts) {
    int i;
    int flag[n_pts];
    for(i=0;i<n_pts;i++) flag[i] = 0;

    for(i=0;i<n_pts;i++) {
            if(tour[i] < 0 || tour[i] >= n_pts)
                    return 0;
            flag[tour[i]] += 1;
    }

    for(i=0;i<n_pts;i++) {
            if(flag[i] != 1) return 0;
    }

    return 1;
}

void read_tsp_data(char *filename, struct point p[MAX_N],int *np, int prec[MAX_N], int *mp) {
  FILE *fp;
  char buff[500];
  int i;

  if ((fp=fopen(filename,"rt")) == NULL) {// 指定ファイルを読み込み用に開く
    fprintf(stderr,"Error: File %s open failed.\n",filename);
    exit(EXIT_FAILURE);
  }   

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // PRECEDENCE_CONSTRAINTS:で始まる行に出会う
	&&(strncmp("PRECEDENCE_CONSTRAINTS:",buff,23)!=0)) ; // まで読み飛ばす. 
  if(strncmp("PRECEDENCE_CONSTRAINTS:",buff,23)==0)  {
    sscanf(buff+24,"%d",mp);
    for(i=0;i<*mp;i++) fscanf(fp,"%d ",&prec[i]);
  } else {
    fprintf(stderr,"Error: There is no precedence constraint in file %s.\n",filename);
    exit(EXIT_FAILURE);
  }

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // DIMENSION で始まる行に出会う
	&&(strncmp("DIMENSION",buff,9)!=0)) ; // まで読み飛ばす. 
      sscanf(buff,"DIMENSION : %d",np);           // 点の数 *np を読み込む

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // NODE_COORD_SECTIONで始まる
          &&(strncmp("NODE_COORD_SECTION",buff,18)!=0)) ; // 行に出会うまで, // 読み飛ばす. 
  for(i=0;i<*np;i++) {                       // i=0 から i=(*np)-1まで
      if(fgets(buff,sizeof(buff),fp)!=NULL) {
          sscanf(buff,"%*d %d %d", &(p[i].x), &(p[i].y)); // i番目の点の座標を読み込む
          p[i].index = i;
          p[i].pos[0] = p[i].x;
          p[i].pos[1] = p[i].y;
          p[i].next = NULL;
          p[i].prev = NULL;
      }
  }                                 

  fclose(fp);
}

void write_tour_data(char *filename, int n, int tour[MAX_N])
{
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
    int i;
    for (i = 0; i < n_pts; i++)
    {
        printf("%3d: %5d %5d\n", pts[i].index, pts[i].x, pts[i].y);
    }
}

void save_tour_if_shortest(struct point pts[], int n_pts, int tour[], int best_tour[], double* min_length)
{
    double length = tour_length(pts, n_pts, tour);
    if (length < *min_length) {
        *min_length = length;
        memcpy(best_tour, tour, sizeof(int) * n_pts);
        sprintf(tourFileName, "tour%08d.dat", ++num);
        write_tour_data(tourFileName, n_pts, tour);
        printf("\n%s: %lf\n", tourFileName, *min_length);
        fflush(stdout);
    }
}


int main(int argc, char *argv[])
{
    int i;
    int n_pts;
    int n_prec;
    struct point pts[MAX_N];
    int prec[MAX_N];   // 順序制約を表現する配列

    if(argc != 2) {
        fprintf(stderr,"Usage: %s <tsp_filename>\n",argv[0]);
        return EXIT_FAILURE;
    }

    read_tsp_data(argv[1], pts, &n_pts, prec, &n_prec);
    //print_prec(prec, n_prec);

    int tour[MAX_N];
    int best_tour[MAX_N];
    double min_length = DBL_MAX;

    struct kdtree* tree = build_kdtree(pts, n_pts);
    struct kdheap* heap = create_kdheap(tree);

    //build_tour_nn_prec(0, pts, n_pts, prec, n_prec, tour, tree, heap);
    //build_tour_ni_prec(pts, n_pts, prec, n_prec, tour, tree);
    build_tour_fi_prec(pts, n_pts, prec, n_prec, tour);
    //print_tour(pts, n_pts, tour);
    save_tour_if_shortest(pts, n_pts, tour, best_tour, &min_length);

    bool success;
    do {
        success = false;
        success |= two_opt_prec(pts, n_pts, prec, n_prec, tour);
        for (i = 0; i < 6; i++) {
            success |= or_opt_prec2(pts, n_pts, prec, n_prec, tour, i);
        }
        save_tour_if_shortest(pts, n_pts, tour, best_tour, &min_length);

    } while(success);
    

    free_kdtree(tree);
    free_kdheap(heap);

    return EXIT_SUCCESS;
}

