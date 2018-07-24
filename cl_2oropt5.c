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
#include <time.h>

#include "point.h"
#include "kdtree.h"
#include "cluster.h"
#include "two_opt_prec.h"
#include "or_opt_prec.h"
#include "shape.h"

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

bool save_tour_if_shortest(struct point pts[], int n_pts, int tour[], int best_tour[], double* min_length)
{
    double length = tour_length(pts, n_pts, tour);
    if (length < *min_length) {
        *min_length = length;
        memcpy(best_tour, tour, sizeof(int) * n_pts);
        sprintf(tourFileName, "tour%08d.dat", ++num);
        write_tour_data(tourFileName, n_pts, tour);
        printf("\n%s: %lf\n", tourFileName, *min_length);
        fflush(stdout);
        return true;
    }
    return false;
}


int main(int argc, char *argv[])
{
    int i, j;
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
    bool in_prec[MAX_N];
    for (i = 0; i < n_pts; i++) in_prec[i] = false;
    for (i = 0; i < n_prec; i++) in_prec[prec[i]] = true;

    int tour[MAX_N];
    int best_tour[MAX_N];
    double min_length = DBL_MAX;

    int clusters[MAX_N];
    int n_clusters[MAX_N];

    build_clusters(pts, n_pts, prec, n_prec, clusters, n_clusters);

    int seed = 0;
    int found = 0;
    int co_fail = 0;
    time_t t1 = time(NULL);
    time_t t2;
    do {
        printf("-");
        fflush(stdout);

        seed++;
        build_tour_cl(pts, n_pts, prec, n_prec, clusters, n_clusters, tour, seed);

        bool success;

        do {
            success = false;

            success = two_opt_prec3(pts, n_pts, prec, n_prec, tour, in_prec);
            for (i = 0; i < 6; i++) {
                success = or_opt_prec5(pts, n_pts, prec, n_prec, tour, i, in_prec);
            }
        } while(success);

        if (save_tour_if_shortest(pts, n_pts, tour, best_tour, &min_length)) {
            found++ ;
            co_fail = 0;
        } else {
            co_fail++;
        }
        t2 = time(NULL);
    } while (t2 - t1 <  10 * 60);

    memcpy(tour, best_tour, sizeof(int) * n_pts);

    co_fail = 0;
    int n_shuffle = 5;
    while (true) {
        for (i = 0; i < n_shuffle; i++) {
            int co_prec = 0;
            int g = (int) (((double)rand()+1) / ((double)RAND_MAX+2) * n_pts);
            if (in_prec[tour[g]]) co_prec++;
            int h = g;
            do {
                h++;
                if (in_prec[tour[h]]) co_prec++;
            } while (h < n_pts && co_prec < 2);

            shuffle(&tour[g], h - g);
        }

        bool success;
        do {
            success = false;

            success = two_opt_prec3(pts, n_pts, prec, n_prec, tour, in_prec);
            for (i = 0; i < 6; i++) {
                success = or_opt_prec5(pts, n_pts, prec, n_prec, tour, i, in_prec);
            }
        } while(success);

        if (save_tour_if_shortest(pts, n_pts, tour, best_tour, &min_length)) {
            n_shuffle = 5;
            co_fail = 0;
        } else {
            co_fail++;
            if (co_fail % 30 == 0) {
                n_shuffle++;
                if (n_shuffle > n_prec) n_shuffle = n_prec;
                memcpy(tour, best_tour, sizeof(int) * n_pts);
            }
        }
    }

    //free_kdtree(tree);
    //free_kdheap(heap);

    return EXIT_SUCCESS;
}

