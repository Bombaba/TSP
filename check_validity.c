#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define INF 100000000 // 無限大の定義
#define MAX_N 10000
#define INTERVAL 30

#define YEAR 2016
//#define PROBLEM_DIR "./in_data1"
#define NUM_OF_PROBLEMS 11
#define NUM_OF_GROUPS 13

struct point {
  int x;
  int y;
};

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
  char buff[500];
  int i;

  if ((fp=fopen(filename,"rt")) == NULL) {// 指定ファイルを読み込み用に開く
    fprintf(stderr,"Error: File %s open failed.\n",filename);
    exit(EXIT_FAILURE);
  }   

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // DIMENSION で始まる行に出会う
	&&(strncmp("DIMENSION",buff,9)!=0)) ; // まで読み飛ばす. 
  sscanf(buff,"DIMENSION : %d",np);           // 点の数 *np を読み込む

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // NODE_COORD_SECTIONで始まる
	&&(strncmp("NODE_COORD_SECTION",buff,18)!=0)) ; // 行に出会うまで, 
                                                        // 読み飛ばす. 
  for(i=0;i<*np;i++) {                       // i=0 から i=(*np)-1まで
    if(fgets(buff,sizeof(buff),fp)!=NULL) 
      sscanf(buff,"%*d %d %d",&(p[i].x),&(p[i].y)); // i番目の点の座標を読み込む
  }                                 

  fclose(fp);
}

void copyarray(int dst[MAX_N], int src[MAX_N], int n) {
  int i;
  for(i=0;i<n;i++) dst[i]=src[i];
}

void shiftarray(int k, int arr[MAX_N], int n) {
  int i, tmp[MAX_N];
  copyarray(tmp,arr,n);
  for (i=0;i<n;i++) arr[(n+i-k)%n]=tmp[i];
}

void printarray(int arr[], int n) {
  int i;

  for(i=0;i<n;i++) printf("%d ",arr[i]);
  printf("\n");
}

int check_validity(int tour[MAX_N], int n) {
  int i,j,start=-1;
  int check[MAX_N];
  int order[MAX_N];

  for(i=0;i<n;i++) check[i]=0;
  for(i=0;i<n;i++) check[tour[i]]+=1; //tour に含まれている都市は何個あるか計数

  for(i=0;i<n;i++) 
    if(check[i]!=1) {
      printf("Error: %d appears %d times.\n",i,check[i]);
      return 0;
    }

  return 1;
}

int read_tour_data(char *filename, int tour[MAX_N], int n){
  FILE *fp;
  int i,nn;

  if ((fp=fopen(filename,"rt")) == NULL) {
    printf("Error: file %s does not exits.\n",filename);
       return 0;
  }
  fscanf(fp,"%d",&nn);
  if(nn!=n) {
    printf("Error: Number of cities in tsp file is %d while that in the tour is %d.\n",n,nn);
       return 0;
  }
  for(i=0; i<n; i++) fscanf(fp,"%d",&tour[i]);
  fclose(fp);

  return 1;
}

int main(int argc, char *argv[]) {
  FILE *fp;
  double tl;
  int i,j;
  char filename[60];
  char dat_filename[60];
  char tsp_filename[60];
  int n;
  int filenamelength;
  int tour[MAX_N];
  struct point p[MAX_N];

  for(i=0;i<MAX_N;i++) tour[i]=-1;

  if(argc!=3) {
    fprintf(stderr,"Usage: %s <dat filename> <tsp filename>\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  filenamelength=strlen(argv[1]);

  if(strncmp(argv[1]+(filenamelength-4),".dat",4)) {
    fprintf(stderr,"1st filename must end with dat.\n");
    exit(EXIT_FAILURE);
  }

  filenamelength=strlen(argv[2]);

  if(strncmp(argv[2]+(filenamelength-4),".tsp",4)) {
    fprintf(stderr,"2nd filename must end with tsp.\n");
    exit(EXIT_FAILURE);
  }

  strcpy(dat_filename,argv[1]);
  //  strncpy(filename,dat_filename,(filenamelength-4));
  //filename[filenamelength-4]='\0';

  sprintf(tsp_filename,"%s",argv[2]);

  read_tsp_data(tsp_filename, p,&n);

  
  if(read_tour_data(dat_filename, tour,n)&&check_validity(tour,n)) {
    printf("Valid dat file. ");
    printf("%lf\n",tour_length(p,n,tour));
    exit(EXIT_SUCCESS);
  }else{
    exit(EXIT_FAILURE);
  }

}


