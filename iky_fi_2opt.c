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

#define MAX_N 10000   // 点の数の最大値
#define INF 100000000 // 無限大の定義
#define EPSILON 0.00000001 //ε 小さい正の値

struct point {
    int index;
    int original_index;
    int x;
    int y;
    int pos[2];
    struct point* next;
    struct point* prev;
};

static inline void swap(int* a, int* b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

static inline void copy_point(struct point* origin, struct point* copy)
{
	copy->index = origin->index;
	copy->original_index = origin->original_index;
	copy->x = origin->x;
	copy->y = origin->y;
	copy->pos[0] = origin->pos[0];
	copy->pos[1] = origin->pos[1];
	copy->next = origin->next;
	copy->prev = origin->prev;
}

static inline void insert(struct point* a, struct point* b, struct point* c)
{
    assert( a->next == b && b->prev == a);
    a->next = c;
    c->next = b;
    c->prev = a;
    b->prev = c;        
}

static inline double metric(const struct point* p, const struct point* q)
{
    return (p->x - q->x) * (p->x - q->x) + (p->y - q->y) * (p->y - q->y);

}

static inline double dist(struct point p, struct point q)
{
    return sqrt((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y));
}

static inline void build_list_from_tour(struct point pts[], int n_pts, int tour[])
{
    int i;
    for (i = 0; i < n_pts-1; i++) {
        pts[tour[i]].next = &pts[tour[i+1]];
        pts[tour[i+1]].prev = &pts[tour[i]];
    }
    pts[tour[n_pts-1]].next = &pts[tour[0]];
    pts[tour[0]].prev = &pts[tour[n_pts-1]];
}

static inline void build_tour_from_list(struct point* start, int n_pts, int tour[])
{
    struct point* p = start;
    int i = 0;
    do {
        tour[i] = p->index;
        i++;
        p = p->next;
    } while (p != start);
    assert(i == n_pts);
}

static inline int check_list_from_tour(struct point pts[], int n_pts, int tour[])
{
	int i;
	int flag[n_pts];
	struct point *current;
	for(i=0;i<n_pts;i++) flag[i] = 0;

	for(i=0;i<n_pts;i++) {
		if(tour[i] < 0 || tour[i] >= n_pts) {
			printf("%d\n", i);
			return 0;
		}
		flag[tour[i]] += 1;
	}

	for(i=n_pts-1;i>=0;i--)
		if(flag[i] != 1) {
			printf("tour %d\n", i);
			return 0;
		}

	for(i=0;i<n_pts;i++) flag[i] = 0;
	current = &pts[tour[0]];
	for(i=0;i<n_pts;i++) {
		flag[current->index] += 1;
		current = current->next;
	}

	for(i=0;i<n_pts;i++)
		if(flag[i] != 1) {
			printf("list %d\n", i);
			return 0;
		}

	if(current == &pts[tour[0]]) return 1;
	return 0;
}

static inline void shuffle(int *array, int n)
{
    if (n > 1) {
        int i;
        for (i = 0; i < n-1; i++) {
            int j = i + rand() / (RAND_MAX / (n - 1) + 1);
            int temp = array[j];
            array[j] = array[i];
            array[i] = temp;
        }
    }
}

void build_tour_fi_prec(struct point pts[], int n_pts,
                        int prec[], int n_prec, int tour[])
{
    int i;

    // Number of points in tour
    int n_in_tour = 0;

    // Index array of the points in tour
    bool* in_tour = (bool*) malloc(sizeof(bool) * n_pts);
    // There is no point in tour yet.
    for (i = 0; i < n_pts; i++) in_tour[i] = false;

    // Build tour from the points of prec.
    struct point* p = &pts[prec[n_prec-1]];
    p->next = &pts[prec[0]];
    p->next->prev = p;
    in_tour[prec[0]] = true;
    for (i = 1; i < n_prec; i++) {
        p = p->next;
        p->next = &pts[prec[i]];
        p->next->prev = p;
        in_tour[prec[i]] = true;
    }
    // The number of points in tour is the same as the number of points of prec.
    n_in_tour = n_prec;

    // Insert other points into the tour.
    while (n_in_tour < n_pts) {
        double max_dist = 0;
        struct point* a;
        struct point* r;

        // Find the farthest point from the tour.
        for (i = 0; i < n_pts; i++) {
            // Continue if `pts[i]` is alread in the tour.
            if (in_tour[i]) continue;

            // `r_tmp` is a point not in the tour.
            struct point* r_tmp = &pts[i];

            double min_dist = DBL_MAX;
            struct point* a_tmp;
            struct point* p_in_tour = &pts[prec[0]];
            // Find the nearest point in tour `a_tmp` from `r_tmp`.
            do {
                double d = metric(&pts[i], p_in_tour);
                if (d < min_dist) {
                    min_dist = d;
                    a_tmp = p_in_tour;
                }
                p_in_tour = p_in_tour->next;
            } while (p_in_tour != &pts[prec[0]]);

            // Remember `a_tmp` and `r_tmp` if the distance between
            // them is the current farthest.
            if (min_dist > max_dist) {
                max_dist = min_dist;
                a = a_tmp;
                r = r_tmp;
            }
        }
        struct point* b = a->prev;
        struct point* c = a->next;

       // Insert the fartest point from the tour `r` into either edge next to `a`.
        if (metric(r, b) - metric(b, a) < metric(r, c) - metric(a, c)) {
            // b->a ==> b->r->a
            insert(b, a, r);
        } else {
            // a->c ==> a->r->c
            insert(a, c, r);
        }
        in_tour[r->index] = true;
        n_in_tour++;
    }

    free(in_tour);
    build_tour_from_list(&pts[prec[0]], n_pts, tour);
}

bool two_opt_prec(struct point pts[], int n_pts,
                  int prec[], int n_prec, int tour[])
{
    bool success = false;
    int i;
    bool is_in_prec[n_pts];
    for (i = 0; i < n_pts; i++) is_in_prec[i] = false;
    for (i = 0; i < n_prec; i++) {
        is_in_prec[prec[i]] = true;
    }

    for (i = 0; i < n_pts-3; i++) {
        int co_prec = 0;
        int a_ix = tour[i];
        int b_ix = tour[i + 1];
        if (is_in_prec[b_ix]) co_prec++;
        double dist_ab = dist(pts[a_ix], pts[b_ix]);

        int j;
        for (j = i+2; j < n_pts-1; j++){
            int c_ix = tour[j];
            int d_ix = tour[j + 1];
            if (is_in_prec[c_ix]) {
                co_prec++;
                if (co_prec >= 2) break;
            }
            double delta = (dist_ab + dist(pts[c_ix], pts[d_ix]))
                           - (dist(pts[a_ix], pts[c_ix]) + dist(pts[b_ix], pts[d_ix]));

            if (delta > 0) {
                success = true;
                int g = i + 1;
                int h = j;
                while (g < h) {
                    swap(tour + g, tour + h);
                    g++;
                    h--;
                }
                i--;
                break;
            }
        }
    }

    return success;
}

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

void save_tour_if_shortest(struct point pts[], int n_pts, int tour[], double* min_length)
{
    double length = tour_length(pts, n_pts, tour);
    if (length < *min_length) {
        *min_length = length;
        //memcpy(best_tour, tour, sizeof(int) * n_pts);
        sprintf(tourFileName, "tour%08d.dat", ++num);
        write_tour_data(tourFileName, n_pts, tour);
        printf("\n%s: %lf\n", tourFileName, *min_length);
    }
}


int main(int argc, char *argv[])
{
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
    double min_length = DBL_MAX;

    build_tour_fi_prec(pts, n_pts, prec, n_prec, tour);
    //print_tour(pts, n_pts, tour);
    save_tour_if_shortest(pts, n_pts, tour, &min_length);

    bool success = false;
    do {
        success = 0;
        success = two_opt_prec(pts, n_pts, prec, n_prec, tour);
        save_tour_if_shortest(pts, n_pts, tour, &min_length);

    } while(success);

    return EXIT_SUCCESS;
}

