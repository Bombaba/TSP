#include <math.h>

#include "shape.h"
#include "kdtree.h"
#include "evaluate.h"

double optimize_thre(struct point pts[], int n_pts, const struct kdtree* tree, double rate)
{
	int i, sum;
	int distribution[1000];
	struct point *current, *nearest;
	struct kdtree* tree_copy;
	double tmp = (n_pts / 1000) * 0.05 + rate;
    if (tree == NULL) {
        tree_copy = build_kdtree(pts, n_pts);
    } else {
        tree_copy = copy_kdtree(tree);
    }

	for(i=0;i<1000;i++) distribution[i] = 0;

	for(i=0;i<n_pts;i++) {
		double dist;
		current = &pts[i];
		nearest = search_nearest(&pts[i], tree_copy);
		dist = sqrt((current->x-nearest->x)*(current->x-nearest->x)+(current->y-nearest->y)*(current->y-nearest->y));
		distribution[(int)dist]++;
	}

	sum = distribution[0];
	for(i=1;i<1000;i++) {
		sum += distribution[i];
		printf("%d ", distribution[i]);
		if((double)sum / n_pts > tmp || sum >= 500) {
			double thre = (tmp - (sum - distribution[i])) / distribution[i] + i;
			printf("\nthre : %lf\n", thre);
			return thre;
		}
	}
	return 1000.0;
}

int reduce_map(struct point pts[], int n_pts, struct point copy_pts[], struct point reduced_pts[], double thre)
{
	int i, num, index;
	struct point list[n_pts*2];
	struct point start, goal;
	struct point *current, *current2;
	struct kdtree *tree = NULL;
	start.index = -1; goal.index = -2;

	if(thre == 0.0) {
		for(i=0;i<n_pts;i++)
			copy_point(&pts[i], &copy_pts[i]);
		return 0;
	}

	if(thre < 0)
		thre = optimize_thre(pts, n_pts, tree, 0.1);

	copy_point(&pts[0], &list[0]);
	list[0].prev = &start; start.next = &list[0];
	for(i=1;i<n_pts;i++) {
		copy_point(&pts[i], &list[i]);
		list[i].prev = &list[i-1];
		list[i-1].next = &list[i];
	}
	list[n_pts-1].next = &goal; goal.prev = &list[n_pts-1];

	printf("made up list\n");

	num = 0;
	current = start.next;
	do {
		current2 = current->next;
		do {
			struct point a = *current, b = *current2;
			if(sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)) < thre) {
				// reduced_ptrsに記録
				copy_point(current, &reduced_pts[num*2]);
				copy_point(current2, &reduced_pts[num*2+1]);

				// nの操作
				struct point *n = &list[n_pts + num++];
				copy_point(current, n);
				n->x = (current->x+current2->x) / 2.0;
				n->y = (current->y+current2->y) / 2.0;
				n->pos[0] = n->x;
				n->pos[1] = n->y;

				current->prev->next = n;
				current->next->prev = n;

				current2->prev->next = current2->next;
				current2->next->prev = current2->prev;


				current = n;
				current2 = current2->next;
			} else {
				current2 = current2->next;
			}
		} while(current2 != &goal);
		current = current->next;
	} while(current->next != &goal);

	current = start.next;
	index = 0;
	do {
		copy_pts[index].index = index;
		copy_pts[index].original_index = current->index;
		copy_pts[index].x = current->x;
		copy_pts[index].y = current->y;
		copy_pts[index].pos[0] = current->x;
		copy_pts[index].pos[1] = current->y;
		copy_pts[index].prev = NULL;
		copy_pts[index].next = NULL;
		index++;
		current = current->next;
	} while(current != &goal);

	return num;
}

void restore_reduced_tour(struct point pts[], struct point reduced_pts[], int n_rpts, int tour[], int copy_tour[], int n_pts)
{
	int i, index = 0;
	struct point list[n_pts*2];
	struct point start, goal;
	struct point *current;
	start.index = -1;
	goal.index = -2;

	if(n_rpts == 0) {
		for(i=0;i<n_pts;i++)
			copy_tour[i] = tour[i];
		return;
	}


	// リストの作成
	copy_point(&pts[tour[0]], &list[0]);
	list[0].index = list[0].original_index;
	start.next = &list[0];
	list[0].prev = &start;
	index++;
	for(i=1;i<n_pts-n_rpts;i++) {
		copy_point(&pts[tour[i]], &list[index]);
		list[index].index = list[index].original_index;
		list[index].prev = &list[index-1];
		list[index-1].next = &list[index];
		index++;
	}
	list[index-1].next = &goal;
	goal.prev = &list[index-1];

	index = 0;
	current = start.next;

	// リストにreduced_ptsの内容を追加
	index = n_pts-n_rpts;
	for(i=n_rpts-1;i>=0;i--) {
		struct point a = reduced_pts[i*2];
		struct point b = reduced_pts[i*2+1];
		copy_point(&a, &list[index]);
		copy_point(&b, &list[index+1]);
		index += 2;
		current = start.next;
		do {
			if(current->index == list[index-2].index) break;
			current = current->next;
		} while(current != &goal);
		if((current->x-list[index-2].x)*(current->x-list[index-2].x)+(current->y-list[index-2].y)*(current->y-list[index-2].y) < 
				(current->x-list[index-1].x)*(current->x-list[index-1].x)+(current->y-list[index-1].y)*(current->y-list[index-1].y)) {
			list[index-2].prev = current->prev;
			list[index-2].next = &list[index-1];
			list[index-1].prev = &list[index-2];
			list[index-1].next = current->next;
			current->prev->next = &list[index-2];
			current->next->prev = &list[index-1];
		} else {
			list[index-1].prev = current->prev;
			list[index-1].next = &list[index-2];
			list[index-2].prev = &list[index-1];
			list[index-2].next = current->next;
			current->prev->next = &list[index-1];
			current->next->prev = &list[index-1];
		}
	}

	current = start.next;
	i = 0;

	do {
		copy_tour[i++] = current->index;
		current = current->next;
	} while(current != &goal);
}

// 都市の場所情報を整形する
// 元となるpoint型配列ptsを参考に，copy_ptsの中に整形後のデータを保存する
// grv_thre : gravitation threshold <- どの距離の都市まで引力の影響を受けるか
// alpha : 移動量にかかる係数
void shape_map(struct point pts[], int n_pts, struct point copy_pts[], double grv_thre, double alpha)
{
	int i, j, n, s_x, s_y, s_m;
	int mass[n_pts];
	struct kdtree *tree = NULL;

	for(i=0;i<n_pts;i++)
		copy_point(&pts[i], &copy_pts[i]);

	if(grv_thre == 0.0 || alpha == 0.0)
		return;

	if(grv_thre < 0) {
		//grv_thre = optimize_thre(pts, n_pts, tree, 0.6);
		double mean, hoge;
		calc_distribution_param(pts, n_pts, tree, &mean, &hoge, &hoge, &hoge);
		grv_thre = (int)mean + 1;
	}

	for(i=0;i<n_pts;i++) {
		mass[i] = 0;
		for(j=0;j<n_pts;j++)
			if(i!=j && sqrt((pts[i].x-pts[j].x)*(pts[i].x-pts[j].x)+(pts[i].y-pts[j].y)*(pts[i].y-pts[j].y)) < grv_thre)
				mass[i] += 1;
	}

	for(i=0;i<n_pts;i++) {
		n = 0;
		s_x = 0;
		s_y = 0;
		s_m = 0;
		for(j=0;j<n_pts;j++) {
			if(i==j) continue;

			double dist = sqrt((pts[i].x-pts[j].x)*(pts[i].x-pts[j].x)+(pts[i].y-pts[j].y)*(pts[i].y-pts[j].y));
			if(dist < grv_thre) {
				s_x += pts[j].x * mass[j];
				s_y += pts[j].y * mass[j];
				s_m += mass[j];
				n++;
			}
		}
		if(n == 0) {
			//printf("n==0\n");
			continue;
		}
		//printf("n=%d\n", n);

		//copy_pts[i].x += alpha * (s_x / (double)s_m - pts[i].x);
		//copy_pts[i].y += alpha * (s_y / (double)s_m - pts[i].y);
		copy_pts[i].x += alpha * (s_x / (double)s_m - pts[i].x) / sqrt(n);
		copy_pts[i].y += alpha * (s_y / (double)s_m - pts[i].y) / sqrt(n);
		copy_pts[i].pos[0] = copy_pts[i].x;
		copy_pts[i].pos[1] = copy_pts[i].y;
	}
}
