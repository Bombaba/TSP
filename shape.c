#include <math.h>

#include "shape.h"

int reduce_map(struct point pts[], int n_pts, struct point copy_pts[], struct point reduced_pts[], double thre)
{
	int i, num, index;
	struct point list[n_pts*2];
	struct point start, goal;
	struct point *current, *current2;
	start.index = -1; goal.index = -2;

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
	start.index = -1;
	goal.index = -2;


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
	struct point *current = start.next;

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
		list[index-2].prev = current->prev;
		list[index-2].next = &list[index-1];
		list[index-1].prev = &list[index-2];
		list[index-1].next = current->next;
		current->prev->next = &list[index-2];
		current->next->prev = &list[index-1];
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

	for(i=0;i<n_pts;i++)
		copy_point(&pts[i], &copy_pts[i]);

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

		copy_pts[i].x += alpha * (s_x / (double)s_m - pts[i].x);
		copy_pts[i].y += alpha * (s_y / (double)s_m - pts[i].y);
		copy_pts[i].pos[0] = copy_pts[i].x;
		copy_pts[i].pos[1] = copy_pts[i].y;
	}
}