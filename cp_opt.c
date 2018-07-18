#include <math.h>

#include "cp_opt.h"
#include "nn.h"

void extract_cp(struct point pts[], int n_pts, 
				struct point c_pts[], struct point p_pts[], int *c_num_ptr, 
				double mean)
{
	int i, j, c_num, p_num;
	int degree[n_pts];
	double population[n_pts], *ptrs[n_pts], popu_memo[n_pts];
	double thre = mean * 1.5;
	//double thre = 200.0;

	for(i=0;i<n_pts;i++) { population[i] = 1.0; degree[i] = 0; }

	// 昼間人口の処理
	for(i=0;i<n_pts;i++) {
		double dist;
		for(j=0;j<n_pts;j++) {
			if(i==j) continue;
			dist = sqrt(metric(&pts[i], &pts[j]));
			if(dist < thre) {
				ptrs[degree[i]] = &population[j];
				degree[i]++;
			}
		}
		population[i] -= 1.0;
		for(j=0;j<degree[i];j++) {
			*ptrs[j] += 1.0 / degree[i];
		}
	}

	// 中心地の抽出
	c_num = 0; p_num = 0;
	for(i=0;i<n_pts;i++) {
		if(degree[i] == 0)
			population[i] = 1.0;
		//printf("%lf ", population[i]);
		if(population[i] >= 1.0) {
			//printf("%d %d %d\n", num, pts[i].x, pts[i].y);
			copy_point(&pts[i], &c_pts[c_num]);
			c_pts[c_num].index = c_num;
			c_num++;
		} else {
			//copy_point(&pts[i], &p_pts[p_num]);
			// 挿入ソートで昼間人口高い順のp_ptsを作る
			if(p_num == 0) {
				popu_memo[0] = population[i];
				copy_point(&pts[i], &p_pts[0]);
			//} else if(popu_memo[p_num-1] < population[i]) {
			} else if(popu_memo[p_num-1] > population[i]) {
				j = p_num;
				do {
					popu_memo[j] = popu_memo[j-1];	
					copy_point(&p_pts[j-1], &p_pts[j]);
					j--;
				//} while(j > 0 && popu_memo[j-1] < population[i]);
				} while(j > 0 && popu_memo[j-1] > population[i]);
				popu_memo[j] = population[i];
				copy_point(&pts[i], &p_pts[j]);
			} else {
				copy_point(&pts[i], &p_pts[p_num]);
			}
			p_num++;
		}
	}
	*c_num_ptr = c_num;
	//for(i=0;i<p_num;i++) printf("%lf ", popu_memo[i]);
	//printf("\n");
	//for(i=0;i<c_num;i++) printf("%d ", c_pts[i].original_index); printf("\n");
	//for(i=0;i<p_num;i++) printf("%d ", p_pts[i].original_index); printf("\n");
	printf("Extracted %d cities\nReduced %d cities\n", c_num, p_num);
}

void restore_cp(struct point c_pts[], int c_num, int c_tour[], struct point p_pts[],
				int n_pts, int tour[])
{
	int i, p_num = n_pts-c_num;
	double dist12, dist13, dist23, min_dist;
	struct point *cur;
	struct point *min_ptr;
	struct point list[n_pts];
	
	// c_ptsのリスト作成
	copy_point(&c_pts[c_tour[0]], &list[0]);
	list[0].prev = &list[c_num-1];
	for(i=1;i<c_num;i++) {
		copy_point(&c_pts[c_tour[i]], &list[i]);
		list[i-1].next = &list[i];
		list[i].prev = &list[i-1];
		//printf("%d ", list[i].original_index);
	}
	list[c_num-1].next = &list[0];
	//printf("\n");

	// p_ptsの順に挿入
	for(i=0;i<p_num;i++) {
		min_dist = 10000000000;

		cur = &list[0];
		do {
			dist12 = metric(cur, cur->next);
			dist13 = metric(cur, &p_pts[i]);
			dist23 = metric(&p_pts[i], cur->next);
			if(dist13+dist23-dist12 < min_dist) {
				min_ptr = cur;
				min_dist = sqrt(dist13)+sqrt(dist23)-sqrt(dist12);
			}
			cur = cur->next;
		} while(cur != &list[0]);

		copy_point(&p_pts[i], &list[c_num+i]);
		list[c_num+i].prev = min_ptr;
		list[c_num+i].next = min_ptr->next;
		min_ptr->next->prev = &list[c_num+i];
		min_ptr->next = &list[c_num+i];
	}

	cur = &list[0];
	i = 0;
	do {
		tour[i] = cur->original_index;
		//printf("%d ", tour[i]);
		cur = cur->next;
		i++;
	} while(cur != &list[0]);
	//printf("\n");
}
