#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_SIZE 1000

int main(int argc, char *argv[])
{
	int i, j, n;
	int flag[MAX_SIZE][MAX_SIZE];

	srand((unsigned)time(NULL));

	n = atoi(argv[1]);
	//scanf("%d%*c", &n);

	for(i=0;i<MAX_SIZE;i++)
		for(j=0;j<MAX_SIZE;j++)
			flag[i][j] = 0;

	printf("NAME : instance%d\n", n);
	printf("COMMENT : %d-city random TSP\n", n);
	printf("TYPE : TSP\n");
	printf("DIMENSION : %d\n", n);
	printf("EDGE_WEIGHT_TYPE : EUC_2D\n");
	printf("NODE_COORD_SECTION\n");

	for(i=0;i<n;i++) {
		int x, y;
		do {
			x = (int)((double)rand() / RAND_MAX * MAX_SIZE);
			y = (int)((double)rand() / RAND_MAX * MAX_SIZE);
			//printf("%d %d\n", x, y); 
		} while(flag[y][x] != 0);
		printf("%d\t%d\t%d\n", i, x, y);
		flag[y][x] = 1;
	}

	printf("EOF");
	
	return 0;
}
