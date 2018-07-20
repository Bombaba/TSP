#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

#include "point.h"
#include "permutation.h"

void perm(int arr[], int n_arr)
{
    int i;
    for (i = n_arr-1; i > 0; i--) {
        if (arr[i-1] >= arr[i]) continue;

        int j = i;
        int k = n_arr-1;

        while (j < k) {
            swap(&arr[j], &arr[k]);
            j++;
            k--;
        }

        for (j = i; j < n_arr; j++) {
            if (arr[j] > arr[i-1]) {
                swap(&arr[i-1], &arr[j]);
                return;
            }
        }

        assert(false);
    }
    
    int j = 0;
    int k = n_arr-1;
    while (j < k) {
        swap(&arr[j], &arr[k]);
        j++;
        k--;
    }
}
            
