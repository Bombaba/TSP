#pragma once

static inline unsigned factorial(int n)
{
    int i;
    unsigned ret = 1;

    for (i = 2; i <= n; i++) {
        ret *= i;
    }

    return i;
}

void perm(int arr[], int n_arr);
