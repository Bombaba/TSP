#pragma once

#include <stdbool.h>

struct point {
    int index;
    int x;
    int y;
    int pos[2];
    struct point* next;
    struct point* prev;
};

