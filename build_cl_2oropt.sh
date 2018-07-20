#!/bin/bash

gcc -Wall cl_2oropt.c kdtree.c cluster.c permutation.c nn_prec.c two_opt_prec.c or_opt_prec.c -lm -o cl_2oropt.out

