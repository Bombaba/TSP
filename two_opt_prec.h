#pragma once

bool two_opt_prec(struct point pts[], int n_pts, int prec[], int n_prec, int tour[]);
bool two_opt_prec2(struct point pts[], int n_pts, int prec[], int n_prec, int tour[]);
bool two_opt_prec3(struct point pts[], int n_pts,
                  int prec[], int n_prec, int tour[],
                  const bool is_in_prec[]);
