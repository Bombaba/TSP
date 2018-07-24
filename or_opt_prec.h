#pragma once

bool or_opt_prec(struct point pts[], int n_pts, int prec[], int n_prec, int tour[], int length);
bool or_opt_prec2(struct point pts[], int n_pts, int prec[], int n_prec, int tour[], int length);
bool or_opt_prec3(struct point pts[], int n_pts, int prec[], int n_prec, int tour[], int length);
bool or_opt_prec4(struct point pts[], int n_pts,
                  int prec[], int n_prec, int tour[],
                  int length1, int length2);
bool or_opt_prec5(struct point pts[], int n_pts,
                  int prec[], int n_prec, int tour[], int length,
                  const bool is_in_prec[]);
