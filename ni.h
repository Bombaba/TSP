struct point* build_tour_ni(struct point pts[], int n_pts, int start,
                            const struct kdtree* tree);

void build_tour_ni_prec(struct point pts[], int n_pts,
                        int prec[], int n_prec,
                        int tour[], const struct kdtree* tree);
