#ifndef FUNCS_H
#define FUNCS_H

#include "str.h"

void compute_p(double, double, double, double, double, double*);

void compute_c_rhoe(double, double, double, double, double, double*);

void compute_c_rhop(double, double, double, double, double, double*);

void compute_rho(double, double, double, double, double, double*);

void compute_e(double, double, double, double, double, double*);

void compute_e_from_prho(double, double, double, double, double, double*);

void compute_pc(double, double, double, double, double, double*, double*);

void cons_2_prim(struct fluid**, struct cons**);

void init(double*, double*, double**);

void advance(double*, double*, double**);

void read_pod(double**, int);

void set_2_zero(int, int, double**);

void from_Q_compute_prim(double*, struct cons**, struct fluid**);

void compute_rhs(struct cons**, struct cons**, struct cons**, double*);

void compute_qmT(double**, double**);


#endif
