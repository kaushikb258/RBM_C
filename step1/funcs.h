#ifndef FUNCS_H
#define FUNCS_H


#include "str.h"

void init(struct fluid **, struct cons **);
void compute_pc(double, double, double, double, double, double*, double*);
void compute_p(double, double, double, double, double, double*);
void compute_c_rhoe(double, double, double, double, double, double*);
void compute_c_rhop(double, double, double, double, double, double*);
void compute_rho(double, double, double, double, double, double*);
void compute_e(double, double, double, double, double, double*);
void compute_e_from_prho(double, double, double, double, double, double*);

void advance(struct fluid **, struct cons **);
void cons_2_prim(struct fluid **, struct cons **);
void set_dF_zero(struct cons **);
void update(struct cons **, struct cons **, struct cons **, struct cons **);
void write_snapshot(FILE*, struct cons **);

#endif
