#ifndef FLUX_H
#define FLUX_H

void compute_dF(struct fluid**, struct cons**);  
void compute_dG(struct fluid**, struct cons**);  

void compute_flux(int, struct fluid, struct fluid, struct fluid, struct fluid, struct fluid, struct cons*);
void Rusanov_flux(int, struct fluid, struct fluid, struct cons*);
void apply_muscl(struct fluid, struct fluid, struct fluid, struct fluid, struct fluid*, struct fluid*);
void muscl_lr(double, double, double, double, double*, double*);
double compute_phi(double, double, double);

#endif
