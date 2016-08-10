#ifndef GRID_H
#define GRID_H

const int imax;
const int jmax;
const int nmax; 

const double dt;
const double dx;
const double dy; 

const double rho0;
const double p0;
const double e0;
const double c0;
const double t0;

const double cp[3];
const double mw[3];
const double Runiv;

const double TF;
const double uF;
const double pF;
const double TO;
const double uO;
const double pO;

const int ncohort;
const int nsteps[4];
const int write_snap[4];


#endif
