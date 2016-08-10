#include <math.h>
#include "consts.h"


// GRID

const int imax = 300;
const int jmax = 100;
const int nmax = 30000;

const double dt = 5.0e-8;
const double dx = 1.0e-3;
const double dy = 1.0e-3; 



// REF VALUES

const double rho0 = 60.0;
const double p0 = 50.0e5;
const double e0 = 2.5e6;
const double c0 = 341.56; //sqrt(1.4*p0/rho0);
const double t0 = 2.928e-6; //dx/c0;



// THERMOD

const double cp[3] = {14.4e3, 918.4, 1956.7};
const double mw[3] = {2.0, 32.0, 18.0};
const double Runiv = 8314.0;


// INFLOW
const double TF = 150.0;
const double uF = 50.0;
const double pF = 100.0e5;
const double TO = 100.0;
const double uO = 20.0;
const double pO = 100.0e5;


// SNAPSHOTS
const int ncohort = 4;
const int nsteps[4] = {2000, 6000, 14000, 30000};


// RBM
const int nmodes = 60;
const int ndof = 300*100*7;  //imax*jmax*7;


