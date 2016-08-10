#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "str.h"
#include "funcs.h"
#include "consts.h"
#include "flux.h"
#include "chemistry.h"
#include "matrix.h"

//-----------------------------------------------------------------------------------------------

void compute_p(double rho, double e, double YH2, double YO2, double YH2O, double *p1)
{

 double cpmix, cvmix, mwmix, gammix;
 double T;

 cpmix = YH2*cp[0] + YO2*cp[1] + YH2O*cp[2];
 mwmix = 1.0/(YH2/mw[0] + YO2/mw[1] + YH2O/mw[2]); 
 cvmix = cpmix - Runiv/mwmix;
 gammix = cpmix/cvmix;
 
 T = e/cvmix;

 double pres = rho*(Runiv/mwmix)*T;
 *p1 = pres; 

 if(pres >= 1.0e8)
 {
  printf("pres very high in compute_p (bar) %f \n",pres/1.0e5);
  printf("T = %f \n",T);
  printf("mass fr = %f %f %f \n",YH2,YO2,YH2O); 
  exit(0);
 }  

}


void compute_c_rhoe(double rho, double e, double YH2, double YO2, double YH2O, double *c1)
{

 double cpmix, cvmix, mwmix, gammix;
 double T;

 cpmix = YH2*cp[0] + YO2*cp[1] + YH2O*cp[2];
 mwmix = 1.0/(YH2/mw[0] + YO2/mw[1] + YH2O/mw[2]); 
 cvmix = cpmix - Runiv/mwmix;
 gammix = cpmix/cvmix;
 
 T = e/cvmix;

 double sound = sqrt(gammix*(Runiv/mwmix)*T); 
 *c1 = sound;  

}


void compute_c_rhop(double rho, double p, double YH2, double YO2, double YH2O, double *c1)
{

 double cpmix, cvmix, mwmix, gammix;

 cpmix = YH2*cp[0] + YO2*cp[1] + YH2O*cp[2];
 mwmix = 1.0/(YH2/mw[0] + YO2/mw[1] + YH2O/mw[2]); 
 cvmix = cpmix - Runiv/mwmix;
 gammix = cpmix/cvmix;
 
 double sound = sqrt(gammix*p/rho); 
 *c1 = sound;  

}


//-----------------------------------------------------------------------------------------------

void compute_rho(double p, double T, double YH2, double YO2, double YH2O, double *rho1)
{

 double mwmix = 1.0/(YH2/mw[0] + YO2/mw[1] + YH2O/mw[2]); 
 
 double r = p/(Runiv/mwmix)/T;
 *rho1 = r;  
  
}



void compute_e(double p, double T, double YH2, double YO2, double YH2O, double *e1)
{

 double cpmix, cvmix, mwmix;

 cpmix = YH2*cp[0] + YO2*cp[1] + YH2O*cp[2];
 mwmix = 1.0/(YH2/mw[0] + YO2/mw[1] + YH2O/mw[2]); 
 cvmix = cpmix - Runiv/mwmix;

 double emix = cvmix*T; 
 *e1 = emix;  
  
}


void compute_e_from_prho(double p, double rho, double YH2, double YO2, double YH2O, double *e1)
{

 double cpmix, cvmix, mwmix;

 cpmix = YH2*cp[0] + YO2*cp[1] + YH2O*cp[2];
 mwmix = 1.0/(YH2/mw[0] + YO2/mw[1] + YH2O/mw[2]); 
 cvmix = cpmix - Runiv/mwmix;

 double T = p/rho/(Runiv/mwmix);

 double emix = cvmix*T; 
 *e1 = emix;  
  
}



//-----------------------------------------------------------------------------------------------


void compute_pc(double rho, double e, double YH2, double YO2, double YH2O, double *p1, double *c1)
{
  
 double pres;
 compute_p(rho,e,YH2,YO2,YH2O,&pres);
 *p1 = pres; 

 double sound;
 compute_c_rhoe(rho,e,YH2,YO2,YH2O,&sound); 
 *c1 = sound; 

}


//-----------------------------------------------------------------------------------------------

void cons_2_prim(struct fluid **prim, struct cons **Q)
{
 
 int i, j;
 double ke;
 double p1, c1; 
 double h2, o2, h2o, sumYk;


 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {

   prim[i][j].rho = Q[i][j].Q1;
   prim[i][j].u = Q[i][j].Q2/Q[i][j].Q1;
   prim[i][j].v = Q[i][j].Q3/Q[i][j].Q1; 
   ke = 0.5*(pow(prim[i][j].u,2) + pow(prim[i][j].v,2));
   prim[i][j].e = Q[i][j].Q4/Q[i][j].Q1 - ke;

   if (prim[i][j].e < 0.0)
   {
    printf("e < 0 %f \n",prim[i][j].e);
    exit(0);
   }

   h2 = Q[i][j].Q5/Q[i][j].Q1;
   o2 = Q[i][j].Q6/Q[i][j].Q1;
   h2o = Q[i][j].Q7/Q[i][j].Q1;

   h2 = (h2 > 0.0)? h2 : 0.0;
   h2 = (h2 < 1.0)? h2 : 1.0;
   o2 = (o2 > 0.0)? o2 : 0.0;
   o2 = (o2 < 1.0)? o2 : 1.0;
   h2o = (h2o > 0.0)? h2o : 0.0;
   h2o = (h2o < 1.0)? h2o : 1.0;

   sumYk = h2 + o2 + h2o;
   prim[i][j].YH2 = h2/sumYk;
   prim[i][j].YO2 = o2/sumYk;
   prim[i][j].YH2O = h2o/sumYk;

   compute_pc(prim[i][j].rho,prim[i][j].e,prim[i][j].YH2,prim[i][j].YO2,prim[i][j].YH2O,&p1,&c1);

   prim[i][j].p = p1;
   prim[i][j].c = c1;

   if(prim[i][j].p <= 0.0)
   {
    printf("p<0 %f \n",prim[i][j].p);
    exit(0);
   }
 

  }
 } 



}


//-----------------------------------------------------------------------------------------------

void init(double *Q, double *Qrbm, double **qm)
{

 int i;
 

 FILE *ff;
 ff = fopen("../step2/octave/qtil0","r");
 for (i = 0; i<nmodes; i++)
 {
  fscanf(ff,"%lf \n",&Qrbm[i]);
 }  
 fclose(ff);


 // read qm
 read_pod(qm,1);  

 // compute Q
 matrix_mult(ndof,nmodes,qm,Qrbm,Q);

}

//-----------------------------------------------------------------------------------------------


void advance(double *Q, double *Qrbm, double **qm)
{
 
 double time = 0.0;
 int i, j, k, n;


 double *rhs;
 rhs = malloc(ndof*sizeof(double));
 for (i=0; i<ndof; i++) rhs[i] = 0.0; 

 double *dQrbm;
 dQrbm = malloc(nmodes*sizeof(double));

  
 struct cons **dF;
 dF = malloc(imax*sizeof(struct cons));
 for (i=0; i<imax; i++)
 {
  dF[i] = malloc(jmax*sizeof(struct cons));
 }
 

 struct cons **dG;
 dG = malloc(imax*sizeof(struct cons));
 for (i=0; i<imax; i++)
 {
  dG[i] = malloc(jmax*sizeof(struct cons));
 }
 

 struct cons **SS;
 SS = malloc(imax*sizeof(struct cons));
 for (i=0; i<imax; i++)
 {
  SS[i] = malloc(jmax*sizeof(struct cons));
 }
 

 struct fluid **prim;
 prim = malloc(imax*sizeof(struct fluid));
 for (i=0; i<imax; i++)
 {
  prim[i] = malloc(jmax*sizeof(struct fluid));
 }
 
 
 double **qmT;
 qmT = malloc(nmodes*sizeof(double));
 for (i=0; i<nmodes; i++)
 {
  qmT[i] = malloc(ndof*sizeof(double));
 }

 
 struct cons **Qc;
 Qc = malloc(imax*sizeof(struct cons));
 for (i=0; i<imax; i++)
 {
  Qc[i] = malloc(jmax*sizeof(struct cons));
 }



 // compute qmT
 compute_qmT(qm,qmT);
 

 // initialize prim
 from_Q_compute_prim(Q,Qc,prim);


  
 int ns = 0;

 
 for (n=0; n<nmax; n++)
 { 

  time += dt;
  if(n%10 == 0) printf("n = %d; time = %f \n",n,time);  

  if (n == nsteps[ns])
  {
   ns++;   
   set_2_zero(ndof,nmodes,qm);   
   read_pod(qm,ns+1);
   compute_qmT(qm,qmT);

   // recompute Qrbm
   for (i=0; i<nmodes; i++) Qrbm[i] = 0.0;
   matrix_mult(nmodes,ndof,qmT,Q,Qrbm);
  }
  
 
  compute_dF(prim,dF);
  compute_dG(prim,dG);
  compute_SS(prim,SS);
  
  compute_rhs(dF,dG,SS,rhs);

  for (i=0; i<nmodes; i++) dQrbm[i] = 0.0;    
  matrix_mult(nmodes,ndof,qmT,rhs,dQrbm);

  for (i=0; i<nmodes; i++) Qrbm[i] += dQrbm[i];


  // compute Q
  for (i=0; i<ndof; i++) Q[i] = 0.0;
  matrix_mult(ndof,nmodes,qm,Qrbm,Q);

  
  from_Q_compute_prim(Q,Qc,prim); 

//--------------------------------


 }
 

  printf("write vtk file \n");
  out_vtk(prim);


  free(rhs);
  free(dQrbm);
  free(dF);
  free(dG);
  free(SS);
  free(qmT);
  free(prim);  
  free(Qc);

}



//-----------------------------------------------------------------------------------------------


void read_pod(double **qm, int ns)
{
  
 char filename[25];
 FILE *ff;
 sprintf(filename, "../step2/octave/pod_modes%d",ns); 
 
 printf("reading POD file: %s \n",filename);

 ff = fopen(filename,"r");
 int i, j; 
  
 for (i=0; i<ndof; i++)
 {
  for (j=0; j<nmodes; j++)
  {
   fscanf(ff,"%lf",&qm[i][j]);
  }
  fscanf(ff,"\n");
 }

 fclose(ff); 

}

//-----------------------------------------------------------------------------------------------

void set_2_zero(int n1, int n2, double **A)
{
 int i, j;
 
 for (i=0; i<n1; i++)
 {
  for (j=0; j<n2; j++)
  {
   A[i][j] = 0.0;
  }
 }
 
}

//-----------------------------------------------------------------------------------------------

void from_Q_compute_prim(double *Q, struct cons **Qc, struct fluid **prim)
{

 int i, j, k;



 k = -1; 
 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {
   k++;
   Qc[i][j].Q1 = Q[k]*(rho0);
   k++;
   Qc[i][j].Q2 = Q[k]*(rho0*c0);
   k++;
   Qc[i][j].Q3 = Q[k]*(rho0*c0);
   k++;
   Qc[i][j].Q4 = Q[k]*(rho0*e0);
   k++;
   Qc[i][j].Q5 = Q[k]*(rho0);
   k++;
   Qc[i][j].Q6 = Q[k]*(rho0);
   k++;
   Qc[i][j].Q7 = Q[k]*(rho0);
  }
 } 

 cons_2_prim(prim,Qc);


}

//-----------------------------------------------------------------------------------------------

void compute_rhs(struct cons **dF, struct cons **dG, struct cons **SS, double *rhs)
{

 int i, j, k;

 k = -1;

 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {
   k++;
   rhs[k] = dt*(-dF[i][j].Q1/dx - dG[i][j].Q1/dy + SS[i][j].Q1)/(rho0);
   k++;
   rhs[k] = dt*(-dF[i][j].Q2/dx - dG[i][j].Q2/dy + SS[i][j].Q2)/(rho0*c0);
   k++;
   rhs[k] = dt*(-dF[i][j].Q3/dx - dG[i][j].Q3/dy + SS[i][j].Q3)/(rho0*c0);
   k++;
   rhs[k] = dt*(-dF[i][j].Q4/dx - dG[i][j].Q4/dy + SS[i][j].Q4)/(rho0*e0);
   k++;
   rhs[k] = dt*(-dF[i][j].Q5/dx - dG[i][j].Q5/dy + SS[i][j].Q5)/(rho0);
   k++;
   rhs[k] = dt*(-dF[i][j].Q6/dx - dG[i][j].Q6/dy + SS[i][j].Q6)/(rho0);
   k++;
   rhs[k] = dt*(-dF[i][j].Q7/dx - dG[i][j].Q7/dy + SS[i][j].Q7)/(rho0);
  }
 }  

}

//-----------------------------------------------------------------------------------------------

void compute_qmT(double **qm, double **qmT)
{

 int i, j;

 for (i=0; i<ndof; i++)
 {
  for (j=0; j<nmodes; j++)
  {
   qmT[j][i] = qm[i][j];
  }
 }

}

//-----------------------------------------------------------------------------------------------

