#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "str.h"
#include "funcs.h"
#include "consts.h"
#include "flux.h"
#include "chemistry.h"

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

void init(struct fluid **prim, struct cons **Q)
{

 int i, j;
 double p1, c1;
 double ke;

 
 // initialize primitive variables

 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {

   prim[i][j].rho = rho0;
   prim[i][j].u = 0.0;
   prim[i][j].v = 0.0;
   prim[i][j].e = cp[2]/1.4*300.0;
   prim[i][j].YH2 = 0.0;
   prim[i][j].YO2 = 0.0; 
   prim[i][j].YH2O = 1.0; 
 
   compute_pc(prim[i][j].rho,prim[i][j].e,prim[i][j].YH2,prim[i][j].YO2,prim[i][j].YH2O,&p1,&c1);
 
   prim[i][j].p = p1;
   prim[i][j].c = c1;
  
  }
 } 

 
 // initialize conservative variables
 
 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {
   Q[i][j].Q1 = prim[i][j].rho;
   Q[i][j].Q2 = prim[i][j].rho*prim[i][j].u;
   Q[i][j].Q3 = prim[i][j].rho*prim[i][j].v;
   ke = 0.5*(pow(prim[i][j].u,2) + pow(prim[i][j].v,2));
   Q[i][j].Q4 = prim[i][j].rho*(prim[i][j].e + ke);
   Q[i][j].Q5 = prim[i][j].rho*prim[i][j].YH2;
   Q[i][j].Q6 = prim[i][j].rho*prim[i][j].YO2;
   Q[i][j].Q7 = prim[i][j].rho*prim[i][j].YH2O;
  }
 } 



 // write init file
 
 FILE *ff;
 char filename[15];
 printf("writing init file  \n"); 
 sprintf(filename, "snapshots/init"); 
 ff = fopen(filename,"w");

 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  { 
     fprintf(ff,"%24.16f\n",Q[i][j].Q1/(rho0));
     fprintf(ff,"%24.16f\n",Q[i][j].Q2/(rho0*c0));  
     fprintf(ff,"%24.16f\n",Q[i][j].Q3/(rho0*c0));
     fprintf(ff,"%24.16f\n",Q[i][j].Q4/(rho0*e0));  
     fprintf(ff,"%24.16f\n",Q[i][j].Q5/(rho0));
     fprintf(ff,"%24.16f\n",Q[i][j].Q6/(rho0));  
     fprintf(ff,"%24.16f\n",Q[i][j].Q7/(rho0));  
  }
 }

 fclose(ff);

}

//-----------------------------------------------------------------------------------------------


void advance(struct fluid **prim, struct cons **Q)
{
 
 double time = 0.0;
 int i, j, n;

  
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
 
  

  cons_2_prim(prim,Q);  

 int ns = 0;
 int nsnap = write_snap[ns];
 char filename[15];
 FILE *ff;

 sprintf(filename, "snapshots/snapshot_u%d",ns+1); 
 ff = fopen(filename,"a");

 
 for (n=0; n<nmax; n++)
 { 

  time += dt;
  printf("n = %d; time = %f \n",n,time);  

  if (n == nsteps[ns])
  {
   ns++;
   nsnap = write_snap[ns];
   fclose(ff);
   sprintf(filename, "snapshots/snapshot_u%d",ns+1); 
   ff = fopen(filename,"a");
  }
  if(n%nsnap == 0)
  {
   write_snapshot(ff,Q); 
  } 
 


  
  compute_dF(prim,dF);
  compute_dG(prim,dG);
  compute_SS(prim,SS);

  update(Q,dF,dG,SS);

  cons_2_prim(prim,Q);     

 }


 fclose(ff); 

}


//-----------------------------------------------------------------------------------------------

void update(struct cons **Q, struct cons **dF, struct cons **dG, struct cons **SS)
{

 int i, j;

 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {

   Q[i][j].Q1 -= dt*(dF[i][j].Q1/dx + dG[i][j].Q1/dy);
   Q[i][j].Q2 -= dt*(dF[i][j].Q2/dx + dG[i][j].Q2/dy);
   Q[i][j].Q3 -= dt*(dF[i][j].Q3/dx + dG[i][j].Q3/dy);
   Q[i][j].Q4 -= dt*(dF[i][j].Q4/dx + dG[i][j].Q4/dy);
   Q[i][j].Q5 -= dt*(dF[i][j].Q5/dx + dG[i][j].Q5/dy);
   Q[i][j].Q6 -= dt*(dF[i][j].Q6/dx + dG[i][j].Q6/dy);  
   Q[i][j].Q7 -= dt*(dF[i][j].Q7/dx + dG[i][j].Q7/dy);

   Q[i][j].Q1 += dt*SS[i][j].Q1;
   Q[i][j].Q2 += dt*SS[i][j].Q2;   
   Q[i][j].Q3 += dt*SS[i][j].Q3;
   Q[i][j].Q4 += dt*SS[i][j].Q4;   
   Q[i][j].Q5 += dt*SS[i][j].Q5;
   Q[i][j].Q6 += dt*SS[i][j].Q6;   
   Q[i][j].Q7 += dt*SS[i][j].Q7;
   
  }
 }

}

//-----------------------------------------------------------------------------------------------

void write_snapshot(FILE* ff, struct cons **Q)
{

 int i, j;

 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  { 
     fprintf(ff,"%24.16f\n",Q[i][j].Q1/(rho0));
     fprintf(ff,"%24.16f\n",Q[i][j].Q2/(rho0*c0));  
     fprintf(ff,"%24.16f\n",Q[i][j].Q3/(rho0*c0));
     fprintf(ff,"%24.16f\n",Q[i][j].Q4/(rho0*e0));  
     fprintf(ff,"%24.16f\n",Q[i][j].Q5/(rho0));
     fprintf(ff,"%24.16f\n",Q[i][j].Q6/(rho0));  
     fprintf(ff,"%24.16f\n",Q[i][j].Q7/(rho0));  
  }
 }

}

//-----------------------------------------------------------------------------------------------


