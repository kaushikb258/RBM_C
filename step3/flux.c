#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "str.h"
#include "flux.h"
#include "consts.h"


//-------------------------------------------------------------

void set_dF_zero(struct cons **dF)
{
 int i, j;

 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {
   dF[i][j].Q1 = 0.0;
   dF[i][j].Q2 = 0.0;  
   dF[i][j].Q3 = 0.0;
   dF[i][j].Q4 = 0.0;  
   dF[i][j].Q5 = 0.0;
   dF[i][j].Q6 = 0.0;  
   dF[i][j].Q7 = 0.0;  
  }
 }  
 
}

//-------------------------------------------------------------

void compute_dF(struct fluid **prim, struct cons **dF)
{
 int i, j;
 struct fluid cell_im2, cell_im1, cell_i, cell_ip1, cell_ip2; 
 struct fluid cell_inj; 
 double ui, vi, Ti, pi, h2i, o2i, h2oi; 
 double rr, cc, ee;
 int jf;

 set_dF_zero(dF); 

 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {
 
   jf = jmax - j;

//-----injector------
  if(i <= 1)
  {

   if(jf >= 20 && jf <= 30)
   {
    // fuel
    ui = uF;
    vi = 0.0;
    pi = pF;
    Ti = TF;
    h2i = 1.0;
    o2i = 0.0; 
   }  

   if(jf >= 31 && jf <= 69)
   {
    // oxidizer
    ui = uO;
    vi = 0.0;
    pi = pO;
    Ti = TO;
    h2i = 0.0;
    o2i = 1.0; 
   }  

   if(jf >= 70 && jf <= 80)
   {
    // fuel
    ui = uF;
    vi = 0.0;
    pi = pF;
    Ti = TF;
    h2i = 1.0;
    o2i = 0.0; 
   }      

   h2oi = 0.0;
 
   compute_rho(pi,Ti,h2i,o2i,h2oi,&rr);
   compute_c_rhop(rr,pi,h2i,o2i,h2oi,&cc);
   compute_e(pi,Ti,h2i,o2i,h2oi,&ee);

   cell_inj.rho = rr; 
   cell_inj.u = ui;
   cell_inj.v = vi;
   cell_inj.e = ee;
   cell_inj.p = pi;
   cell_inj.c = cc;
   cell_inj.YH2 = h2i;
   cell_inj.YO2 = o2i;  
   cell_inj.YH2O = h2oi;

  }
//-----injector------   
   
   cell_i = prim[i][j];
 
   if(i > 0)
   {
    cell_im1 = prim[i-1][j];
   }
   if(i > 1)
   {
    cell_im2 = prim[i-2][j];
   }
   if(i < imax-1)  
   {
    cell_ip1 = prim[i+1][j];
   }
   if(i < imax-2)
   {
    cell_ip2 = prim[i+2][j];
   }

   
  if(i ==0 || i ==1)
  {
   if(jf < 20 || jf > 80)
   {
    // wall
    if(i == 0)
    {
     cell_im1 = prim[i][j];
     cell_im1.u = -prim[i][j].u;

     cell_im2 = prim[i+1][j];
     cell_im2.u = -prim[i+1][j].u;
    }
    else if(i == 1)
    {
     cell_im2 = prim[i-1][j];
     cell_im2.u = -prim[i-1][j].u;
    }
   }
   else
   {
    // injector
    if(i==0 || i ==1)
    {
     cell_im2 = cell_inj;
    }
    if(i == 0)
    {
     cell_im1 = cell_inj;
    }
   }
  }


   else if(i == imax-1 || i == imax-2)
   {
    // supersonic outflow BC
    cell_ip2 = prim[imax-1][j];
    if(i == imax-1)
    {
     cell_ip1 = prim[imax-1][j];
    }
   }   
   

   // compute flux  
   compute_flux(1,cell_im2,cell_im1,cell_i,cell_ip1,cell_ip2,&dF[i][j]); 

  }
 }


}  


//-------------------------------------------------------------

void compute_dG(struct fluid **prim, struct cons **dG)
{
 int i, j;
 struct fluid cell_jm2, cell_jm1, cell_j, cell_jp1, cell_jp2; 


 set_dF_zero(dG); 

 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {

   cell_j = prim[i][j];

   if(j > 0)
   {
    cell_jm1 = prim[i][j-1];
   }
   if(j > 1)
   {
    cell_jm2 = prim[i][j-2];
   }
   if (j < jmax-1)
   {
    cell_jp1 = prim[i][j+1];
   }
   if (j < jmax-2)
   {
    cell_jp2 = prim[i][j+2];
   }

  
   if (j == 0)
   {
    // wall BC
    cell_jm2 = prim[i][j+1];
    cell_jm2.v = -prim[i][j+1].v; 
    cell_jm1 = prim[i][j];
    cell_jm1.v = -prim[i][j].v; 
   }
   else if(j == 1)
   {
    // wall BC
    cell_jm2 = prim[i][j-1];
    cell_jm2.v = -prim[i][j-1].v;
   } 
   else if (j == jmax-1) 
   {
    // wall BC
    cell_jp1 = prim[i][j];
    cell_jp1.v = -prim[i][j].v;    
    cell_jp2 = prim[i][j-1];
    cell_jp2.v = -prim[i][j-1].v;       
   }
   else if(j == jmax-2)
   {    
    // wall BC
    cell_jp2 = prim[i][j+1];
    cell_jp2.v = -prim[i][j+1].v;
   } 
   
   
   // compute flux  
   compute_flux(2,cell_jp2,cell_jp1,cell_j,cell_jm1,cell_jm2,&dG[i][j]); 

  }
 }


}  


//-------------------------------------------------------------


void compute_flux(int index, struct fluid cell_im2, struct fluid cell_im1, struct fluid cell_i, struct fluid cell_ip1, struct fluid cell_ip2, struct cons *d_F)
{
 struct fluid left, right;
 struct cons flux_l, flux_r; 
 //struct cons diff; 

 // i+1/2
 apply_muscl(cell_im1,cell_i,cell_ip1,cell_ip2,&left,&right);
 Rusanov_flux(index,left,right,&flux_r);

 // i-1/2
 apply_muscl(cell_im2,cell_im1,cell_i,cell_ip1,&left,&right);
 Rusanov_flux(index,left,right,&flux_l);

 // flux difference
 d_F->Q1 = flux_r.Q1 - flux_l.Q1;
 d_F->Q2 = flux_r.Q2 - flux_l.Q2;
 d_F->Q3 = flux_r.Q3 - flux_l.Q3;
 d_F->Q4 = flux_r.Q4 - flux_l.Q4;
 d_F->Q5 = flux_r.Q5 - flux_l.Q5;
 d_F->Q6 = flux_r.Q6 - flux_l.Q6;
 d_F->Q7 = flux_r.Q7 - flux_l.Q7;

 //d_F = &diff;

}


//-------------------------------------------------------------

void Rusanov_flux(int index, struct fluid left, struct fluid right, struct cons *flux)
{
  
 double vel1, vel2, cmax, nx, ny, ke;
 double fluxl, fluxr, ql, qr, consl, consr, energyl, energyr;

 if (index == 1)
 {
  nx = 1.0;
  ny = 0.0;
  ql = left.u;
  qr = right.u;
 }
 else if (index == 2)
 {
  nx = 0.0;
  ny = 1.0;
  ql = left.v;
  qr = right.v;
 }
 else
 {
  printf("wrong index %d \n",index);
  exit(0);
 } 


 vel1 = abs(ql) + left.c;
 vel2 = abs(qr) + right.c;
 cmax = (vel1 >= vel2)? vel1 : vel2;
 

 // continuity
 consl = left.rho;
 consr = right.rho;
 fluxl = left.rho*ql;
 fluxr = right.rho*qr;
 flux->Q1 = 0.5*(fluxl+fluxr) - 0.5*cmax*(consr-consl);

 // x-momentum
 consl = left.rho*left.u;
 consr = right.rho*right.u;
 fluxl = left.rho*left.u*ql + nx*left.p;
 fluxr = right.rho*right.u*qr + nx*right.p;
 flux->Q2 = 0.5*(fluxl+fluxr) - 0.5*cmax*(consr-consl);

 // y-momentum
 consl = left.rho*left.v;
 consr = right.rho*right.v;
 fluxl = left.rho*left.v*ql + ny*left.p;
 fluxr = right.rho*right.v*qr + ny*right.p;
 flux->Q3 = 0.5*(fluxl+fluxr) - 0.5*cmax*(consr-consl);

 // energy
 ke = 0.5*(pow(left.u,2) + pow(left.v,2));
 energyl = left.e + ke;
 ke = 0.5*(pow(right.u,2) + pow(right.v,2));
 energyr = right.e + ke;

 consl = left.rho*energyl;
 consr = right.rho*energyr;
 fluxl = left.rho*ql*energyl + ql*left.p;
 fluxr = right.rho*qr*energyr + qr*right.p;
 flux->Q4 = 0.5*(fluxl+fluxr) - 0.5*cmax*(consr-consl);

 // H2
 consl = left.rho*left.YH2;
 consr = right.rho*right.YH2;
 fluxl = left.rho*left.YH2*ql;
 fluxr = right.rho*right.YH2*qr;
 flux->Q5 = 0.5*(fluxl+fluxr) - 0.5*cmax*(consr-consl);

 // O2
 consl = left.rho*left.YO2;
 consr = right.rho*right.YO2;
 fluxl = left.rho*left.YO2*ql;
 fluxr = right.rho*right.YO2*qr;
 flux->Q6 = 0.5*(fluxl+fluxr) - 0.5*cmax*(consr-consl);

 // H2O
 consl = left.rho*left.YH2O;
 consr = right.rho*right.YH2O;
 fluxl = left.rho*left.YH2O*ql;
 fluxr = right.rho*right.YH2O*qr;
 flux->Q7 = 0.5*(fluxl+fluxr) - 0.5*cmax*(consr-consl);

}

//-------------------------------------------------------------

void apply_muscl(struct fluid uim1, struct fluid ui, struct fluid uip1, struct fluid uip2, struct fluid *uiphL, struct fluid *uiphR)
{
 
 double left, right;
 double sumyk;

 // rho 
 muscl_lr(uim1.rho,ui.rho,uip1.rho,uip2.rho,&left,&right);
 uiphL->rho = left;
 uiphR->rho = right; 
 
 // u 
 muscl_lr(uim1.u,ui.u,uip1.u,uip2.u,&left,&right);
 uiphL->u = left;
 uiphR->u = right; 
 
 // v 
 muscl_lr(uim1.v,ui.v,uip1.v,uip2.v,&left,&right);
 uiphL->v = left;
 uiphR->v = right; 

 // p 
 muscl_lr(uim1.p,ui.p,uip1.p,uip2.p,&left,&right);
 uiphL->p = left;
 uiphR->p = right; 

 // YH2 
 muscl_lr(uim1.YH2,ui.YH2,uip1.YH2,uip2.YH2,&left,&right);
 uiphL->YH2 = left;
 uiphR->YH2 = right; 

 // YO2 
 muscl_lr(uim1.YO2,ui.YO2,uip1.YO2,uip2.YO2,&left,&right);
 uiphL->YO2 = left;
 uiphR->YO2 = right; 

 // YH2O 
 muscl_lr(uim1.YH2O,ui.YH2O,uip1.YH2O,uip2.YH2O,&left,&right);
 uiphL->YH2O = left;
 uiphR->YH2O = right; 


 uiphL->YH2 = (uiphL->YH2 > 0.0)? uiphL->YH2 : 0.0;
 uiphL->YH2 = (uiphL->YH2 < 1.0)? uiphL->YH2 : 1.0;
 uiphL->YO2 = (uiphL->YO2 > 0.0)? uiphL->YO2 : 0.0;
 uiphL->YO2 = (uiphL->YO2 < 1.0)? uiphL->YO2 : 1.0;
 uiphL->YH2O = (uiphL->YH2O > 0.0)? uiphL->YH2O : 0.0;
 uiphL->YH2O = (uiphL->YH2O < 1.0)? uiphL->YH2O : 1.0;

 sumyk = uiphL->YH2 + uiphL->YO2 + uiphL->YH2O;
 uiphL->YH2 = uiphL->YH2/sumyk;
 uiphL->YO2 = uiphL->YO2/sumyk;
 uiphL->YH2O = uiphL->YH2O/sumyk;

 sumyk = uiphR->YH2 + uiphR->YO2 + uiphR->YH2O;
 uiphR->YH2 = uiphR->YH2/sumyk;
 uiphR->YO2 = uiphR->YO2/sumyk;
 uiphR->YH2O = uiphR->YH2O/sumyk;

 double e1;
 compute_e_from_prho(uiphL->p,uiphL->rho,uiphL->YH2,uiphL->YO2,uiphL->YH2O,&e1);
 uiphL->e = e1;

 double c1;
 compute_c_rhoe(uiphL->rho,uiphL->e,uiphL->YH2,uiphL->YO2,uiphL->YH2O,&c1);
 uiphL->c = c1;

 compute_e_from_prho(uiphR->p,uiphR->rho,uiphR->YH2,uiphR->YO2,uiphR->YH2O,&e1);
 uiphR->e = e1;

 compute_c_rhoe(uiphR->rho,uiphR->e,uiphR->YH2,uiphR->YO2,uiphR->YH2O,&c1);
 uiphR->c = c1;

}

//-------------------------------------------------------------

void muscl_lr(double uim1, double ui, double uip1, double uip2, double *uiphL, double *uiphR)
{

 double phi;
 const double kappa = 1.0/3.0;
 double duimh, duiph, duip3h, t1, t2;

 duimh = ui - uim1;
 duiph = uip1 - ui;
 duip3h = uip2 - uip1;

 phi = compute_phi(uim1,ui,uip1);
 t1 = (1.0 - kappa)*duimh;
 t2 = (1.0 + kappa)*duiph;
 *uiphL = ui + phi/4.0*(t1 + t2);

 phi = compute_phi(ui,uip1,uip2);
 t1 = (1.0 - kappa)*duip3h;
 t2 = (1.0 + kappa)*duiph;
 *uiphR = uip1 - phi/4.0*(t1 + t2);   

}

//-------------------------------------------------------------


double compute_phi(double uim1, double ui, double uip1)
{
 double ri, phi;
 
 if(ui == uip1)
 {
  phi = 0.0;
 }  
 else
 {
  ri = (ui-uim1)/(uip1-ui);
  phi = 2.0*ri/(1.0 + ri*ri);
 }

 return phi;
}

//-------------------------------------------------------------

