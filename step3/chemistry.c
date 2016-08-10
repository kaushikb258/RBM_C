#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "str.h"
#include "consts.h"
#include "chemistry.h"

void compute_SS(struct fluid **prim, struct cons **SS)
{

 double Arr = 5.5e12;
 double Ea = 8.0e3;
 double Qheat = 8.314*9800.0;
 int i, j;
 double CH2, CO2, RR, mwmix, T; 

 set_dF_zero(SS);



 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {

   SS[i][j].Q1 = 0.0;
   SS[i][j].Q2 = 0.0;
   SS[i][j].Q3 = 0.0;
   SS[i][j].Q4 = 0.0;
   SS[i][j].Q5 = 0.0;
   SS[i][j].Q6 = 0.0;
   SS[i][j].Q7 = 0.0;

   CH2 = prim[i][j].rho*prim[i][j].YH2/mw[0]/1000.0;
   CO2 = prim[i][j].rho*prim[i][j].YO2/mw[1]/1000.0;
   mwmix = 1.0/(prim[i][j].YH2/mw[0] + prim[i][j].YO2/mw[1] + prim[i][j].YH2O/mw[2]);
   T = prim[i][j].p/(Runiv/mwmix)/prim[i][j].rho;

   if(T >= 3000.0)
   {
    printf("T very high in chemistry %f \n",T);
    printf("i, j %d %d \n",i,j);
    printf(" mass fr %f %f %f \n ", prim[i][j].YH2, prim[i][j].YO2, prim[i][j].YH2O);
    printf("rho = %f \n",prim[i][j].rho);
    printf("p = %f \n",prim[i][j].p);
    exit(0);
   } 

   RR = Arr*exp(-Ea/8.314/T)*(CH2*CH2)*(CO2);

   SS[i][j].Q5 = -2.0*RR*mw[0]*1000.0;
   SS[i][j].Q6 = -RR*mw[1]*1000.0;
   SS[i][j].Q7 = 2.0*RR*mw[2]*1000.0;
   SS[i][j].Q4 = Qheat*(2.0*RR)*(1.0e6); 
  }
 }

}  
