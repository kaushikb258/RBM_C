#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "write_vtk.h"
#include "consts.h"
#include "str.h"



void out_vtk(struct fluid **fl)
{

 FILE *f;

 int kmax = 1;

 f = fopen("output.vtk", "w+");

 char c1[100] = "# vtk DataFile Version 2.0";
 fprintf(f,"%s\n",c1);
 char title[100] = "CFD H2-O2 flame";
 fprintf(f,"%s\n",title);
 char c2[100] = "ASCII";
 fprintf(f,"%s\n",c2);
 char c3[100] = "DATASET STRUCTURED_GRID";
 fprintf(f,"%s\n",c3);
 char c4[100] = "DIMENSIONS";
 char s_imax[5];
 snprintf(s_imax,5,"%d",imax);
 char s_jmax[5];
 snprintf(s_jmax,5,"%d",jmax);
 char s_kmax[5];
 snprintf(s_kmax,5,"%d",kmax);
 fprintf(f,"%s \t %s \t\t %s \t\t %s \n",c4,s_imax,s_jmax,s_kmax);
 char c5[100] = "POINTS";
 char s_node_num[10];
 snprintf(s_node_num,10,"%d",imax*jmax);
 char c6[100] = "double";
 fprintf(f,"%s \t %s \t %s \n",c5,s_node_num,c6);
 
 double coord[imax*jmax*3];
 int i, j, n;

 n = 0;
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  coord[n] = ((double) (i) + 0.5)*dx; 
  coord[n+1] = (double)(jmax)*dy - ((double) (j) + 0.5)*dy; 
  coord[n+2] = 0.0; //z = 0 as 2D grid
  n += 3;
 }
}

 n = 0;
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  fprintf(f,"%f %s %f %s %f \n ",coord[n]," ",coord[n+1]," ",coord[n+2]); 
  n += 3;
 }
 }
 
 char s_cells[10]; 
 snprintf(s_cells,10,"%d",(imax-1)*(jmax-1));
 char c7[100] = "CELL_DATA ";
 fprintf(f,"%s \t %s \n",c7,s_cells); 
 char c8[100] = "POINT_DATA ";
 fprintf(f,"%s \t %s \n",c8,s_node_num);
 
//-------------

 char c9[100] = "SCALARS density double "; 
 fprintf(f,"%s \n",c9); 
 char c10[100] = "LOOKUP_TABLE default ";
 fprintf(f,"%s \n",c10); 
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  fprintf(f,"%f %s",fl[i][j].rho," ");
 } 
 }
 fprintf(f,"\n");

//-------------

 char c11[100] = "SCALARS u double "; 
 fprintf(f,"%s \n",c11); 
 char c12[100] = "LOOKUP_TABLE default ";
 fprintf(f,"%s \n",c12); 
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  fprintf(f,"%f %s",fl[i][j].u," ");
 } 
 }
 fprintf(f,"\n");

//-------------

 char c13[100] = "SCALARS v double "; 
 fprintf(f,"%s \n",c13); 
 char c14[100] = "LOOKUP_TABLE default ";
 fprintf(f,"%s \n",c14); 
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  fprintf(f,"%f %s",fl[i][j].v," ");
 } 
 }
 fprintf(f,"\n");

//-------------

 char c15[100] = "SCALARS pressure double "; 
 fprintf(f,"%s \n",c15); 
 char c16[100] = "LOOKUP_TABLE default ";
 fprintf(f,"%s \n",c16); 
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  fprintf(f,"%f %s",fl[i][j].p," ");
 } 
 }
 fprintf(f,"\n");

//-------------
 
 char c17[100] = "SCALARS YH2 double "; 
 fprintf(f,"%s \n",c17); 
 char c18[100] = "LOOKUP_TABLE default ";
 fprintf(f,"%s \n",c18); 
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  fprintf(f,"%f %s",fl[i][j].YH2," ");
 } 
 }
 fprintf(f,"\n");

//-------------

 char c19[100] = "SCALARS YO2 double "; 
 fprintf(f,"%s \n",c19); 
 char c20[100] = "LOOKUP_TABLE default ";
 fprintf(f,"%s \n",c20); 
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  fprintf(f,"%f %s",fl[i][j].YO2," ");
 } 
 }
 fprintf(f,"\n");

//-------------

 char c21[100] = "SCALARS YH2O double "; 
 fprintf(f,"%s \n",c21); 
 char c22[100] = "LOOKUP_TABLE default ";
 fprintf(f,"%s \n",c22); 
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  fprintf(f,"%f %s",fl[i][j].YH2O," ");
 } 
 }
 fprintf(f,"\n");

//-------------

 double T[imax][jmax];
 double mwmix, cpmix, cvmix;

 for (i=0; i<imax; i++)
 {
  for (j=0; j<jmax; j++)
  {
   mwmix = 1.0/(fl[i][j].YH2/mw[0] + fl[i][j].YO2/mw[1] + fl[i][j].YH2O/mw[2]);
   cpmix = fl[i][j].YH2*cp[0] + fl[i][j].YO2*cp[1] + fl[i][j].YH2O*cp[2];
   cvmix = cpmix - Runiv/mwmix;
   T[i][j] = fl[i][j].e/cvmix;  
  }
 }
 

 char c23[100] = "SCALARS Temperature double "; 
 fprintf(f,"%s \n",c23); 
 char c24[100] = "LOOKUP_TABLE default ";
 fprintf(f,"%s \n",c24); 
 for (j=0; j<jmax; j++)
 {
 for (i=0; i<imax; i++)
 {
  fprintf(f,"%f %s",T[i][j]," ");
 } 
 }
 fprintf(f,"\n");
 

//-------------


 fclose(f);

}
