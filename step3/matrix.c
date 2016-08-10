#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

void matrix_mult(int n1, int n2, double **first, double *second, double *mult)
{
 
 int i, k;
 double sum = 0.0;

 for (i = 0; i < n1; i++)
 {
   for (k = 0; k < n2; k++)
   {
    sum = sum + first[i][k]*second[k];
   }
   mult[i] = sum;
   sum = 0.0;  
 }
 

}
