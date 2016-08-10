#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "consts.h"
#include "str.h"
#include "write_vtk.h"

main()
{
 
 double *Q;
 Q = malloc(ndof*sizeof(double));
 
 double *Qrbm;
 Qrbm = malloc(nmodes*sizeof(double));

 int i;
 double **qm;
 qm = malloc(ndof*sizeof(double));
 for (i=0; i<ndof; i++)
 {
  qm[i] = malloc(nmodes*sizeof(double));
 }


 printf("initializing \n");
 init(Q,Qrbm,qm);

 
 printf("advance in time \n");
 advance(Q,Qrbm,qm);

 
 free(Q);
 free(Qrbm);
 free(qm);

}
