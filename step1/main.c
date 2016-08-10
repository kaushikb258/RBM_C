#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "consts.h"
#include "str.h"
#include "write_vtk.h"

main()
{
 
 int i;

 struct fluid **prim;
 prim = malloc(imax*sizeof(struct fluid));
 for (i=0; i<imax; i++)
 {
  prim[i] = malloc(jmax*sizeof(struct fluid));
 }
 
 struct cons **Q;
 Q = malloc(imax*sizeof(struct cons));
 for (i=0; i<imax; i++)
 {
  Q[i] = malloc(jmax*sizeof(struct cons));
 }


 printf("initializing \n");
 init(prim,Q);

 
 printf("advance in time \n");
 advance(prim,Q);

 printf("write vtk file \n");
 out_vtk(prim);

}
