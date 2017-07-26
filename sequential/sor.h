#ifndef __SOR_H_
#define __SOR_H_

void sor(
         double omg,
         double dx,
         double dy,
         double dP,
         int imax,
         int jmax,
         double ***P,
         double ***RS,
         double *res,
         int **Flag,
         int k,
         int mode
        );


#endif

