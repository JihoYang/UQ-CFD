#ifndef __SOR_H_
#define __SOR_H_

void sor(
         double omg,
         double dx,
         double dy,
         double dP,
         int il, 
         int ir,
         int jb, 
         int jt,
         int imax,
         int jmax,
         int rank_l, 
         int rank_r, 
         int rank_b, 
         int rank_t,
         double ***P,
         double ***RS,
         double *res,
         int **Flag,
         double *bufSend,
         double *bufRecv,
         int k,
         int mode
        );


#endif

