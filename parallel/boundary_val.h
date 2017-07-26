#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

void boundaryvalues(
                    int il, 
                    int ir,
                    int jb, 
                    int jt,
                    int imax,
                    int jmax,
                    double ***U,
                    double ***V,
                    int wl,
                    int wr,
                    int wt,
                    int wb,
                    int **Flag,
                    int k,
                    int mode
                   );

void spec_boundary_val(
                       char *problem, 
                       int il, 
                       int ir, 
                       int jb, 
                       int jt, 
                       int imax, 
                       int jmax, 
                       double ***U, 
                       double ***V, 
                       double ***P, 
                       double Re, 
                       double xlength, 
                       double ylength, 
                       double dP,
                       int k,
                       int mode
                      );

#endif

