#ifndef __UVP_H__
#define __UVP_H__

void calculate_fg(
                  double Re,
                  double GX,
                  double GY,
                  double alpha,
                  double dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double ***U,
                  double ***V,
                  double ***F,
                  double ***G,
                  int **Flag,
                  double lambda,
                  int k,
                  int mode,
                  double **F_sg1,
                  double **F_sg2,
                  double **G_sg1,
                  double **G_sg2
                );


void calculate_rs(
                  double dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double ***F,
                  double ***G,
                  double ***RS,      
                  int k,
                  int mode
                );


void calculate_dt(
                  double Re,
                  double tau,
                  double *dt,
                  double dx,
                  double dy,
                  double ***U,
                  double ***V,
                  int imax,
                  int jmax
                );


void calculate_uv(
                  double dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double ***U,
                  double ***V,
                  double ***F,
                  double ***G,
                  double ***P,
                  int **Flag,
                  int k,
                  int mode
                );

#endif

