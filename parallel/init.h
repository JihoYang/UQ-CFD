#ifndef __INIT_H_
#define __INIT_H_

#define CUT_CELL 1

#define C_F         16
#define C_B         0
#define B_O         8
#define B_W         4
#define B_S         2
#define B_N         1  
#define B_SO        10  
#define B_SW        6   
#define B_NO        9
#define B_NW        5

#define NO_SLIP     1
#define FREE_SLIP   2
#define OUTFLOW     3

int read_parameters( 
                    const char *szFileName,
                    double *Re,
                    double *UI,
                    double *VI,
                    double *PI,
                    double *GX,
                    double *GY,
                    double *t_end,
                    double *xlength,
                    double *ylength,
                    double *dt,
                    double *dx,
                    double *dy,
                    int  *imax,
                    int  *jmax,
                    double *alpha,
                    double *omg,
                    double *tau,
                    int  *itermax,
                    double *eps,
                    double *dt_value,
                    double *dp,
                    int *wl,
                    int *wr,
                    int *wt,
                    int *wb,
                    int *iproc,
                    int *jproc,
                    int *mode
                   );


void init_uvp(
              double UI,
              double VI,
              double PI,
              int il, 
              int ir,
              int jb, 
              int jt,
              double ***U,
              double ***V,
              double ***P,
              int mode
             );


void init_flag(
               char *problem,
               int imax,
               int jmax,
               int il, 
               int ir,
               int jb, 
               int jt,
               int **Flag
              );

#endif

