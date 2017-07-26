#include "helper.h"
#include "init.h"
#include <string.h>

int read_parameters( 
                    const char *szFileName,    //Name of the file
                    double *Re,                //Reynolds number   
                    double *UI,                //Initial velocity x-direction
                    double *VI,                //Initial velocity y-direction
                    double *PI,                //Pressure
                    double *GX,                //Gravitation x-direction
                    double *GY,                //Gravitation y-direction
                    double *t_end,             //End time
                    double *xlength,           //Length of the domain x-dir
                    double *ylength,           //Length of the domain y-dir
                    double *dt,                //Time step
                    double *dx,                //Length of a cell x-direction
                    double *dy,                //Length of a cell y-direction
                    int  *imax,                //Number of cells x-direction
                    int  *jmax,                //Number of cells y-direction
                    double *alpha,             //Upwind differencing factor
                    double *omg,               //Relaxation factor 
                    double *tau,               //Safety factor for time step
                    int  *itermax,             //Maximum number of SOR iteration
                    double *eps,               //Accuracy bound for pressure
		            double *dt_value,          //Time for output 
		            double *dP,                //Pressure difference
	                int *wl,                   //Left wall
                    int *wr,                   //Right wall
                    int *wt,                   //Top wall
                    int *wb,                   //Bottom wall
                    int *mode                  //Number of modes for Stochastic Galerkin (P)
                   ){

   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );

   READ_INT   ( szFileName, *itermax );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *dP );
   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );

   READ_INT (szFileName, *wl);
   READ_INT (szFileName, *wr);
   READ_INT (szFileName, *wt);
   READ_INT (szFileName, *wb);

   READ_INT (szFileName, *mode);

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   return 1;

}


void init_uvp(
              double UI,
              double VI,
              double PI,
              double ***U,
              double ***V,
              double ***P,
              int imax,
              int jmax,
              int mode
            ){

	init_matrix_3d(U, 0, imax,     0, jmax + 1, 0, mode, UI);
	init_matrix_3d(V, 0, imax + 1, 0, jmax,     0, mode, VI);
	init_matrix_3d(P, 0, imax + 1, 0, jmax + 1, 0, mode, PI);

}


void init_flag(
               char *problem,
               int imax,
               int jmax,
               int **Flag
              ){

    char image_filename[50];
    int i,j;
    int **pic;
    
    strcpy(image_filename, problem);
    strcat(image_filename, ".pgm");

    if(strcmp(problem, "driven_cavity") == 0){
   
        init_imatrix(Flag, 1, imax, 1, jmax, C_F);
         
    } 

    else{

        pic = read_pgm(image_filename);

        //Pixel values: 1 for fluid, 0 for obstacle
        for(i = 1; i <= imax; i++){

            for(j = 1; j <= jmax; j++){

                Flag[i][j] = max(pic[i+1][j] * B_O + pic[i-1][j] * B_W + pic[i][j-1] * B_S + pic[i][j+1] * B_N, pic[i][j] * C_F);

            }

        }

    }

}

