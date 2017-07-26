#include "helper.h"
#include "init.h"
#include <string.h>
#include <mpi.h>

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
                    int *iproc,                //Number of decomposed domain in x-direction
                    int *jproc,                //Number of decomposed domain in y-direction
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
   
   READ_INT (szFileName, *iproc);
   READ_INT (szFileName, *jproc);

   READ_INT (szFileName, *mode);

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   return 1;

}


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
            ){

    printf("Flow field initialisation started\n");
    printf("\n");

	init_matrix_3d(U, il - 2, ir + 1, jb - 1, jt + 1, 0, mode, UI);

    printf("U initialised\n");
    printf("\n");

	init_matrix_3d(V, il - 1, ir + 1, jb - 2, jt + 1, 0, mode, VI);

    printf("V initialised\n");
    printf("\n");
    
	init_matrix_3d(P, il - 1, ir + 1, jb - 1, jt + 1, 0, mode, PI);

    printf("P initialised\n");
    printf("\n");

}


void init_flag(
               char *problem,
               int imax,
               int jmax,
               int il, 
               int ir,
               int jb, 
               int jt,
               int **Flag
              ){

    char image_filename[50];
    int i,j;
    int **pic;
    
    strcpy(image_filename, problem);
    strcat(image_filename, ".pgm");

    if(strcmp(problem, "driven_cavity") == 0){
   
        init_imatrix(Flag, il-1, ir+1, jb-1, jt+1, C_F);                                                                                                                         
         
    } 

    else{

        pic = read_pgm(image_filename);

        //Convert pixel values to 1 and 0, if necessary
        if (pic[ir][jt] > 1 || pic[il][jb] > 1){

            for(i = il; i <= ir; i++){
     
                for(j = jb; j <= jt; j++){
     
                    if(pic[i][j] ==  255){
     
                        pic[i][j] = 1;
     
                    }
     
                    else{
     
                        pic[i][j] = 0;
     
                    }
     
                } 
     
            }

        }

        //Pixel values: 1 for fluid, 0 for obstacle
        for(i = max(il - 1, 1); i <= min(ir + 1, imax); i++){

            for(j = max(jb - 1, 1); j <= min(jt + 1, jmax); j++){

                Flag[i][j] = max(pic[i+1][j] * B_O + pic[i-1][j] * B_W + pic[i][j-1] * B_S + pic[i][j+1] * B_N, pic[i][j] * C_F);

            }

        }

    }

}

