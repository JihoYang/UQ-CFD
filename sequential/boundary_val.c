#include "init.h"
#include "helper.h"
#include <string.h>
#include <stdio.h>

void boundaryvalues(
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
                   ){

	int i, j;
	
	//Left boundary
    switch(wl){

        case NO_SLIP:

            for(j = 0; j <= jmax; j++){

                U[k][0][j] = 0;
                V[k][0][j] = -V[k][1][j];

            }

            break;

        case FREE_SLIP:

            for(j = 0; j <= jmax; j++){

                U[k][0][j] = 0;
                V[k][0][j] = V[k][1][j];

            }

            break;

        case OUTFLOW:

            for(j = 0; j <= jmax; j++){

                U[k][0][j] = U[k][1][j];
                V[k][0][j] = V[k][1][j];

            }

            break;

    }

	//Right boundary 
    switch(wr){

        case NO_SLIP:

            for(j = 0; j <= jmax; j++){

                U[k][imax][j] = 0;
                V[k][imax+1][j] = -V[k][imax][j];

            }

            break;

        case FREE_SLIP:

            for(j = 0; j <= jmax; j++){

                U[k][imax][j] = 0;
                V[k][imax+1][j] = V[k][imax][j];

            }

            break;

        case OUTFLOW:

            for(j = 0; j <= jmax; j++){

                U[k][imax][j] = U[k][imax-1][j];
                V[k][imax+1][j] = V[k][imax][j];

            }

            break;

    }
	
	//Top boundary
    switch(wt){

        case NO_SLIP:

            for(i = 0; i <= imax; i++){

                U[k][i][jmax+1] = -U[k][i][jmax];
                V[k][i][jmax] = 0;

            }

            break;

        case FREE_SLIP:

            for(i = 0; i <= imax; i++){

                U[k][i][jmax+1] = U[k][i][jmax];
                V[k][i][jmax] = 0;

            }

            break;

        case OUTFLOW:

            for(i = 0; i <= imax; i++){

                U[k][i][jmax+1] = U[k][i][jmax];
                V[k][i][jmax] = V[k][i][jmax-1];

            }

            break;

    }

	//Bottom boundary
    switch(wb){

        case NO_SLIP:

            for(i = 0; i <= imax; i++){

                U[k][i][0] = -U[k][i][1];
                V[k][i][0] = 0;

            }

            break;

        case FREE_SLIP:

            for(i = 0; i <= imax; i++){

                U[k][i][0] = U[k][i][1];
                V[k][i][0] = 0;

            }

            break;

        case OUTFLOW:

            for(i = 0; i <= imax; i++){

                U[k][i][0] = U[k][i][1];
                V[k][i][0] = V[k][i][1];

            }

            break;

    }

	//Obstacles 
	for(i = 1; i <= imax; i++){

		for(j = 1; j <= jmax; j++){
		
			switch(Flag[i][j]){

				case B_N:

					V[k][i][j] = 0;
					U[k][i][j] = -U[k][i][j+1];

					break;

				case B_S:

					V[k][i][j-1] = 0;
					U[k][i][j] = -U[k][i][j-1];

					break;

				case B_O:

					U[k][i][j] = 0;
					V[k][i][j] = -V[k][i+1][j];

					break;

				case B_W:

					U[k][i-1][j] = 0;
					V[k][i][j] = -V[k][i-1][j];

					break;

				case B_NW:

					V[k][i][j] = 0;
					U[k][i-1][j] = 0;
					U[k][i][j] = -U[k][i][j+1];

					break;

				case B_NO:

					U[k][i][j] = 0;
					V[k][i][j] = 0;

					break;

				case B_SW:

					U[k][i-1][j] = 0;
					V[k][i][j-1] = 0;
					U[k][i][j] = -U[k][i][j-1];
					V[k][i][j] = -V[k][i-1][j];

					break;

				case B_SO:

					U[k][i][j] = 0;
					V[k][i][j-1] = 0;
					V[k][i][j] = -V[k][i+1][j];

					break;

			}
			
		}

	}
	
}

void spec_boundary_val(
                       char *problem, 
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
                      ){

	int i, j;
    double U_wall = 1;

    //Inflow BC
	if(dP == 0){

		if(strcmp(problem, "flow_over_step") == 0){

            if (k == 0){

                for (j = jmax / 2 + 1; j <= jmax; j++){

                    U[k][0][j] = 1;
                    V[k][0][j] = 0;

                }

            }

            else{

                for (j = jmax / 2 + 1; j <= jmax; j++){

                    U[k][0][j] = 0;
                    V[k][0][j] = 0;

                }

            }


		} 

        else if(strcmp(problem, "karman_vortex_street") == 0){

            if (k == 0){

                for (j = 0; j <= jmax; j++){

                    U[k][0][j] = 1;
                    V[k][0][j] = 0;

                }

            }

            else{

                for (j = 0; j <= jmax; j++){

                    U[k][0][j] = 0;
                    V[k][0][j] = 0;

                }

            }


		} 

        else if(strcmp(problem, "plane_shear_flow") == 0){
    
            if (k == 0){

                for(j = 0; j <= jmax; j++){

                    U[k][0][j] = 1;
                    V[k][0][j] = 0;

                }

            }

            else{

                 for(j = 0; j <= jmax; j++){

                    U[k][0][j] = 0;
                    V[k][0][j] = 0;

                }

            }

		} 
    
        else if(strcmp(problem, "driven_cavity") == 0){

            if (k == 0){

                for(i = 0; i <= imax; i++){

                    U[k][i][jmax + 1] = 2 * U_wall - U[k][i][jmax];
                    V[k][i][jmax] = 0;

                }

            }

            else{

                for(i = 0; i <= imax; i++){

                    U[k][i][jmax + 1] = - U[k][i][jmax];
                    V[k][i][jmax] = 0;

                }

            }

		} 

	}
	
}

