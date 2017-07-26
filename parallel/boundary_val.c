#include "init.h"
#include "helper.h"
#include <string.h>
#include <stdio.h>

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
                   ){

	int i, j;
	
	//Left boundary
	if(il == 1){

		switch(wl){

			case NO_SLIP:

				V[k][0][jb-1] = -V[k][1][jb-1];

				for(j = jb; j <= jt; j++){

					U[k][0][j] = 0;
					V[k][0][j] = -V[k][1][j];

				}

				V[k][0][jt+1] = -V[k][1][jt+1];

				break;

			case FREE_SLIP:

				V[k][0][jb-1] = V[k][1][jb-1];

				for(j = jb; j <= jt; j++){

					U[k][0][j] = 0;
					V[k][0][j] = V[k][1][j];

				}

				V[k][0][jt+1] = V[k][1][jt+1];

				break;

			case OUTFLOW:

				V[k][0][jb-1] = V[k][1][jb-1];

				for(j = jb; j <= jt; j++){

					U[k][0][j] = U[k][1][j];
					V[k][0][j] = V[k][1][j];

				}

				V[k][0][jt+1] = V[k][1][jt+1];

				break;

		}

	}

	//Right boundary 
	if(ir == imax){

		switch(wr){

			case NO_SLIP:

				V[k][imax+1][jb-1] = -V[k][imax][jb-1];

				for(j = jb; j <= jt; j++){

					U[k][imax][j] = 0;
					V[k][imax+1][j] = -V[k][imax][j];

				}

				V[k][imax+1][jt+1] = -V[k][imax][jt+1];

				break;

			case FREE_SLIP:

				V[k][imax+1][jb-1] = V[k][imax][jb-1];

				for(j = jb; j <= jt; j++){

					U[k][imax][j] = 0;
					V[k][imax+1][j] = V[k][imax][j];

				}

				V[k][imax+1][jt+1] = V[k][imax][jt+1];

				break;

			case OUTFLOW:

				V[k][imax+1][jb-1] = V[k][imax][jb-1];

				for(j = jb; j <= jt; j++){

					U[k][imax][j] = U[k][imax-1][j];
					V[k][imax+1][j] = V[k][imax][j];

				}

				V[k][imax+1][jt+1] = V[k][imax][jt+1];

				break;

		}

	}
	
	//Top boundary
	if(jt == jmax){

		switch(wt){

			case NO_SLIP:

				U[k][il-1][jmax+1] = -U[k][il-1][jmax];

				for(i = il; i <= ir; i++){

					U[k][i][jmax+1] = -U[k][i][jmax];
					V[k][i][jmax] = 0;

				}

				U[k][ir+1][jmax+1] = -U[k][ir+1][jmax];

				break;

			case FREE_SLIP:

				U[k][il-1][jmax+1] = U[k][il-1][jmax];

				for(i = il; i <= ir; i++){

					U[k][i][jmax+1] = U[k][i][jmax];
					V[k][i][jmax] = 0;

				}

				U[k][ir+1][jmax+1] = U[k][ir+1][jmax];

				break;

			case OUTFLOW:

				U[k][il-1][jmax+1] = U[k][il-1][jmax];

				for(i = il; i <= ir; i++){

					U[k][i][jmax+1] = U[k][i][jmax];
					V[k][i][jmax] = V[k][i][jmax-1];

				}

				U[k][ir+1][jmax+1] = U[k][ir+1][jmax];

				break;

		}

	}

	//Bottom boundary
	if(jb == 1){

		switch(wb){

			case NO_SLIP:

				U[k][il-1][0] = -U[k][il-1][1];

				for(i = il; i <= ir; i++){

					U[k][i][0] = -U[k][i][1];
					V[k][i][0] = 0;

				}

				U[k][ir+1][0] = -U[k][ir+1][1];

				break;

			case FREE_SLIP:

				U[k][il-1][0] = U[k][il-1][1];

				for(i = il; i <= ir; i++){

					U[k][i][0] = U[k][i][1];
					V[k][i][0] = 0;

				}

				U[k][ir+1][0] = U[k][ir+1][1];

				break;

			case OUTFLOW:

				U[k][il-1][0] = U[k][il-1][1];

				for(i = il; i <= ir; i++){

					U[k][i][0] = U[k][i][1];
					V[k][i][0] = V[k][i][1];

				}

				U[k][ir+1][0] = U[k][ir+1][1];

				break;

		}

	}
	
	//Obstacles 
	for(i = il; i <= ir; i++){

		for(j = jb; j <= jt; j++){
		
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
                      ){

	int i, j;
    double U_wall = 1;

    //Inflow BC
	if(dP == 0){

		if(strcmp(problem, "flow_over_step") == 0){

			if(il == 1){

                if (k == 0){

                    for (j = jmax / 2 + 1; j <= jt; j++){

                        U[k][0][j] = 1;
                        V[k][0][j] = 0;

                    }

                }

                else{

                    for (j = jmax / 2 + 1; j <= jt; j++){

                        U[k][0][j] = 0;
                        V[k][0][j] = 0;

                    }

                }

			}

		} 

        else if(strcmp(problem, "karman_vortex_street") == 0){

			if(il == 1){

                if (k == 0){

                    for (j = jb; j <= jt; j++){

                        U[k][0][j] = 1;
                        V[k][0][j] = 0;

                    }

                }

                else{

                    for (j = jb; j <= jt; j++){

                        U[k][0][j] = 0;
                        V[k][0][j] = 0;

                    }

                }

			}

		} 

        else if(strcmp(problem, "plane_shear_flow") == 0){

			if(il == 1){
    
                if (k == 0){

                    for(j = jb; j <= jt; j++){

                        U[k][0][j] = 1;
                        V[k][0][j] = 0;

                    }

                }

                else{

                     for(j = jb; j <= jt; j++){

                        U[k][0][j] = 0;
                        V[k][0][j] = 0;

                    }

                }

			}

		} 
    
        else if(strcmp(problem, "driven_cavity") == 0){

			if(jt == jmax){

                if (k == 0){

                    for(i = max(1, il-1); i <= min(imax, ir+1); i++){

                        U[k][i][jmax + 1] = 2 * U_wall - U[k][i][jmax];
                        V[k][i][jmax] = 0;

                    }

                }

                else{

                    for(i = max(1, il-1); i <= min(imax, ir+1); i++){

                        U[k][i][jmax + 1] = - U[k][i][jmax];
                        V[k][i][jmax] = 0;

                    }

                }

			}

		} 

	}
	
}

