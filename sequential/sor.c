#include "sor.h"
#include "init.h"
#include <math.h>

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
        ){

    int i,j;
    int C_F_count = 0;

    double rloc;
    double coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));

    //SOR iteration
    for(i = 1; i <= imax; i++){

        for(j = 1; j <= jmax; j++){

            if(Flag[i][j] == C_F){

                P[k][i][j] = (1.0 - omg) * P[k][i][j]
                           + coeff*((P[k][i+1][j] + P[k][i-1][j]) / (dx * dx) + (P[k][i][j+1] + P[k][i][j-1]) / (dy * dy) - RS[k][i][j]);

            }

        }

	}
	
	//Compute the residual 
	rloc = 0;

	for(i = 1; i <= imax; i++){

		for(j = 1; j <= jmax; j++){

			if(Flag[i][j] == C_F){

				rloc += ( (P[k][i+1][j]-2.0*P[k][i][j]+P[k][i-1][j])/(dx*dx) + ( P[k][i][j+1]-2.0*P[k][i][j]+P[k][i][j-1])/(dy*dy) - RS[k][i][j])*
						( (P[k][i+1][j]-2.0*P[k][i][j]+P[k][i-1][j])/(dx*dx) + ( P[k][i][j+1]-2.0*P[k][i][j]+P[k][i][j-1])/(dy*dy) - RS[k][i][j]);

				C_F_count++;

			}

		}

	}
	
	//Set residual
    rloc = rloc / C_F_count;
    rloc = sqrt(rloc);
	*res = rloc;

	//Neumann BC for lower and upper walls
    for(i = 1; i <= imax; i++){

        P[k][i][0] = P[k][i][1];
        P[k][i][jmax+1] = P[k][i][jmax];

    }

    //Dirichlet BC for pressure
	if(dP != 0){

        //Left wall
        if (k == 0){

            for(j = 0; j <= jmax; j++){

                P[k][0][j] = 2 * dP - P[k][1][j]; 

            } 

        }

        else{

            for(j = 0; j <= jmax; j++){

                P[k][0][j] = - P[k][1][j]; 

            } 

        }

        //Right wall
        for(j = 0; j <= jmax; j++){

            P[k][imax+1][j] = -P[k][imax][j];

        } 


	} 

    //Neumann BC for pressure 
    else{

        //Left and right wall
        for(j = 1; j <= jmax; j++){

            P[k][0][j] = P[k][1][j];
            P[k][imax+1][j] = P[k][imax][j];

        }

	}
	
	//Pressure BC for obstacles
	for(i = 1; i <= imax; i++){

		for(j = 1; j <= jmax; j++){

			switch(Flag[i][j]){

				case B_N:

					P[k][i][j] = P[k][i][j+1];

					break;

				case B_S:

					P[k][i][j] = P[k][i][j-1];

					break;

				case B_O:

					P[k][i][j] = P[k][i+1][j];

					break;

				case B_W:

					P[k][i][j] = P[k][i-1][j];

					break;

				case B_NW:

					P[k][i][j] = 0.5 * (P[k][i-1][j] + P[k][i][j+1]);

					break;

				case B_NO:

					P[k][i][j] = 0.5 * (P[k][i+1][j] + P[k][i][j+1]);

					break;

				case B_SW:

					P[k][i][j] = 0.5 * (P[k][i-1][j] + P[k][i][j-1]);

					break;

				case B_SO:

					P[k][i][j] = 0.5 * (P[k][i+1][j] + P[k][i][j-1]);

					break;

			}

		}

	}
	
}

