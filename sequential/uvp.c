#include "uvp.h"
#include <math.h>
#include <stdlib.h>
#include "helper.h"
#include "init.h"
#include "stochastic_galerkin.h"

//Adaptive time stepping
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
                ){

	int i, j;

	double u_max = 0, v_max = 0;
	double dt_min;

    //Find the maximum velocity
	if(tau > 0){

		for (i = 1; i <= imax; i++){

			for (j = 1; j <= jmax; j++){

				if (fabs(U[0][i][j]) > u_max){

					u_max = fabs(U[0][i][j]);

                }

				if (fabs(V[0][i][j]) > v_max){

					v_max = fabs(V[0][i][j]);

                }

			}

		}

        //Time step criterion
		dt_min = fmin(dy / v_max, dx / u_max);
		dt_min = fmin(Re * 0.5 * 1 / (1 / (dx * dx) + 1 / (dy * dy)), dt_min);

		*dt = tau * dt_min;

	}
	
}

//Calculate the F(n) and G(n) terms
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
                 ){

	int i, j, l, m;
    double C_klm;
    double C_kl1;
	
	//Left and right wall BC for F(n)
    for(j = 1; j <= jmax; j++){

        F[k][0][j] = U[k][0][j];
		F[k][imax][j] = U[k][imax][j];

    }

    //Boundary conditions for inner cell obstacles
    for (i = 1; i <= imax - 1; i++) {

        for (j = 1; j <= jmax; j++) {

            if(Flag[i][j] == B_W) {

                F[k][i-1][j] = U[k][i-1][j];

            }

            else if(Flag[i][j] == B_O) {

                F[k][i][j] = U[k][i][j];

            }

            else if(Flag[i][j] == B_NO){

                F[k][i][j] = U[k][i][j];

            }

            else if(Flag[i][j] == B_NW){

                F[k][i-1][j] = U[k][i-1][j];

            }

            else if(Flag[i][j] == B_SO){

                F[k][i][j] = U[k][i][j];

            }

            else if(Flag[i][j] == B_SW){

                F[k][i-1][j] = U[k][i-1][j];

            }

        }

    }
    
    //Compute the 1st new term for F
    C_kl1 = compute_C_kl1(k,0);
    
    for (i = 1; i <= imax - 1; i++){
        
        for (j = 1; j <= jmax; j++){
            
            F_sg1[i][j] = C_kl1 * ((U[0][i+1][j] - 2 * U[0][i][j] + U[0][i-1][j]) / (dx * dx) + (U[0][i][j+1] - 2 * U[0][i][j] + U[0][i][j-1]) / (dy * dy));
            
        }
        
    }
    
    for (l = 1; l <= mode; l++){
        
        C_kl1 = compute_C_kl1(k,l);
        
        for (i = 1; i <= imax - 1; i++){
            
            for (j = 1; j <= jmax; j++){
                
                F_sg1[i][j] = F_sg1[i][j] + C_kl1 * ((U[l][i+1][j] - 2 * U[l][i][j] + U[l][i-1][j]) / (dx * dx)
                                                     + (U[l][i][j+1] - 2 * U[l][i][j] + U[l][i][j-1]) / (dy * dy));
                
            }
            
        }
        
    }
    
    //Compute the 2nd new term for F
    C_klm = compute_C_klm(k,0,0);
    
    for (i = 1; i <= imax - 1; i++){
        
        for (j = 1; j <= jmax; j++){
            
            F_sg2[i][j] = C_klm * (
                                   1 / dx * (((U[0][i][j] + U[0][i+1][j]) * (U[0][i][j] + U[0][i+1][j]) / 4)
                                             - ((U[0][i-1][j] + U[0][i][j]) * (U[0][i-1][j] + U[0][i][j]) / 4))
                                   + alpha / dx * (fabs(U[0][i][j] + U[0][i+1][j]) * (U[0][i][j] - U[0][i+1][j]) / 4
                                                   - fabs(U[0][i-1][j] + U[0][i][j]) * (U[0][i-1][j] - U[0][i][j]) / 4)
                                   + 1 / dy * ((V[0][i][j] + V[0][i+1][j]) * (U[0][i][j] + U[0][i][j+1]) / 4
                                               - (V[0][i][j-1] + V[0][i+1][j-1]) * (U[0][i][j-1] + U[0][i][j]) / 4)
                                   + alpha / dy * (fabs(V[0][i][j] + V[0][i+1][j]) * (U[0][i][j] - U[0][i][j+1]) / 4
                                                   - fabs(V[0][i][j-1] + V[0][i+1][j-1]) * (U[0][i][j-1] - U[0][i][j]) / 4)
                                   
                                   );
            
        }
        
    }
    
    for (l = 0; l <= mode; l++){
        
        for (m = 1; m <= mode; m++){
            
            C_klm = compute_C_klm(k,l,m);
            
            for (i = 1; i <= imax - 1; i++){
                
                for (j = 1; j <= jmax; j++){
                    
                    F_sg2[i][j] = F_sg2[i][j] + C_klm * (
                                                         1 / dx * (((U[l][i][j] + U[l][i+1][j]) * (U[m][i][j] + U[m][i+1][j]) / 4)
                                                                   - ((U[l][i-1][j] + U[l][i][j]) * (U[m][i-1][j] + U[m][i][j]) / 4))
                                                         + alpha / dx * (fabs(U[l][i][j] + U[l][i+1][j]) * (U[m][i][j] - U[m][i+1][j]) / 4
                                                                         - fabs(U[l][i-1][j] + U[l][i][j]) * (U[m][i-1][j] - U[m][i][j]) / 4)
                                                         + 1 / dy * ((V[l][i][j] + V[l][i+1][j]) * (U[m][i][j] + U[m][i][j+1]) / 4
                                                                     - (V[l][i][j-1] + V[l][i+1][j-1]) * (U[m][i][j-1] + U[m][i][j]) / 4)
                                                         + alpha / dy * (fabs(V[l][i][j] + V[l][i+1][j]) * (U[m][i][j] - U[m][i][j+1]) / 4
                                                                         - fabs(V[l][i][j-1] + V[l][i+1][j-1]) * (U[m][i][j-1] - U[m][i][j]) / 4)
                                                         );
                    
                }
                
            }
            
        }
        
    }
   
    //Compute F(n)
	for (i = 1; i <= imax - 1; i++){

		for (j = 1; j <= jmax; j++){

            if(Flag[i][j] == C_F){

                F[k][i][j] = U[k][i][j] + dt * (
                                                + 1 / Re * ((U[k][i+1][j] - 2 * U[k][i][j] + U[k][i-1][j]) / (dx * dx) + (U[k][i][j+1] - 2 * U[k][i][j] + U[k][i][j-1]) / (dy * dy))
                                                + lambda / Re * F_sg1[i][j]
                                                - F_sg2[i][j]
                                                + GX
                                                );

			}

		}

	}

    //Bottom and top wall BC for G(n)
	for (i = 1; i <= imax; i++){

		G[k][i][0] = V[k][i][0];
		G[k][i][jmax] = V[k][i][jmax];

	}

    //Boundary conditions for inner cell obstacles
    for (i = 1; i <= imax; i++) {

        for (j = 1; j <= jmax - 1; j++) {

            if(Flag[i][j] == B_N) {

                G[k][i][j] = V[k][i][j];

            }

            else if(Flag[i][j] == B_S) {

                G[k][i][j-1] = V[k][i][j-1];

            }

            else if(Flag[i][j] == B_NO) {

                G[k][i][j] = V[k][i][j];

            }

            else if(Flag[i][j] == B_NW) {

                G[k][i][j] = V[k][i][j];

            }

            else if(Flag[i][j] == B_SO) {

                G[k][i][j-1] = V[k][i][j-1];

            }

            else if(Flag[i][j] == B_SW) {

                G[k][i][j-1] = V[k][i][j-1];

            }
        
        }

    }
    
    //Compute the 1st new term for G
    C_kl1 = compute_C_kl1(k,0);
    
    for (i = 1; i <= imax; i++){
        
        for (j = 1; j <= jmax - 1; j++){
            
            G_sg1[i][j] = C_kl1 * ((V[0][i+1][j] - 2 * V[0][i][j] + V[0][i-1][j]) / (dx * dx) + (V[0][i][j+1] - 2 * V[0][i][j] + V[0][i][j-1]) / (dy * dy));
            
        }
        
    }
    
    for (l = 1; l <= mode; l++){
        
        C_kl1 = compute_C_kl1(k,l);
        
        for (i = 1; i <= imax; i++){
            
            for (j = 1; j <= jmax - 1; j++){
                
                G_sg1[i][j] = G_sg1[i][j] + C_kl1 * ((V[l][i+1][j] - 2 * V[l][i][j] + V[l][i-1][j]) / (dx * dx)
                                                     + (V[l][i][j+1] - 2 * V[l][i][j] + V[l][i][j-1]) / (dy * dy));
                
            }
            
        }
        
    }
    
    //Compute the 2nd new term for G
    C_klm = compute_C_klm(k,0,0);
    
    for (i = 1; i <= imax; i++){
        
        for (j = 1; j <= jmax - 1; j++){
            
            G_sg2[i][j] = C_klm * (
                                   1 / dy * (((V[0][i][j] + V[0][i][j+1]) * (V[0][i][j] + V[0][i][j+1]) / 4)
                                             - ((V[0][i][j-1] + V[0][i][j]) * (V[0][i][j-1] + V[0][i][j]) / 4))
                                   + alpha / dy * (fabs(V[0][i][j] + V[0][i][j+1]) * (V[0][i][j] - V[0][i][j+1]) / 4
                                                   - fabs(V[0][i][j-1] + V[0][i][j]) * (V[0][i][j-1] - V[0][i][j]) / 4)
                                   + 1 / dx * ((V[0][i][j] + V[0][i+1][j]) * (U[0][i][j] + U[0][i][j+1]) / 4
                                               - (V[0][i-1][j] + V[0][i][j]) * (U[0][i-1][j] + U[0][i-1][j+1]) / 4)
                                   + alpha / dx * (fabs(V[0][i][j] + V[0][i+1][j]) * (U[0][i][j] - U[0][i][j+1]) / 4
                                                   - fabs(V[0][i-1][j] - V[0][i][j]) * (U[0][i-1][j] - U[0][i-1][j+1]) / 4)
                                   );
            
        }
        
    }
    
    for (l = 0; l <= mode; l++){
        
        for (m = 1; m <= mode; m++){
            
            C_klm = compute_C_klm(k,l,m);
            
            for (i = 1; i <= imax; i++){
                
                for (j = 1; j <= jmax - 1; j++){
                    
                    G_sg2[i][j] = G_sg2[i][j] + C_klm * (
                                                         1 / dy * (((V[l][i][j] + V[l][i][j+1]) * (V[m][i][j] + V[m][i][j+1]) / 4)
                                                                   - ((V[l][i][j-1] + V[l][i][j]) * (V[m][i][j-1] + V[m][i][j]) / 4))
                                                         + alpha / dy * (fabs(V[l][i][j] + V[l][i][j+1]) * (V[m][i][j] - V[m][i][j+1]) / 4
                                                                         - fabs(V[l][i][j-1] + V[l][i][j]) * (V[m][i][j-1] - V[m][i][j]) / 4)
                                                         + 1 / dx * ((V[m][i][j] + V[m][i+1][j]) * (U[l][i][j] + U[l][i][j+1]) / 4
                                                                     - (V[m][i-1][j] + V[m][i][j]) * (U[l][i-1][j] + U[l][i-1][j+1]) / 4)
                                                         + alpha / dx * (fabs(V[m][i][j] + V[m][i+1][j]) * (U[l][i][j] - U[l][i][j+1]) / 4
                                                                         - fabs(V[m][i-1][j] - V[m][i][j]) * (U[l][i-1][j] - U[l][i-1][j+1]) / 4)
                                                         );
                    
                }
                
            }
            
        }
        
    }

    //Calculate G(n)
	for (i = 1; i <= imax; i++){

		for (j = 1; j <= jmax - 1; j++){

            if(Flag[i][j] == C_F){

                G[k][i][j] = V[k][i][j] + dt * (
                                                + 1 / Re * ((V[k][i+1][j] - 2 * V[k][i][j] + V[k][i-1][j]) / (dx * dx) + (V[k][i][j+1] - 2 * V[k][i][j] + V[k][i][j-1]) / (dy * dy))
                                                + lambda / Re * G_sg1[i][j]
                                                - G_sg2[i][j]
                                                + GY
                                                );

            }

		}

	}
	
}

//Calculate the RHS of the Pressure Poisson Equation
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
                 ){

	int i, j;

	for (i = 1; i <= imax; i++){

		for (j = 1; j <= jmax; j++){

			RS[k][i][j] = 1 / dt * ((F[k][i][j] - F[k][i-1][j]) / dx + (G[k][i][j] - G[k][i][j-1]) / dy);

		}

	}

}

//Update the velocities
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
                ){

	int i, j;
    
    //Update the U velocities
	for (i = 1; i <= imax - 1; i++){

		for (j = 1; j <= jmax; j++){

			if(Flag[i][j] == C_F){

				U[k][i][j] = F[k][i][j] - dt * (P[k][i+1][j] - P[k][i][j]) / dx;

			} 

            else{

				U[k][i][j] = 0;

			}

		}

	}
	
    //Update the V velocities
	for (i = 1; i <= imax; i++){

		for (j = 1; j <= jmax - 1; j++){

			if(Flag[i][j] == C_F){

				V[k][i][j] = G[k][i][j] - dt * (P[k][i][j+1] - P[k][i][j]) / dy;

			} 

            else{

				V[k][i][j] = 0;

			}

		}

	}
	
}

