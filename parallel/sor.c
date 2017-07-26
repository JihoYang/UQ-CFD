#include "sor.h"
#include "init.h"
#include <math.h>
#include <mpi.h>
#include "parallel.h"

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
        ){

  int i,j;
  int C_F_count = 0, C_F_total;

  double rloc;
  double coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));

	//BC for lower wall
	if(jb == 1){

		for(i = il; i <= ir; i++){

			P[k][i][0] = P[k][i][1];

		}

	}

	//BC for upper wall 
	if(jt == jmax){

		for(i = il; i <= ir; i++){

			P[k][i][jmax+1] = P[k][i][jmax];

		}

	}
	
    //Dirichlet BC for pressure
	if(dP != 0){

        //Left wall
		if(il == 1){

            if (k == 0){

                for(j = jb; j <= jt; j++){

                    P[k][0][j] = 2 * dP - P[k][1][j]; 

                } 

            }

            else{

                for(j = jb; j <= jt; j++){

                    P[k][0][j] = - P[k][1][j]; 

                } 

            }

		}

        //Right wall
		if(ir == imax){

			for(j = jb; j <= jt; j++){

				P[k][imax+1][j] = -P[k][imax][j];

			} 

		}

	} 

    //Neumann BC for pressure 
    else{

        //Left wall
		if(il == 1){

			for(j = jb; j <= jt; j++){

				P[k][0][j] = P[k][1][j];

			}

		}

        //Right wall
		if(ir == imax){

			for(j = jb; j <= jt; j++){

				P[k][imax+1][j] = P[k][imax][j];

			}

		}

	}
	
	//Pressure BC for obstacles
	for(i = il - 1; i <= ir + 1; i++){

		for(j = jb - 1; j <= jt + 1; j++){

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

	//SOR iteration
	for(i = il; i <= ir; i++){

		for(j = jb; j <= jt; j++){

			if(Flag[i][j] == C_F){

				P[k][i][j] = (1.0 - omg) * P[k][i][j]
						     + coeff*((P[k][i+1][j] + P[k][i-1][j]) / (dx * dx) + (P[k][i][j+1] + P[k][i][j-1]) / (dy * dy) - RS[k][i][j]);

			}

		}

	}
	
	//Communicate pressure values
	pressure_comm(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, k, mode);

	//Compute the residual 
	rloc = 0;

	for(i = il; i <= ir; i++){

		for(j = jb; j <= jt; j++){

			if(Flag[i][j] == C_F){

				rloc += ( (P[k][i+1][j]-2.0*P[k][i][j]+P[k][i-1][j])/(dx*dx) + ( P[k][i][j+1]-2.0*P[k][i][j]+P[k][i][j-1])/(dy*dy) - RS[k][i][j])*
						( (P[k][i+1][j]-2.0*P[k][i][j]+P[k][i-1][j])/(dx*dx) + ( P[k][i][j+1]-2.0*P[k][i][j]+P[k][i][j-1])/(dy*dy) - RS[k][i][j]);

				C_F_count++;

			}

		}

	}

	MPI_Allreduce(&rloc, res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&C_F_count, &C_F_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	//Set residual
	*res = sqrt(*res/C_F_total);

}

