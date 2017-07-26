//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2D Finite Difference Stochastic-Galerkin CFD solver                     //
//                                                                          //
//  Developer:      Jiho Yang (MEng)                                        //
//                  M.Sc. student, Computational Science & Engineering      //
//                                                                          //
//  Institution:    Technische Universitat Munchen                          //
//                  Department of Informatics                               //
//                  Chair of Scientific Computing                           //
//                                                                          //
//  Final update date: 08/08/2016                                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// NEEDS THROUGH CHECK - DIVERGES AT THE MOMENT

//Note that, unlike the 3D case, index convention for 3d matrices are as follows: [polynomial order][x][y]

//Check helper functions for memory allocation (segmentation error)

#include <stdio.h>
#include <sys/stat.h>
#include <dirent.h>
#include "helper.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.c"
#include "visual.h"
#include "parallel.h"
#include "stochastic_galerkin.h"

int main(int argn, char **args){

    //CPU time variables
    clock_t start_t, end_t, total_t; 
    
    //Fluid solver variables
    double ***U, ***V, ***P, ***F, ***G, ***RS;
    double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value, dP;
    int imax, jmax, itermax, it;
    double res = 0, t = 0, n = 0;
    int wl, wr, wt, wb;
    int **Flag;
    
    //Stochastic Galerkin variables
    double **F_sg1, **F_sg2, **G_sg1, **G_sg2;
    double nu_0, nu_1, lambda;
    int k, mode;
    
    //Directory variables
    char problem[60];
    char parameters_filename[60];
    char output_dirname[60];
    char old_output_filename[128];
    struct dirent *old_outputfile;
    DIR *output_dir;
    
    //Parallelisation variables
    int iproc, jproc, myrank, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, omg_i, omg_j, num_proc;
    double min_dt;
    double *bufSend, *bufRecv;

    //Start CPU time measurement
    start_t = clock();
    
    //Set up MPI environment
    MPI_Init(&argn, &args);
    
    //Get total number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    //Get number of rank (taskID of each process)
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);  
    
    //Read the name of the problem from the command line arguments
    if(argn > 1){
        
        strcpy(problem, args[1]);
        
    }
    
    else{
        
        printf("Please provide the name of the problem!\n e.g. Run ./sim <problem> for a problem based on problem.dat\n");
        
        MPI_Finalize();
        
        return 1;
        
    }
    
    //Generate input filename based on problem name
    strcpy(parameters_filename, problem);
    strcat(parameters_filename, ".dat");
    
    //Create a problem specific directory for saving the output
    strcpy(output_dirname, problem);
    strcat(output_dirname, "/");
    strcat(output_dirname, problem);
    mkdir(problem, 0777);
    output_dir = opendir(problem);
    
    //Delete existing files in the output folder
    while((old_outputfile = readdir(output_dir))){
        
        sprintf(old_output_filename, "%s/%s", problem, old_outputfile -> d_name);
        remove(old_output_filename);
        
    }
    
    //Read the problem configuration
    read_parameters(parameters_filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &dP, &wl, &wr, &wt, &wb, &iproc, &jproc, &mode);

    //Check if number of processes is appropriate and, if necessary, make changes accordingly
    check_np (&iproc, &jproc, num_proc, myrank);

    //Synchronise processes
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Determine subdomain and neighbours for each process
    init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &rank_l, &rank_r, &rank_b, &rank_t, &omg_i, &omg_j, num_proc);

    if (myrank == 0){
 
        printf("\n");
        printf("Simulation Start\n");
        printf("\n");
                                                                                                                                                                                 
    }
    
    //Define viscosity values (nu_0: average viscosity, nu_1: standard deviation of viscosity)
    nu_0 = 1/Re;
    nu_1 = nu_0 *0.3;
    lambda = nu_1 / nu_0;
    
    //Allocate memory
    U = matrix_3d(il - 2, ir + 1, jb - 1, jt + 1, 0, mode);
    V = matrix_3d(il - 1, ir + 1, jb - 2, jt + 1, 0, mode);
    P = matrix_3d(il - 1, ir + 1, jb - 1, jt + 1, 0, mode);
    F = matrix_3d(il - 2, ir + 1, jb - 1, jt + 1, 0, mode);
    G = matrix_3d(il - 1, ir + 1, jb - 2, jt + 1, 0, mode);
    RS= matrix_3d(il, ir, jb, jt, 0, mode);
    Flag = imatrix(il - 1, ir + 1, jb - 1, jt + 1);
    
    F_sg1 = matrix(il - 2, ir + 1, jb - 1, jt + 1);
    F_sg2 = matrix(il - 2, ir + 1, jb - 1, jt + 1);
    G_sg1 = matrix(il - 1, ir + 1, jb - 2, jt + 1);
    G_sg2 = matrix(il - 1, ir + 1, jb - 2, jt + 1);

    printf("U[%d][%d][%d] = %f\n", mode, ir, jt, U[mode][ir][jt]);
    printf("\n");
    
    //Allocate memory for buffers
    bufSend = malloc(max(ir - il + 3, jt - jb + 3) * sizeof(double));
    bufRecv = malloc(max(ir - il + 3, jt - jb + 3) * sizeof(double));
    
    //Initialize lower part of the domain with UI = 0 for the flow_over_step problem
    if(strcmp(problem, "flow_over_step") == 0){
        
        init_matrix_3d(U, il, ir, jb, min(jmax/2, jt), 0, mode, 0);
        
    }

    printf("il = %d,    ir = %d,    jb = %d,    jt = %d\n", il, ir, jb, jt);
    printf("\n");
    
    //Initialise flag field
    init_flag(problem, imax, jmax, il, ir, jb, jt, Flag);

    printf("flag initialised\n");
    printf("\n");
    
    //Initialise U, V, P
    init_uvp(UI, VI, PI, il, ir, jb, jt, U, V, P, mode);

    printf("flow field initialised\n");
    printf("\n");
    
    //Loop over time
    while(t <= t_end){
        
        //Compute dt
        calculate_dt(Re, tau, &dt, dx, dy, il, ir, jb, jt, U, V);

        printf("time step computed\n");
        printf("\n");

        printf("dt = %f\n", dt);

        MPI_Allreduce(&dt, &min_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        
        dt = min_dt;

        //Iterate over each polynomial
        for (k = 0; k <= mode; k++){

            printf("polynomial iteration started\n");
            
            //Set wall boundary conditions
            boundaryvalues(il, ir, jb, jt, imax, jmax, U, V, wl, wr, wt, wb, Flag, k, mode);

            printf("boundaryvalues done\n");
            
            //Set inflow boundary conditions
            spec_boundary_val(problem, il, ir, jb, jt, imax, jmax, U, V, P, Re, xlength, ylength, dP, k, mode);
            
            //Compute F(n) and G(n)
            calculate_fg(Re, GX, GY, alpha, dt, dx, dy, il, ir, jb, jt, imax, jmax, U, V, F, G, Flag, lambda, k, mode, F_sg1, F_sg2, G_sg1, G_sg2);
            
            //Compute the RHS of Pressure Poisson Equation
            calculate_rs(dt, dx, dy, il, ir, jb, jt, imax, jmax, F, G, RS, k, mode);
            
            //Perform SOR iterations
            it = 0;
            res = 1e6;
            
            while(it < itermax && res > eps){
                
                sor(omg, dx, dy, dP, il, ir, jb, jt, imax, jmax, rank_l, rank_r, rank_b, rank_t, P, RS, &res, Flag, bufSend, bufRecv, k, mode);
                it++;
                
            }
            
            //Update velocities
            calculate_uv(dt, dx, dy, il, ir, jb, jt, imax, jmax, U, V, F, G, P, Flag, k, mode);
            
            //Communicate velocity values
            uv_comm(U, V, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, k, mode);
            
            //Free the F_sg1,2 and G_sg1,2 matrices
            free_matrix(F_sg1, il - 2, ir + 1, jb - 1, jt + 1);
            free_matrix(F_sg2, il - 2, ir + 1, jb - 1, jt + 1);
            free_matrix(G_sg1, il - 1, ir + 1, jb - 2, jt + 1);
            free_matrix(G_sg2, il - 1, ir + 1, jb - 2, jt + 1);
            
            //Re-allocate memory for F_sg1,2 and G_sg1,2 matrices
            F_sg1 = matrix(il - 2, ir + 1, jb - 1, jt + 1);
            F_sg2 = matrix(il - 2, ir + 1, jb - 1, jt + 1);
            G_sg1 = matrix(il - 1, ir + 1, jb - 2, jt + 1);
            G_sg2 = matrix(il - 1, ir + 1, jb - 2, jt + 1);
            
            printf("Polynomial order k = %i finished for t = %f     |     Residual = %f\n", k, t, res);
            
        }
        
        t = t + dt;
        n++;

        if (myrank == 0){

            printf("======================================================================================================\n");
            printf(" Residual = %f     |      Timestep = %f     |       Number of SOR Iteration = %d\n", res, t, it);
            printf("======================================================================================================\n");

        }
        
        //Export solutions
        if((int) n % (int) dt_value == 0){
            
            write_vtkFile(output_dirname, myrank, n, xlength, ylength, il, ir, jb, jt, imax, jmax, dx, dy, U, V, P, mode);
            
        }
        
    }
    
    //Close the output folder
    closedir(output_dir);
    
    if(myrank == 0){
        
        printf("Please find the output in the folder \"%s\".\n", problem);
        
    }
    
    //Free allocated memory
    free_matrix_3d(U, il - 2, ir + 1, jb - 1, jt + 1, 0, mode);
    free_matrix_3d(V, il - 1, ir + 1, jb - 2, jt + 1, 0, mode);
    free_matrix_3d(P, il - 1, ir + 1, jb - 1, jt + 1, 0, mode);
    free_matrix_3d(F, il - 2, ir + 1, jb - 1, jt + 1, 0, mode);
    free_matrix_3d(G, il - 1, ir + 1, jb - 2, jt + 1, 0, mode);
    free_matrix_3d(RS, il, ir, jb, jt, 0, mode);

    free_matrix(F_sg1, il - 2, ir + 1, jb - 1, jt + 1);
    free_matrix(F_sg2, il - 2, ir + 1, jb - 1, jt + 1);
    free_matrix(G_sg1, il - 1, ir + 1, jb - 2, jt + 1);
    free_matrix(G_sg2, il - 1, ir + 1, jb - 2, jt + 1);
    free_imatrix(Flag, il - 1, ir + 1, jb - 1, jt + 1);
    free(bufSend);
    free(bufRecv);
    
    //Synchronise all processes	
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Finalise MPI
    MPI_Finalize();

    //Finish CPU time measurement
    end_t = clock();
    total_t = (long double)(end_t - start_t) / CLOCKS_PER_SEC;

    if (myrank == 0){
 
        printf("\n");
        printf("\n");
        printf("Total time taken by CPU: %lu\n", total_t  );
        printf("\n");
        printf("Exiting of the program...\n");
        printf("\n");

    }
        
    return 0;
    
}

