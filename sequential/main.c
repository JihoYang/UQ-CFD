//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  2D Finite Difference Stochastic-Galerkin CFD solver                     //
//                                                                          //
//  Developer:      Jiho Yang (MEng)                                        //
//                  M.Sc. student, Computational Science & Engineering      //
//                                                                          //
//  Institution:    Technische Universität München                          //
//                  Department of Computer Science                          //
//                  Chair of Scientific Computing                           //
//                                                                          //
//  Final update date: 12/08/2016                                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

//Note that, unlike the 3D case, index convention for 3d matrices are as follows: [polynomial order][x][y]

//Gaussian distribution used for Reynolds number (and hence viscosity). No other distributions are possible but can be done with minor change in C_klm functions

//ISSUES: Diverges for some Re > 100
//        Problems with obstacles diverges for Re > 10

//        Particularly high residual observed for cases with obstacles

#include <stdio.h>
#include <sys/stat.h>
#include <dirent.h>
#include "helper.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.c"
#include "visual.h"
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

    //Start CPU time measurement
    start_t = clock();
    
    //Read the name of the problem from the command line arguments
    if(argn > 1){
        
        strcpy(problem, args[1]);
        
    }
    
    else{
        
        printf("Please provide the name of the problem!\n e.g. Run ./sim <problem> for a problem based on problem.dat\n");
        
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
    read_parameters(parameters_filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &dP, &wl, &wr, &wt, &wb, &mode);
 
    printf("\n");
    printf("Simulation Start\n");
    printf("\n");
    
    //Define viscosity values (nu_0: average viscosity, nu_1: standard deviation of viscosity)
    nu_0 = 1 / Re;
    nu_1 = nu_0 * 0.1;
    lambda = nu_1 / nu_0;
    
    //Allocate memory
    U = matrix_3d(0, imax,     0, jmax + 1, 0, mode);
    V = matrix_3d(0, imax + 1, 0, jmax,     0, mode);
    P = matrix_3d(0, imax + 1, 0, jmax + 1, 0, mode);
    F = matrix_3d(0, imax,     0, jmax + 1, 0, mode);
    G = matrix_3d(0, imax + 1, 0, jmax,     0, mode);
    RS= matrix_3d(0, imax + 1, 0, jmax + 1, 0, mode);

    Flag = imatrix(0, imax + 1, 0, jmax + 1);
    
    F_sg1 = matrix(0, imax,     0, jmax + 1);
    F_sg2 = matrix(0, imax,     0, jmax + 1);
    G_sg1 = matrix(0, imax + 1, 0, jmax);
    G_sg2 = matrix(0, imax + 1, 0, jmax);

    //Initialize lower part of the domain with UI = 0 for the flow_over_step problem
    if(strcmp(problem, "flow_over_step") == 0){
        
        init_matrix_3d(U, 0, imax, 0, jmax/2, 0, mode, 0);
        
    }
    
    //Initialise flag field
    init_flag(problem, imax, jmax, Flag);
    
    //Initialise U, V, P
    init_uvp(UI, VI, PI, U, V, P, imax, jmax, mode);
    
    //Loop over time
    while(t <= t_end){
        
        //Compute dt
        calculate_dt(Re, tau, &dt, dx, dy, U, V, imax, jmax);

        //Iterate over each polynomial
        for (k = 0; k <= mode; k++){
            
            //Set wall boundary conditions
            boundaryvalues(imax, jmax, U, V, wl, wr, wt, wb, Flag, k, mode);
            
            //Set inflow boundary conditions
            spec_boundary_val(problem, imax, jmax, U, V, P, Re, xlength, ylength, dP, k, mode);
            
            //Compute F(n) and G(n)
            calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, Flag, lambda, k, mode, F_sg1, F_sg2, G_sg1, G_sg2);
            
            //Compute the RHS of Pressure Poisson Equation
            calculate_rs(dt, dx, dy, imax, jmax, F, G, RS, k, mode);
            
            //Perform SOR iterations
            it = 0;
            res = 1e6;
            
            while(it < itermax && res > eps){
                
                sor(omg, dx, dy, dP, imax, jmax, P, RS, &res, Flag, k, mode);
                it++;
                
            }
            
            //Update velocities
            calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, Flag, k, mode);
            
            //Free the F_sg1,2 and G_sg1,2 matrices
            free_matrix(F_sg1, 0, imax,     0, jmax + 1);
            free_matrix(F_sg2, 0, imax,     0, jmax + 1);
            free_matrix(G_sg1, 0, imax + 1, 0, jmax);
            free_matrix(G_sg2, 0, imax + 1, 0, jmax);
            
            //Re-allocate memory for F_sg1,2 and G_sg1,2 matrices
            F_sg1 = matrix(0, imax,     0, jmax + 1);
            F_sg2 = matrix(0, imax,     0, jmax + 1);
            G_sg1 = matrix(0, imax + 1, 0, jmax);
            G_sg2 = matrix(0, imax + 1, 0, jmax);
       
            printf("Polynomial order k = %i finished for t = %f     |     Residual = %f     |    Number of SOR Iteration = %d\n", k, t, res, it);
            
        }
        
        t = t + dt;
        n++;

        printf("\n");
        printf("======================================================================================================\n");
        printf(" Residual = %f     |      Timestep = %f     |       Number of SOR Iteration = %d\n", res, t, it);
        printf("======================================================================================================\n");
        printf("\n");
        
        //Export solutions
        if((int) n % (int) dt_value == 0){
            
            write_vtkFile(output_dirname, n, xlength, ylength, imax, jmax, dx, dy, U, V, P, mode);
            
        }
        
    }
    
    //Close the output folder
    closedir(output_dir);
    
    printf("Please find the output in the folder \"%s\".\n", problem);
    
    //Free allocated memory
    free_matrix_3d(U,  0, imax,     0, jmax + 1, 0, mode);
    free_matrix_3d(V,  0, imax + 1, 0, jmax,     0, mode);
    free_matrix_3d(P,  0, imax + 1, 0, jmax + 1, 0, mode);
    free_matrix_3d(F,  0, imax,     0, jmax + 1, 0, mode);
    free_matrix_3d(G,  0, imax + 1, 0, jmax,     0, mode);
    free_matrix_3d(RS, 0, imax + 1, 0, jmax + 1, 0, mode);
    
    free_imatrix(Flag, 0, imax + 1, 0, jmax + 1);
    free_matrix(F_sg1, 0, imax,     0, jmax + 1);
    free_matrix(F_sg2, 0, imax,     0, jmax + 1);
    free_matrix(G_sg1, 0, imax + 1, 0, jmax);
    free_matrix(G_sg2, 0, imax + 1, 0, jmax);

    //Finish CPU time measurement
    end_t = clock();
    total_t = (long double)(end_t - start_t) / CLOCKS_PER_SEC;
 
    printf("\n");
    printf("\n");
    printf("Total time taken by CPU: %lu\n", total_t  );
    printf("\n");
    printf("Exiting of the program...\n");
    printf("\n");
        
    return 0;
    
}

