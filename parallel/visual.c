#include "helper.h"
#include "visual.h"
#include <stdio.h>

void write_vtkFile(
                   const char *szProblem,
                   int rank,
                   int timeStepNumber,
                   double xlength,
                   double ylength,
                   int il, 
                   int ir,
                   int jb, 
                   int jt,
                   int imax,
                   int jmax,
                   double dx,
                   double dy,
                   double ***U,
                   double ***V,
                   double ***P,
                   int mode
                  ){
              
    int i, j, k;

    char szFileName[80];
    FILE *fp=NULL;

    sprintf(szFileName, "%s_%i.%i.vtk", szProblem, rank, timeStepNumber);
    fp = fopen( szFileName, "w");

    if(fp == NULL){

        char szBuff[80];

        sprintf( szBuff, "Failed to open %s", szFileName );
        ERROR( szBuff );

        return;

    }
  
    write_vtkHeader(fp, ir - il + 1, jt - jb + 1, dx, dy);
    write_vtkPointCoordinates(fp, il, ir, jb, jt, dx, dy);

    fprintf(fp,"POINT_DATA %i \n", (ir-il+2)*(jt-jb+2) );

    for (k = 0; k <= mode; k++){

        fprintf(fp,"\n");
        fprintf(fp, "VECTORS velocity_%i float\n", k);

        for(j = jb - 1; j <= jt; j++){

            for(i = il - 1; i <= ir; i++){

                fprintf(fp, "%f %f 0\n", (U[k][i][j] + U[k][i][j+1]) * 0.5, (V[k][i][j] + V[k][i+1][j]) * 0.5);

            }

        }

    }

    fprintf(fp,"\n");
    fprintf(fp,"CELL_DATA %i \n", ((ir -il + 1)*(jt - jb + 1)));
    
    for (k = 0; k <= mode; k++){

        fprintf(fp, "SCALARS pressure_%i float 1 \n", k); 
        fprintf(fp, "LOOKUP_TABLE default \n");

        for(j = jb; j <= jt; j++){

            for(i = il; i <= ir; i++){

                fprintf(fp, "%f\n", P[k][i][j]);

            }

        }

    }

    if(fclose(fp)){

        char szBuff[80];

        sprintf( szBuff, "Failed to close %s", szFileName );
        ERROR( szBuff );

    }

}


void write_vtkHeader( 
                     FILE *fp, 
                     int imax, 
                     int jmax, 
                     double dx, 
                     double dy
                    ){

    if(fp == NULL){

        char szBuff[80];

        sprintf( szBuff, "Null pointer in write_vtkHeader" );
        ERROR( szBuff );

        return;

    }

    fprintf(fp,"# vtk DataFile Version 2.0\n");
    fprintf(fp,"2D Finite Difference Stochastic Galerkin solver output (written by Jiho Yang) \n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"\n");	
    fprintf(fp,"DATASET STRUCTURED_GRID\n");
    fprintf(fp,"DIMENSIONS  %i %i 1 \n", imax+1, jmax+1);
    fprintf(fp,"POINTS %i float\n", (imax+1)*(jmax+1) );
    fprintf(fp,"\n");

}


void write_vtkPointCoordinates( 
                               FILE *fp, 
                               int il, 
                               int ir, 
                               int jb, 
                               int jt, 
                               double dx, 
                               double dy
                              ){

    double originX = 0.0;  
    double originY = 0.0;

    int i = 0;
    int j = 0;

    for(j = jb - 1; j <= jt; j++){

        for(i = il - 1; i <= ir; i++){

            fprintf(fp, "%f %f 0\n", originX + (i * dx), originY +(j * dy));

        }

    }

}


