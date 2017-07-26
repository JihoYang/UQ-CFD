#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#define MASTER 0

#define     RIGHT   1
#define     LEFT    2
#define     UP      3
#define     DOWN    4

#define     VAR_U           1
#define     VAR_V           2
#define     VAR_P           3


void init_parallel(
                   int iproc,
                   int jproc,
                   int imax,
                   int jmax,
                   int *myrank,
                   int *il,
                   int *ir,
                   int *jb,
                   int *jt,
                   int *rank_l,
                   int *rank_r,
                   int *rank_b,
                   int *rank_t,
                   int *omg_i,
                   int *omg_j,
                   int num_proc
                  );


int check_np(
             int *iproc,
             int *jproc,
             int num_proc,
             int myrank
            );


void correct_ijproc(
                    int *iproc,
                    int *jproc,
                    int num_proc
                   );


void copy_to_buffer(
                    double ***matrix,
                    double *buffer,
                    int il,
                    int ir,
                    int jb,
                    int jt,
                    int k,
                    int mode
                   );


void copy_from_buffer(
                      double ***matrix,
                      double *buffer,
                      int il,
                      int ir,
                      int jb,
                      int jt,
                      int k,
                      int mode
                     );


void communicate(
                 double ***matrix,
                 int il,
                 int ir,
                 int jb,
                 int jt,
                 int rank_l,
                 int rank_r,
                 int rank_b,
                 int rank_t,
                 int direction,
                 int variable,
                 double *bufSend,
                 double *bufRecv,
                 int k,
                 int mode
                );


void pressure_comm(
                   double ***P,
                   int  il,
                   int ir,
                   int jb,
                   int jt,
                   int rank_l,
                   int rank_r,
                   int rank_b,
                   int rank_t,
                   double *bufSend,
                   double *bufRecv,
                   int k,
                   int mode
                  );


void uv_comm(
             double ***U,
             double ***V,
             int il,
             int ir,
             int jb,
             int jt,
             int rank_l,
             int rank_r,
             int rank_b,
             int rank_t,
             double *bufSend,
             double *bufRecv,
             int k,
             int mode
            );


int get_rank_location(
                      int omg_i,
                      int omg_j,
                      int iproc,
                      int jproc
                     );


//Produces a stderr text output 
void Program_Message(char *txt);

//Produces a stderr textoutput and synchronize all processes 
void Programm_Sync(char *txt);

//All processes will produce a text output, be synchronized and finished 
void Programm_Stop(char *txt);


