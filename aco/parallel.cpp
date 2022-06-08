/**
 * @file parallel.c
 * @author Patricia
 * @brief File containing general functions about the parallelization of the 
 * algorithm with MPI.
 */

#include <signal.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <time.h>
#include "aco.hpp"

MPI_Request request[100][100];
MPI_Request Srequest;
MPI_Request SRrequest;
MPI_Status status;

int mpi_id;
int NPROC;

double **comm_solution;
double *best_solution_to_send;

void init_parallel ( void )
{
 
    MPI_Init(NULL,NULL);
 
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
    MPI_Comm_size(MPI_COMM_WORLD, &NPROC);

}


/**
 Routine to start recv. of communications from Colonies
 **/

void startComm ( void )
{
    int i;
    
    /* Allocation of the communication buffers */
    if ((comm_solution = (double **) malloc(sizeof(double *)*NPROC))==NULL) {
        printf("Comm Solution. Out of memory, exit.");
        exit(1);
    }
    for ( i = 0 ; i < NPROC ; i++ )
        comm_solution[i] = (double *) malloc(sizeof(double)*(n+1));
    
    if ((best_solution_to_send = (double *) malloc(sizeof(double)*(n+1)))==NULL) {
        printf("Best solution to send. Out of memory, exit.");
        exit(1);
    }

    
    /* Prepare to receive the best solutions found in
     the rest of the colonies */
    for (i=0 ; i<NPROC ; i++)
        if (i!=mpi_id)
   		MPI_Irecv(&comm_solution[i][0], n+1, MPI_DOUBLE, i,
        	      3000, MPI_COMM_WORLD, &request[i][ntry]);
    
}


/**
 Routine to send the best solution found to the rest of Colonies
 **/
void sendBest ( void )
{
    int     i;
    
    for ( i = 0; i < n; i++)
        best_solution_to_send[i] = (double) best_so_far_ant->solution[i];
    best_solution_to_send[n] = best_so_far_ant->score;

 // printf(">>>>> >>>> >>>>> ntry %d iter %d ID: %d, isend %f \n",ntry,iteration,mpi_id,best_so_far_ant->score);
   
   /* Send best solution found to the rest of the colonies */
    
    for( i = 0 ; i < NPROC ; i++ )
       if ( i != mpi_id ){
           MPI_Isend(best_solution_to_send, n+1, MPI_DOUBLE, i,
                     3000, MPI_COMM_WORLD, &Srequest);
        }
}

/**
 Routine to check if there are pending messages from Colonies
 **/
void listen()
{
    int         i, j;
    int         listen;
    
    /* Loop to listen best solutions from colonies */
    for( i = 0; i < NPROC ; i++) {
        
      if( i != mpi_id ){
          
        listen = 1;
          
        while ( listen == 1 ) {
            MPI_Test(&request[i][ntry], &listen, &status);
        
            if ( listen ) {
              //printf(" ntry %d iter %d ID: %d Recibido best: %f, del %d <<<<<<<<   <<<<<<<   <<<<<<<<\n",ntry,iteration,
	      	//		mpi_id,comm_solution[status.MPI_SOURCE][n], status.MPI_SOURCE);
            
              if ( comm_solution[status.MPI_SOURCE][n] < best_so_far_ant->score) {
	                
			best_so_far_ant->score = comm_solution[status.MPI_SOURCE][n];
        	        for ( j = 0 ; j < n ; j++ )
                	    best_so_far_ant->solution[j] = (int) comm_solution[status.MPI_SOURCE][j];
        
        		if ( mmas_flag ) {
          			 trail_max = 1. / ( (rho) * best_so_far_ant->score );
           	 	 	trail_min = trail_max / ( 2. * n );
           		 	trail_0 = trail_max;
       		 	}

                	if ( report ) fprintf(report,"%f \t %f\n",best_so_far_ant->score,elapsed_time(REAL));
                	if ( report_iter ) fprintf(report_iter,"%f \t %d\n",best_so_far_ant->score,iteration);
              }
             
            /* Prepare to receive more solutions */
               MPI_Irecv(&comm_solution[status.MPI_SOURCE][0], n+1, MPI_DOUBLE,
                      status.MPI_SOURCE, 3000, MPI_COMM_WORLD, &request[i][ntry]);
            
            }
        }
      }
    }
    
}

/**
 Routine to clear up pending messages in the asynchronous version
 **/
void exit_parallel ( void )
{
    int     i, listen, flag;
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* Reception of lost pending messages */
    for( i = 0 ; i < NPROC ; i++ ) {
        if( i != mpi_id ){
            listen = 1;
	    while ( listen ) {
	   	 MPI_Test(&request[i][ntry], &listen, &status);
           	 if ( listen ){
                    //printf("Mensaje perdido de %d  ID; %d \n",status.MPI_SOURCE,mpi_id);
                    MPI_Irecv(&comm_solution[status.MPI_SOURCE][0], n+1, MPI_DOUBLE,
                            status.MPI_SOURCE, 3000, MPI_COMM_WORLD, &request[i][ntry]);
                 }
	    }
	MPI_Cancel(&request[i][ntry]);
        }

    }	
    MPI_Barrier(MPI_COMM_WORLD);

    for ( i = 0 ; i < NPROC ; i++ ) 
        free( comm_solution[i] );
    free(comm_solution);
    free(best_solution_to_send);
}
       

/***************************  Final MPI-Report ************************************/

void init_parallel_report( void )
/*
 FUNCTION:       prepare report file,
 INPUT:          none
 OUTPUT:         none
 COMMENTS:
 */
{
    char temp_buffer[LINE_BUF_LEN];

    if (mpi_id == 0 ) { 
    	sprintf(temp_buffer,"conv_report.%dx%d",NPROC,omp_get_max_threads());
    	report = fopen(temp_buffer, "w");
    	sprintf(temp_buffer,"conv_report_iter.%dx%d",NPROC,omp_get_max_threads());
    	report_iter = fopen(temp_buffer, "w");
    	sprintf(temp_buffer,"final_report.%dx%d",NPROC,omp_get_max_threads());
    	final_report = fopen(temp_buffer, "w");
    	sprintf(temp_buffer,"results_report.%dx%d",NPROC,omp_get_max_threads());
    	results_report = fopen(temp_buffer, "w");
    }

}


/**
 Routine to write one summarize mpi report
 **/
void write_parallel_report ( void )
{
    int     i, best_com_id;
    double    best_com_score, com_score, best_com_time, com_time;
    int     com_iter, best_com_iter;
    
    best_com_score = best_local_score;
    best_com_iter = best_iteration;
    best_com_time = best_time;
    best_com_id = mpi_id;

    if( mpi_id == 0 ) {
        for( i=1 ; i<NPROC ; i++ ) {
              MPI_Recv(&com_score, 1, MPI_DOUBLE, i, 2000, MPI_COMM_WORLD, &status);
              MPI_Recv(&com_iter, 1, MPI_INT, i, 2000, MPI_COMM_WORLD, &status);
              MPI_Recv(&com_time, 1, MPI_DOUBLE, i, 2000, MPI_COMM_WORLD, &status);
            
              if ( com_score < best_com_score ) {
                    best_com_score = com_score;
                    best_com_iter = com_iter;
                    best_com_time = com_time;
                    best_com_id=i;
              }
	    
        }
 	if (final_report)
        	fprintf(final_report, "Try: %d \t Best: %.12f\t Iterations: %d\t Time %f\t Tot.time %f \t nproc: %d\n",
                	ntry, best_com_score, best_com_iter, best_com_time, elapsed_time( REAL ), best_com_id);
                fprintSolution(best_so_far_ant->solution);
		fflush(final_report);
		fflush(report_iter);
		fflush(report);
		fflush(results_report);
    } else {
        MPI_Isend(&best_com_score, 1, MPI_DOUBLE, 0, 2000, MPI_COMM_WORLD,&SRrequest);
        MPI_Isend(&best_com_iter, 1, MPI_INT, 0, 2000, MPI_COMM_WORLD,&SRrequest);
        MPI_Isend(&best_com_time, 1, MPI_DOUBLE, 0, 2000, MPI_COMM_WORLD,&SRrequest);
    }
   
}

