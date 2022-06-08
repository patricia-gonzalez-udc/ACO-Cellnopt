/**
 * Ant Colony Optimization for the CellNopt
 *
 * @file cellnopt_parallel_aco.cpp
 * @author patricia.gonzalez@udc.es
 * @brief File contains aco main procedures
 *
 */

#ifdef OMP
   #include <omp.h>
#endif
#ifdef MPI2
   #include <mpi.h>
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "vector_funcs.hpp"
#include "data.hpp"
#include "compute_score_t1.hpp"
#include "aco.hpp"

using std::cout;
using std::endl;

int ntry;


/*************************    ACO procedures  *****************************/

int termination_condition( void )
/*    
      FUNCTION:       checks whether termination condition is met 
      INPUT:          none
      OUTPUT:         0 if condition is not met, number neq 0 otherwise
      (SIDE)EFFECTS:  none
*/
{
  return ( ((iteration >= max_iters) || (elapsed_time( REAL ) >= max_time)) ||
	  (best_so_far_ant->score <= optimal));
}



void construct_solutions( double(*obj_function)(const Data&,std::vector<int>&,double,double,int),
                         const Data& user_data)
/*    
      FUNCTION:       manage the solution construction phase
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution  
*/
{
    int k;        /* counter variable */
    double time1, time2;

//    time1 = elapsed_time(REAL);
    #pragma omp parallel for private(k)
    for ( k = 0 ; k < n_ants ; k++ ) {
        std::vector<int> bs(n,0);
        int step = 0;
        while ( step < n ) {
            select_gate( &ant[k], step);
            bs[step] = ant[k].solution[step];
            step++;
        }
        
        /* compute scores */
        ant[k].score = obj_function(user_data, bs, sizeFac, NAFac, timeIndex);
    }
//    time2 = elapsed_time(REAL);
//    printf("Time construction %f\n",time2-time1);
}


void init_ants( double(*obj_function)(const Data&,std::vector<int>&,double,double,int),
                         const Data& user_data)
/*
 FUNCTION:       manage the solution construction phase
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution
 */
{
    int k;        /* counter variable */
    double time1, time2;
    
    //    time1 = elapsed_time(REAL);
    /* solution for ants initialized randomly */
    #pragma omp parallel for private(k)
    for ( k = 0 ; k < n_ants ; k++ ) {
        std::vector<int> bs(n,0);
        int rnd;
        for ( int j = 0 ; j < n ; j++ ) {
            rnd = (int) round( ran01( &seed ) ); /* random number 0 or 1 */
            ant[k].solution[j] = rnd;
            bs[j] = rnd;
        }
        ant[k].score = obj_function(user_data, bs, sizeFac, NAFac, timeIndex);
    }
    //    time2 = elapsed_time(REAL);
    //    printf("Time construction %f\n",time2-time1);
}


void init_cellnopt_aco( void )
/*    
      FUNCTION: initilialize variables appropriately when starting a trial
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{
    int  i;

    /* Allocate ants */
    allocate_ants();

    /* Initialize variables concerning statistics etc. */
    iteration    = 1;
    best_iteration = 1;
    restart_best = 1;
    n_restarts = 0;
    best_so_far_ant->score = INFTY;
    
    start_timers();
    best_time = 0.0;
    time_used = elapsed_time( REAL );
    time_passed = time_used;
   
    /* allocate pheromone matrix */
    pheromone = generate_double_matrix( n, 2 );
    total = generate_double_matrix( n, 2 );

    /* Initialize pheromone trails */
    if (mmas_flag) {
        trail_max = 1. / ( (rho) * 0.5 );
        trail_min = trail_max / ( 2. * n );
        trail_0 = trail_max;
        init_pheromone_trails( trail_0 );
    }
    else {
        trail_0 = 0;
        init_pheromone_trails( trail_0 );
    }

    if (report) fprintf(report,"******** Try: %d **********\n",ntry);
    if (report_iter) fprintf(report_iter,"******** Try: %d **********\n",ntry);
   
    /* Allocate communication buffers */
    startComm ( );
 
}

void exit_cellnopt_aco( void )
/*
 FUNCTION: end trial
 INPUT:    none
 OUTPUT:   none
 COMMENTS: none
 */
{
 
    int     i;
    
#ifdef MPI2
    write_parallel_report ( );
    exit_parallel ( );
#else
    write_report ( );
#endif

    free( pheromone );
    free( total );
    for ( i = 0 ; i < n_ants ; i++ ) 
        free( ant[i].solution );
    free( ant );
    free( best_so_far_ant->solution );
    for ( i = 0 ; i < 2 ; i++ ) 
        free( offsprings[i].solution );
    free( offsprings );

}
    
void update_statistics( void )
/*    
      FUNCTION:       manage some statistical information about the trial, especially
                      if a new best solution (best-so-far or restart-best) is found and
                      adjust some parameters if a new best solution is found
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min 
                      and trail_max used by MMAS may be updated
*/
{

    int i, iteration_best_ant;

    iteration_best_ant = find_best(); /* iteration_best_ant is a global variable */

    if ( ant[iteration_best_ant].score < best_so_far_ant->score ) {
        
        time_used = elapsed_time( REAL ); /* best sol found after time_used */
        copy_from_to( &ant[iteration_best_ant], best_so_far_ant );
        
	/* Send best to the rest of the Colonies */
#ifdef MPI2
        sendBest ( );
#endif

        if ( report ) fprintf(report,"%f \t %f\n",best_so_far_ant->score,elapsed_time(REAL));
        if ( report_iter ) fprintf(report_iter,"%f \t %d\n",best_so_far_ant->score,iteration);
        
	best_iteration = iteration;
        restart_best = iteration;
        best_time = time_used;
        best_local_score = best_so_far_ant->score;

         if ( mmas_flag ) {
            trail_max = 1. / ( (rho) * best_so_far_ant->score );
            trail_min = trail_max / ( 2. * n );
            trail_0 = trail_max;
        }

    }
   
    if ( mmas_flag && (iteration - restart_best > restart_iters) ) {
        /* MAX-MIN Ant System was the first ACO algorithm to use
         pheromone trail re-initialisation as implemented
         here. Other ACO algorithms may also profit from this mechanism.
         */
	n_restarts++;
        
	init_pheromone_trails( trail_0 );
        compute_total_information();
        restart_best = iteration;
        restart_time = elapsed_time( REAL );
    }
    
}


void pheromone_trail_update( void )  
/*    
      FUNCTION:       manage global pheromone trail update for the ACO algorithms
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  pheromone trails are evaporated and pheromones are deposited 
 
*/
{
    int     k;
  
    /* Listen promising solutions from other Colonies */
#ifdef MPI2
    listen ( );
#endif

    /* Simulate the pheromone evaporation of all pheromones */
    evaporation();

    /* Apply the pheromone deposit for different ACO variants */
    if (mmas_flag) mmas_update();
    else if (eas_flag) eas_update();
    else as_update();

    /* Check pheromone trail limits for MMAS */
    if ( mmas_flag )
        check_pheromone_trail_limits();
    
  /* Compute combined information pheromone times heuristic info after
     the pheromone update  */
    compute_total_information();

}

void as_update( void )
/*
 FUNCTION:       manage global pheromone deposit for Ant System
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  all ants deposit pheromones on matrix "pheromone"
 */
{
    int   k;
    
    for ( k = 0 ; k < n_ants ; k++ )
        global_update_pheromone( &ant[k] );
    
}



void eas_update( void )
/*
 FUNCTION:       manage global pheromone deposit for Elitist Ant System
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  best ant so far deposit more pheromone
 */
{
    int   k;
    
    for ( k = 0 ; k < n_ants ; k++ )
        global_update_pheromone( &ant[k] );

    global_update_pheromone_weighted( best_so_far_ant, 2 );
    
}


void mmas_update( void )
/*
 FUNCTION:       manage global pheromone deposit for MAX-MIN Ant System
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  either the iteration-best or the best-so-far ant deposit pheromone
 on matrix "pheromone"
 */
{
    
    int iteration_best_ant;
   
    if ( iteration % u_gb ) {
        iteration_best_ant = find_best();
        global_update_pheromone( &ant[iteration_best_ant] );
    }
    else {
        global_update_pheromone( best_so_far_ant );
    }
    
    
    if ( ( iteration - restart_best ) < (int)(restart_iters/10) )
        u_gb = 10;
    else if ( (iteration - restart_best) < (int)(restart_iters/2) )
        u_gb = 5;
    else if ( (iteration - restart_best) < (int)(restart_iters/1.3) )
        u_gb = 3;
    else if ( (iteration - restart_best) < restart_iters )
        u_gb = 2;
    else
        u_gb = 1;

    
}


/* Genetic Algoritm as local improvement procedure for ACO*/
/* Only best ant per iteration is mated with best so far ant are replaced */
void GA_single (double(*obj_function)(const Data&,std::vector<int>&,double,double,int),
                        const Data& user_data)
/*
 FUNCTION:       GA to improve the ants solutions
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  best ant in the iteration is mate with the best ant so far
                    and the offsprings replace it if its score is improved
 */
{
    
    int     j, iteration_best_ant;
    int     *crossover_mask;
    int     best_os;
    std::vector<int> bs(n,0);

    crossover_mask  = (int*) calloc(n, sizeof(int));
    
    iteration_best_ant = find_best();
    
    /* crossover mask */
    for ( j = 0 ; j < n ; j++ )
         crossover_mask[j]  = round(ran01( &seed )); /* random number between 0 .. 1 */

    for ( j = 0 ; j < n ; j++ ) {
        if ( crossover_mask[j] ) {
            offsprings[0].solution[j] = ant[iteration_best_ant].solution[j];
            offsprings[1].solution[j] = best_so_far_ant->solution[j];
        }
	else {
            offsprings[1].solution[j] = ant[iteration_best_ant].solution[j];
            offsprings[0].solution[j] = best_so_far_ant->solution[j];
	}
    }
  
    /* evaluate offsprings */
    for ( j =0 ; j < n ; j++)
          bs[j] = (offsprings[0].solution[j]);
    offsprings[0].score = obj_function(user_data, bs, sizeFac, NAFac, timeIndex);
    for ( j =0 ; j < n ; j++)
          bs[j] = (offsprings[1].solution[j]);
    offsprings[1].score = obj_function(user_data, bs, sizeFac, NAFac, timeIndex);
    
    if( offsprings[0].score < offsprings[1].score ) best_os = 0;
    else best_os = 1;
        
    /* if offspring improves best iteration solution ->  replace */
    if ( offsprings[best_os].score < ant[iteration_best_ant].score ) {
        copy_from_to( &offsprings[best_os], &ant[iteration_best_ant] );
    }

    free(crossover_mask);
}

/* Genetic Algoritm as local improvement procedure for ACO*/
/* All ants are mated */
void GA_improvement(double(*obj_function)(const Data&,std::vector<int>&,double,double,int),
                        const Data& user_data)
/*
 FUNCTION:       GA to improve the ants solutions
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  ants are modified if crossover improves their scores
 */
{
 
    int     i, j, k;
    int     crossover_point;
    int     best_os;
    std::vector<int> bs(n,0);
    
    for ( i = 0 ; i < n_ants ; i++ ) {
        /* select ant to mate */
        k = (int) (ran01( &seed ) * (double) n_ants); /* random number between 0 .. n_ants-1 */
        
	/* crossover */
        crossover_point = (int) (ran01( &seed ) * (double) n); /* random number between 0 .. n-1 */
        
        for ( j = 0 ; j < n ; j++ ) {
            if ( j < crossover_point ) {
                offsprings[0].solution[j] = ant[i].solution[j];
                offsprings[1].solution[j] = ant[k].solution[j];
            }
            else {
                offsprings[1].solution[j] = ant[i].solution[j];
                offsprings[0].solution[j] = ant[k].solution[j];
            }
        }
	
        /* evaluate offsprings */
        for ( j =0 ; j < n ; j++)
                bs[j] = (offsprings[0].solution[j]);
        offsprings[0].score = obj_function(user_data, bs, sizeFac, NAFac, timeIndex);
        for ( j =0 ; j < n ; j++)
                bs[j] = (offsprings[1].solution[j]);
        offsprings[1].score = obj_function(user_data, bs, sizeFac, NAFac, timeIndex);
        
	if( offsprings[0].score < offsprings[1].score ) best_os = 0;
        else best_os = 1;
        
        /* if offsprings improves parents replace */
        if ( offsprings[best_os].score < ant[i].score ) {
            copy_from_to( &offsprings[best_os], &ant[i] );
        }
        else if ( offsprings[best_os].score < ant[k].score ) {
            copy_from_to( &offsprings[best_os], &ant[k] );
        }
        
        if ( offsprings[!best_os].score < ant[k].score ) {
            copy_from_to( &offsprings[!best_os], &ant[k] );
        }
        
    }

}


/*************************    ACO main    *********************************/

/*
 Here is a "template" for the optimization function.
 */
double aco_algorithm(int n_vars,
                     double(*obj_function)(const Data&,std::vector<int>&,double,double,int),
                     const Data& user_data){
    
    int     k;
    double  score;
    double  time1, time2;
   
    n = n_vars;
    
    init_cellnopt_aco();
    
    /* First iteration */
    init_ants(obj_function, user_data);
    update_statistics();
    pheromone_trail_update();
    
    /* next iterations */
    while ( !termination_condition() ) {
        iteration++;
        construct_solutions(obj_function, user_data);
        if (ga_flag) GA_single(obj_function, user_data);
        update_statistics();
        pheromone_trail_update();
    }
    
    score = best_so_far_ant->score;
    exit_cellnopt_aco();
    return(score);
}


/*************************    MAIN    *********************************/
/*
 This function currently requires a single hdf5 input file that describes the model.
 Later, the optimization variables could be also inputs (or hard coded in the code).
 */


int main(int argc, char **argv) {
    
    int  np, rank;
    
    /** MPI Initialization **/
#ifdef MPI2    
    init_parallel ( );
    init_parallel_report ( );
#else
    init_report ( );
#endif
    
    set_default_parameters ( );
    read_parameters ( );
    print_parameters ( );

    if (argc < 1) throw std::runtime_error("Need an input file (hdf5 description of the model");
    
    // variable `model` stores the description of the model,
    // all data that needed for cost function computation
    auto model = Data(std::string(argv[1]));
    double opt_results;
    double n_optim_variables = model.nReacs;
    
    for ( ntry = 0 ; ntry < max_tries ; ntry++ ) {
        opt_results = aco_algorithm(n_optim_variables, &compute_score_t1, model);
    }

#ifdef MPI2    
    MPI_Finalize();
#endif

}




