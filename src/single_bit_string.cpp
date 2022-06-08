#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include "vector_funcs.hpp"
#include "data.hpp"
#include "compute_score_t1.hpp"

using std::cout;
using std::endl;


/* 
Here is a "template" for the optimization function.
This only calls the objective function with a vector of parameters and 
returns the objective function value. 
*/
double optim_algorithm(int n_vars,
                        double(*obj_function)(const Data&,std::vector<int>&,double,double,int),
                        const Data& user_data){

    // bs: bitString, vector of integers {0, 1}.
    std::vector<int> bs(n_vars,1);

    double score = obj_function(user_data, bs, 0.0001, 1, 2);
    return score;
}



/*
This function currently requires a single hdf5 input file that describes the model. 
Later, the optimization variables could be also inputs (or hard coded in the code).
*/


int main(int argc, char **argv) {

    if (argc < 1) throw std::runtime_error("Need an input file (hdf5 description of the model");

    // variable `model` stores the description of the model, 
    // all data that needed for cost function computation
    auto model = Data(std::string(argv[1]));
    
    std::vector<int> bs;
    double opt_results;

    // number of optimization variables for the optimizer
    double n_optim_variables = model.nReacs;

    opt_results = optim_algorithm(n_optim_variables, &compute_score_t1, model);

    cout << opt_results << endl;

}





