# aco-cellnopt

## Purpose
This repo is hosting the development of the C++-based training of CellNOpt using parallel ACO 

## Content

* `/aco` folder contains the aco code
* `/benchmark` folder contains the case studies
  * `/toymodel` simple network model with 12 nodes, 16 interactions
    * `ToyDataMMB.csv` experimental data for the model,
    * `ToyPKNMMB.sif` interaction network for the model,
    * `export_toymodel_to_hdf5.R` creates a CellNOpt-based model in R and exports to hdf5-file,
    * `toymodel.h5` hdf5 representation of the model,
    * `calibrate_toymodel_R.R` calibrattes the CellNOpt model using genetic algorithm,
    * `toymodel_calibration_R.txt` calibration results including time, optimal solution and objective function. 
* `/R` utility functions for the R code
  * `export_model_to_hdf5.R` exports a CellNopt-based model to hdf5 format
* `/src` C++ code
  * `compute_score_t1.hpp` calls `simulator` and `get_fit` functions
  * `data.hpp` class that reads the HDF5 files generated in R
  * `get_fit.hpp` function that calculates the score 
  * `simulator.hpp` simulates given a bit string
  * `vector_funcs.hpp` helper functions for vector and matricies
* `CMakeLists.txt` for building C++ code
* `Makefile` for running tests

## Execute

To run a toymodel, execute: 
```
make aco_toymodel
```

