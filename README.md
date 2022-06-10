# ACO-Cellnopt

## Purpose

This repo is hosting the development of the C++-based training of CellNOpt using parallel ACO 

## Content

* `/aco` folder contains the aco code
* `/benchmark` folder contains the case studies  
* `/R` utility functions for the R code
* `/src` C++ code
* `CMakeLists.txt` for building C++ code
* `Makefile` for running tests

## Compile 

To compile and run a toymodel, execute (sequential): 
```
make aco_toymodel
```
## Parallel execution

`CMakeLists.tx` should be modified to be
adapted to the computing platform and base 
software. The code has been tested in 
different infrastructures with different C 
compilers and different MPI libraries (impi, 
openmpi).

`parallel_job.sh` is an example script 
to launch  a parallel job with SLURM. 
Note that it should be adapted to the 
infrastructure at hand. 

Use:
```
sbatch ./parallel_job.sh
```
