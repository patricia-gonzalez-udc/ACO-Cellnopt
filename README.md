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

To compile and run a toymodel, execute: 
```
make aco_toymodel
```
## Parallel execution

`parallel_job.sh` is a script to launch 
a parallel job with SLURM. Note that it 
should be adapted to the infrastructure 
at hand.

Use:
```
sbatch ./parallel_job.sh
```
