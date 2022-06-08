#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("CellNOptR", quietly = TRUE))
    BiocManager::install("CellNOptR")
if (!requireNamespace("rhdf5", quietly = TRUE))
    BiocManager::install("rhdf5")


library(CellNOptR)
library(rhdf5)

source("R/export_model_to_hdf5.R")

### 1. prepare toy model

CNOlistToy <- CNOlist("./benchmark/toymodel/ToyDataMMB.csv")
ToyModel <- readSIF("./benchmark/toymodel/ToyPKNMMB.sif")

ToyModel_prep <- preprocessing(CNOlistToy, ToyModel, expansion=TRUE, compression=TRUE, cutNONC=TRUE, verbose=FALSE)


# export
export_model_to_hdf5(CNOlistToy,ToyModel_prep,"toymodel.h5")
file.copy("./toymodel.h5",to = "benchmark/toymodel/toymodel.h5")
file.remove("./toymodel.h5")


