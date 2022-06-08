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

# prepare toy model

CNOlistToy <- CNOlist("./benchmark/toymodel/ToyDataMMB.csv")
ToyModel <- readSIF("./benchmark/toymodel/ToyPKNMMB.sif")

ToyModel_prep <- preprocessing(CNOlistToy, ToyModel, expansion=TRUE, compression=TRUE, cutNONC=TRUE, verbose=FALSE)

sink("./benchmark/toymodel/toymodel_calibration_R.txt")

print("Calibration time:")
system.time({
    ga_results <- gaBinaryT1(CNOlist = CNOlistToy,
                             model = ToyModel_prep,sizeFac = 0.0001,NAFac = 1,verbose = FALSE)
})

print("Optimal bitstring:")
print(ga_results$bString)

print("Optimal score:")
print(ga_results$bScore)

sink()