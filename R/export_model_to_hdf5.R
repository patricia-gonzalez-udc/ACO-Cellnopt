#' export CellNOpt data and model to hdf5 file
#' 
#' exports the model and data structure to an HDF5 file
#' 
#' @param cnolist CNOlist data structure, see CellNOptR::CNOlist()
#' @param model CellNOpt pkn model object, use CellNOptR::readSIF() and optionally CellNOptR::preprocessing(). 
#' @param h5_model_file file name to save the data to. 
#' 
export_model_to_hdf5 <- function(cnolist,model,h5_model_file){
    
    simList = CellNOptR::prep4sim(model)
    indexList = CellNOptR::indexFinder(cnolist, model)
    
    ### save all data generated
    h5_model <- h5_model_file
    rhdf5::h5createFile(h5_model)
    
    # 1. export cnolist data structure:
    h5createGroup(h5_model,group = "cnolist")
    
    ## these are matrices
    h5write(t(cnolist@cues), h5_model,"cnolist/cues")
    h5write(t(cnolist@inhibitors), h5_model,"cnolist/inhibitors")
    h5write(t(cnolist@stimuli), h5_model,"cnolist/stimuli")
    ## vector: 
    h5write(cnolist@timepoints, h5_model,"cnolist/timepoints")
    
    # signals
    # TODO: time index should be flexible
    h5createGroup(h5_model,group = "cnolist/signals")
    h5write(t(cnolist@signals[[1]]), h5_model,"cnolist/signals/t0")
    h5write(t(cnolist@signals[[2]]), h5_model,"cnolist/signals/t1")
    
    # variances
    # TODO: time index should be flexible
    h5createGroup(h5_model,group = "cnolist/variances")
    h5write(t(cnolist@variances[[1]]), h5_model,"cnolist/variances/t0")
    h5write(t(cnolist@variances[[2]]), h5_model,"cnolist/variances/t1")
    
    
    # 2. export the network data structure:
    h5createGroup(h5_model,group = "model")
    
    ## these are scalars: 
    h5write(length(model$reacID), h5_model,"model/reacIDNum")
    h5write(length(model$namesSpecies), h5_model,"model/namesSpeciesNum")
    h5write(length(model$speciesCompressed), h5_model,"model/speciesCompressedNum")
    
    ## matrices: 
    h5write(t(model$interMat), h5_model,"model/interMat")
    h5write(t(model$notMat), h5_model,"model/notMat")
    
    ## !! two variables that are varying length of lists are not exported. 
    # I think we wont need them - they contain information about the compression. 
    
    # 3. auxilary data
    #
    
    h5createGroup(h5_model,group = "simList")
    
    # matrices
    h5write(t(simList$finalCube), h5_model,"simList/finalCube")
    h5write(t(simList$ixNeg), h5_model,"simList/ixNeg")
    h5write(t(simList$ignoreCube), h5_model,"simList/ignoreCube")
    
    # vectors
    h5write(simList$maxIx, h5_model,"simList/maxIx")
    # constant
    h5write(simList$maxInput, h5_model,"simList/maxInput")
    
    
    h5createGroup(h5_model,group = "indexList")
    
    # vectors
    h5write(indexList$signals, h5_model,"indexList/signals")
    h5write(indexList$stimulated, h5_model,"indexList/stimulated")
    h5write(indexList$inhibited, h5_model,"indexList/inhibited")
    
    H5close()
    
}
