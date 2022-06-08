
#ifndef COMPUTESCORET1
#define COMPUTESCORET1

#include <algorithm>
#include <iostream>

#include "data.hpp"
#include "get_fit.hpp"
#include "simulator.hpp"



double compute_score_t1(const Data& data, 
                        std::vector<int>& bString, 
                        double sizeFac=0.0001, 
                        double NAFac=1, 
                        int timeIndex=2) {

    auto interMat = subset_cols(data.interMat, bString);
    
    auto finalCube = subset_rows(data.finalCube, bString);
    auto ixNeg = subset_rows(data.ixNeg, bString);
    auto ignoreCube = subset_rows(data.ignoreCube, bString);
    auto maxIx = subset_vector(data.maxIx, bString);

    auto nStimuli = data.nStimuli;
    auto nInhibitors = data.nInhibitors;
    auto nCond = data.nCond;
    auto nReacsCut = std::count(bString.begin(), bString.end(), 1);
    auto nReacs = data.nReacs;
    auto nSpecies = data.nSpecies;
    auto nMaxInputs = finalCube[0].size();

    auto indexSignals = data.indexSignals;
    auto indexStimuli = data.indexStimuli;
    auto indexInhibitors = data.indexInhibitors;
    auto nSignals = data.nSignals;

    auto valueInhibitors = data.valueInhibitors;
    auto valueStimuli = data.valueStimuli;

    auto cnolist0 = data.cnolistSignalsT0;
    auto cnolist1 = data.cnolistSignalsT1;

    int nInTot = 0;
    for (const auto& row : data.interMat) {
        nInTot += std::count(row.begin(), row.end(), -1);
    }

    auto simResults = simulator(nStimuli, nInhibitors, nCond, nReacsCut, nSpecies, nSignals, 
                                nMaxInputs, finalCube, ixNeg, ignoreCube, maxIx, indexSignals, 
                                indexStimuli, indexInhibitors, valueInhibitors, valueStimuli, 1);

    auto simResultsT0 = simulator(nStimuli, nInhibitors, nCond, nReacsCut, nSpecies, nSignals, 
                                  nMaxInputs, finalCube, ixNeg, ignoreCube, maxIx, indexSignals, 
                                  indexStimuli, indexInhibitors, valueInhibitors, valueStimuli, 0);

    auto simResultsT0Cut = select_cols(simResultsT0, indexSignals);
    auto simResultsCut = select_cols(simResults, indexSignals);
    double score = get_fit(nCond, nSignals, simResultsT0Cut, simResultsCut, 
                           cnolist0, cnolist1, interMat, sizeFac, NAFac, nInTot);
    
	return score;
}

#endif


