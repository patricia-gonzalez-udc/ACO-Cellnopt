
#ifndef DATA_HPP
#define DATA_HPP

#include <highfive/H5Easy.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

#include <string>
#include <vector>
#include "vector_funcs.hpp"

using std::vector;

class Data {

public:

    vector<vector<int>> interMat;

    vector<vector<int>> notMat;
    
    int nSpecies;

    int nSpeciesCompressed;

    int nReacs;

    int nStimuli;

    int nInhibitors;

    int nSignals;

    int nCond;

    int nMaxInputs;

    vector<vector<int>> valueInhibitors;

    vector<vector<int>> valueStimuli;

    vector<int> indexStimuli;

    vector<int> indexInhibitors;

    vector<int> indexSignals;

    vector<vector<int>> finalCube;

    vector<vector<int>> ixNeg;

    vector<vector<int>> ignoreCube;

    vector<int> maxIx;

    vector<vector<double>> cnolistSignalsT0;

    vector<vector<double>> cnolistSignalsT1;

    Data(const std::string& filename) {
        H5Easy::File file(filename, H5Easy::File::ReadOnly);

        interMat = H5Easy::load<vector<vector<int>>>(file, "model/interMat");

        notMat = H5Easy::load<vector<vector<int>>>(file, "model/notMat");

        nSpecies = H5Easy::load<int>(file, "model/namesSpeciesNum");

        nSpeciesCompressed = H5Easy::load<int>(file, "model/speciesCompressedNum");

        nReacs = H5Easy::load<int>(file, "model/reacIDNum");
         
        indexStimuli = H5Easy::load<vector<int>>(file, "indexList/stimulated");
        indexStimuli -= 1;
        nStimuli = indexStimuli.size();

        indexInhibitors = H5Easy::load<vector<int>>(file, "indexList/inhibited");
        indexInhibitors -= 1;
        nInhibitors = indexInhibitors.size();

        indexSignals = H5Easy::load<vector<int>>(file, "indexList/signals");
        indexSignals -= 1;
        nSignals = indexSignals.size();

        finalCube = H5Easy::load<vector<vector<int>>>(file, "simList/finalCube");
        finalCube -= 1;

        ignoreCube = H5Easy::load<vector<vector<int>>>(file, "simList/ignoreCube");

        ixNeg = H5Easy::load<vector<vector<int>>>(file, "simList/ixNeg");

        nMaxInputs = finalCube[0].size();

        maxIx = H5Easy::load<vector<int>>(file, "simList/maxIx");
        maxIx -= 1;

        valueInhibitors = H5Easy::load<vector<vector<int>>>(file, "cnolist/inhibitors");

        valueStimuli = H5Easy::load<vector<vector<int>>>(file, "cnolist/stimuli");
        nCond = valueStimuli.size();

        cnolistSignalsT0 = H5Easy::load<vector<vector<double>>>(file, "cnolist/signals/t0");
        cnolistSignalsT1 = H5Easy::load<vector<vector<double>>>(file, "cnolist/signals/t1");
    }
};

#endif
