
#ifndef GET_FIT
#define GET_FIT

#include <algorithm>
#include <cmath>

double get_fit (int nCond, 
                int nSignals, 
                const std::vector<std::vector<double>>& simResT0,
                const std::vector<std::vector<double>>& simResT1, 
                const std::vector<std::vector<double>>& cnolist0, 
                const std::vector<std::vector<double>>& cnolist1, 
                const std::vector<std::vector<int>>& interMatCut, 
                double sizeFac, 
                double NAFac,
                int nInTot) {

    int nInputs = 0;
    for (const auto& row : interMatCut) {
        nInputs += std::count(row.begin(), row.end(), -1);
    }

    int NAs = 0;
    for (const auto& row : simResT1) {
        for (const auto& elem : row) {
            if (std::isnan(elem)) {
                NAs += 1;
            }
        }
    }

    double NAPen = NAFac * NAs;
    int nDataPts = cnolist1.size() * cnolist1[0].size();

    double deviationPen = 0;
    int nDataP = 0;
    double r;

    for (unsigned i = 0; i < nCond; i++) {
	    for (unsigned j = 0; j < nSignals; j++) {
	        r =  simResT0[i][j] - cnolist0[i][j];
            if (!std::isnan(r)){
                deviationPen += r*r;
            }

            r = (simResT1[i][j] - cnolist1[i][j]);
            if (!std::isnan(r)){
                deviationPen += r*r;
            }
            if (!std::isnan(cnolist1[i][j])){
                nDataP += 1;
            }
        }
    }
    deviationPen *= 0.5;

    double sizePen = (double)(nDataPts*sizeFac*nInputs)/nInTot;
    double score = (deviationPen + NAPen + sizePen) / double(nDataP);

    return score;

}

#endif
