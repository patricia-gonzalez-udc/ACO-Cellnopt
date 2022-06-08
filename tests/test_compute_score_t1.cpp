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

std::vector<int> split(const std::string& input, char separator)
{
    std::istringstream ss(input);
    std::string token;

    std::vector<int> output;
    while(getline(ss, token, separator)) {
		output.push_back(std::stoi(token));
    }

    return output;
}

int main(int argc, char **argv) {

    if (argc < 4) throw std::runtime_error("Need three input files");

    auto d = Data(std::string(argv[1]));
    
    std::fstream fin;
    std::vector<int> bs;
    std::vector<double> scores;

    fin.open(argv[2]);
    if (fin.is_open()) {
        std::string temp;
        while (std::getline(fin, temp)) {
            bs = split(temp, ',');
            scores.push_back(compute_score_t1(d, bs));
        }
        fin.close();
    }

    std::cout << "Scores: \t" << scores;

    std::vector<double> r_scores;
    fin.open(argv[3]);
    if (fin.is_open()) {
        std::string temp;
        while (std::getline(fin, temp)) {
            r_scores.push_back(std::stod(temp));
        }
    }
    std::cout << "R scores: \t" << r_scores;
    
    scores -= r_scores;
    std::cout << "Diffs: \t" << scores;
    
}

