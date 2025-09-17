#pragma once 
#include "Options.h"
#include <random>

class EulerMonteCarloEngine { 
    private: 
        std::mt19937 rng; 
        std::normal_distribution<> norm; 

    public:
        EulerMonteCarloEngine(unsigned seed=42); 

        double simulate(const Option& opt, double S0, double r, double d, double vol, int nPaths, int nSteps); 
}; 
