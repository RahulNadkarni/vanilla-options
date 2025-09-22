#pragma once 
#include "Options.h"
#include <random> 

class MonteCarloEngine { 
    private:
        std::mt19937 rng; 
        std::normal_distribution<> norm; 

    public: 
        MonteCarloEngine(unsigned seed = 42); 

        double simulate(const Option& opt, double S0, double r, double d, double vol, int nPaths); 
        double delta_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4); 
        double gamma_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4); 
        double vega_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4); 
        double rho_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps);
        double theta_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps);
        double gamma_pathwise(double S0, double K, double r, double d, double vol, double T, int nPaths);
        double delta_pathwise(double S0, double K, double r, double d, double vol, double T, int nPaths); 
        double delta_likelihood(double s0, double K, double r, double d, double vol, double T, int nPaths); 
}; 