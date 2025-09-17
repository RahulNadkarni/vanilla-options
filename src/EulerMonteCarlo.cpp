#include "EulerMonteCarlo.h"
#include <cmath> 
#include <iostream> 

EulerMonteCarloEngine::EulerMonteCarloEngine(unsigned seed) : rng(seed), norm(0.0, 1.0) {}
double EulerMonteCarloEngine::simulate(const Option& opt, double S0, double r, double d, double vol, int nPaths, int nSteps) {
    double T = opt.expiry; 
    double dt = T /nSteps; 

    double sum = 0.0; 
    double sumSq = 0.0; 
    for (int i = 0; i < nPaths; i++){
        double S = S0; 
        for (int j = 0; j < nSteps; j++) {
            double Z = norm (rng); 
            S += (r-d) * S * dt + vol * S * std::sqrt(dt) * Z;
        }
        double payoff = opt.payoff(S); 
        double discPayoff = std::exp(-r * T) * payoff; 

        sum += discPayoff; 
        sumSq += discPayoff * discPayoff; 
    }
    double mean = sum / nPaths; 
    double var = (sumSq / nPaths) - (mean * mean); 
    double stderr = std::sqrt(var/nPaths); 

    std::cout << "Euler MC Price: " << mean << " ± " << stderr << " (steps= " << nSteps <<")" <<std::endl; 
    return mean; 
}