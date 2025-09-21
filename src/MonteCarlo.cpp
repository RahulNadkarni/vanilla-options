#include "MonteCarlo.h"
#include <cmath>
#include <iostream>

MonteCarloEngine::MonteCarloEngine(unsigned seed) : rng(seed), norm(0.0,1.0) {}

double MonteCarloEngine::simulate(const Option& opt, double S0, double r, double d, double vol, int nPaths) {
    double T = opt.expiry;
    double sum = 0.0, sumSq = 0.0;

    for (int i = 0; i < nPaths; i++) {
        double Z = norm(rng);
        double ST = S0 * std::exp((r - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        double payoff = opt.payoff(ST);
        double discPayoff = std::exp(-r * T) * payoff;

        sum += discPayoff;
        sumSq += discPayoff * discPayoff;
    }

    double mean = sum / nPaths;
    double var = (sumSq / nPaths) - (mean * mean);
    double stderr = std::sqrt(var / nPaths);

    std::cout << "MC Price: " << mean << " ± " << stderr << std::endl;
    return mean;
}

double MonteCarloEngine::delta_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps){
    double base = simulate(opt, S0, r, d, vol, nPaths); 
    double bumped = simulate(opt, S0 + eps, r,d, vol, nPaths); 
    return (bumped-base);
}

double MonteCarloEngine::vega_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps){
    double base = simulate(opt, S0, r, d, vol, nPaths); 
    double bumped = simulate(opt, S0, r, d, vol + eps, nPaths); 
    return (bumped-base) / eps; 
}

double MonteCarloEngine::gamma_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps){
    double up = simulate(opt, S0 + eps, r, d, vol, nPaths); 
    double mid = simulate(opt, S0, r, d, vol, nPaths); 
    double down = simulate(opt, S0 -  eps, r, d, vol, nPaths); 
    return (up - 2*mid + down)/(eps*eps); 
}

double MonteCarloEngine::delta_pathwise(double S0, double K, double r, double d, double vol, double T, int nPaths){
    double sum = 0.0; 
    for (int i = 0; i < nPaths; ++i){
        double Z = norm(rng); 
        double ST = S0 * std::exp((r-d-.5*vol*vol) * T + vol * std::sqrt(T) * Z); 
        double payoff_deriv = (ST > K > ST/S0 : 0.0); 
        sum+= std::exp(-r*T) * payoff_deriv;
    }

    return sum/nPaths; 
}

double MonteCarloEngine::delta_likelihood(double S0, double K, double r, double d, double vol, double T, int nPaths){
    double sum = 0.0; 
    for (int i = 0; i < nPaths; i++){
        double Z = norm(rng); 
        double ST = S0 * std::exp((r-d-.5*vol*vol) * T + vol*std::sqrt(T)*Z); 
        double payoff = std::max(ST-K, 0.0); 
        double score = (std::log(ST/S0) - (r-d-.5*vol*vol) *T) / (vol * vol * T * S0); 
        sum += std::exp(-r * T) * payoff * score; 
    }

    return sum /nPaths; 
}
