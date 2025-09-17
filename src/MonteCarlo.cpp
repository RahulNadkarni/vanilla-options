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