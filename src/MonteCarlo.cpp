#include "MonteCarlo.h"
#include <cmath>
#include <iostream>
#include <algorithm>

MonteCarloEngine::MonteCarloEngine(unsigned seed) : rng(seed), norm(0.0, 1.0) {}

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

double MonteCarloEngine::delta_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps) {
    double T = opt.expiry;
    double sum_base = 0.0, sum_bumped = 0.0;
    
    for (int i = 0; i < nPaths; i++) {
        double Z = norm(rng);
        
        double ST_base = S0 * std::exp((r - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        double payoff_base = opt.payoff(ST_base);
        sum_base += std::exp(-r * T) * payoff_base;
        
        double ST_bumped = (S0 + eps) * std::exp((r - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        double payoff_bumped = opt.payoff(ST_bumped);
        sum_bumped += std::exp(-r * T) * payoff_bumped;
    }
    
    double base = sum_base / nPaths;
    double bumped = sum_bumped / nPaths;
    return (bumped - base) / eps;
}

double MonteCarloEngine::vega_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps) {
    double T = opt.expiry;
    double sum_base = 0.0, sum_bumped = 0.0;
    
    for (int i = 0; i < nPaths; i++) {
        double Z = norm(rng);
        
        double ST_base = S0 * std::exp((r - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        double payoff_base = opt.payoff(ST_base);
        sum_base += std::exp(-r * T) * payoff_base;
        
        double vol_bumped = vol + eps;
        double ST_bumped = S0 * std::exp((r - d - 0.5 * vol_bumped * vol_bumped) * T + vol_bumped * std::sqrt(T) * Z);
        double payoff_bumped = opt.payoff(ST_bumped);
        sum_bumped += std::exp(-r * T) * payoff_bumped;
    }
    
    double base = sum_base / nPaths;
    double bumped = sum_bumped / nPaths;
    return (bumped - base) / eps;
}
double MonteCarloEngine::gamma_pathwise(double S0, double K, double r, double d, double vol, double T, int nPaths) {
    double sum = 0.0;
    
    for (int i = 0; i < nPaths; i++) {
        double Z = norm(rng);
        double ST = S0 * std::exp((r - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        
        if (ST > K) {  
            double h = S0 * 0.01; 
            double weight = std::exp(-0.5 * std::pow((ST - K) / h, 2)) / (h * std::sqrt(2 * M_PI));
            double gamma_contribution = weight * (ST / (S0 * S0));
            sum += std::exp(-r * T) * gamma_contribution;
        }
    }
    
    return sum / nPaths;
}

double MonteCarloEngine::gamma_fd(const Option& opt, double S0, double r, double d,
                                  double vol, int nPaths, double h_frac) {
    double T = opt.expiry;
    double h = h_frac * S0;     
    double sum_up = 0.0, sum_mid = 0.0, sum_down = 0.0;
    for (int i = 0; i < nPaths; i++) {
        double Z = norm(rng);
        
        double ST_mid = S0 * std::exp((r - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        sum_mid += std::exp(-r * T) * opt.payoff(ST_mid);
    
        double ST_up = (S0 + h) * std::exp((r - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        sum_up += std::exp(-r * T) * opt.payoff(ST_up);
        
        double ST_down = (S0 - h) * std::exp((r - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        sum_down += std::exp(-r * T) * opt.payoff(ST_down);
    }
    
    double price_up   = sum_up / nPaths;
    double price_mid  = sum_mid / nPaths;
    double price_down = sum_down / nPaths;
    
    return (price_up - 2.0 * price_mid + price_down) / (h * h);
}




double MonteCarloEngine::delta_pathwise(double S0, double K, double r, double d, double vol, double T, int nPaths) {
    double sum = 0.0;
    for (int i = 0; i < nPaths; ++i) {
        double Z = norm(rng);
        double ST = S0 * std::exp((r - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        double payoff_deriv = (ST > K) ? (ST / S0) : 0.0;
        sum += std::exp(-r * T) * payoff_deriv;
    }
    return sum / nPaths;
}

double MonteCarloEngine::delta_likelihood(double S0, double K, double r, double d, double vol, double T, int nPaths) {
    double sum = 0.0;
    for (int i = 0; i < nPaths; i++) {
        double Z = norm(rng);
        double ST = S0 * std::exp((r - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        double payoff = std::max(ST - K, 0.0);
        double score = Z / (S0 * vol * std::sqrt(T));
        sum += std::exp(-r * T) * payoff * score;
    }
    return sum / nPaths;
}

double MonteCarloEngine::rho_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps) {
    double T = opt.expiry;
    double sum_base = 0.0, sum_bumped = 0.0;
    
    for (int i = 0; i < nPaths; i++) {
        double Z = norm(rng);
        
        // Base case
        double ST_base = S0 * std::exp((r - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        double payoff_base = opt.payoff(ST_base);
        sum_base += std::exp(-r * T) * payoff_base;
        
        // Bumped case
        double r_bumped = r + eps;
        double ST_bumped = S0 * std::exp((r_bumped - d - 0.5 * vol * vol) * T + vol * std::sqrt(T) * Z);
        double payoff_bumped = opt.payoff(ST_bumped);
        sum_bumped += std::exp(-r_bumped * T) * payoff_bumped;
    }
    
    double base = sum_base / nPaths;
    double bumped = sum_bumped / nPaths;
    return (bumped - base) / eps;
}

double MonteCarloEngine::theta_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps) {
    double time_eps = eps * 0.01; 
    
    if (opt.expiry <= time_eps) return 0.0; 
    
    double T_base = opt.expiry;
    double T_shorter = opt.expiry - time_eps; 
    
    double sum_base = 0.0, sum_shorter = 0.0;
    
    for (int i = 0; i < nPaths; i++) {
        double Z = norm(rng);
        
        double ST_base = S0 * std::exp((r - d - 0.5 * vol * vol) * T_base + vol * std::sqrt(T_base) * Z);
        double payoff_base = opt.payoff(ST_base);
        sum_base += std::exp(-r * T_base) * payoff_base;
        
        double ST_shorter = S0 * std::exp((r - d - 0.5 * vol * vol) * T_shorter + vol * std::sqrt(T_shorter) * Z);
        double payoff_shorter = opt.payoff(ST_shorter);
        sum_shorter += std::exp(-r * T_shorter) * payoff_shorter;
    }
    
    double base = sum_base / nPaths;
    double shorter = sum_shorter / nPaths;
    
    return (shorter - base) / time_eps;
}