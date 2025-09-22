#include <iostream>
#include <fstream>
#include "BlackScholes.h"
#include "Options.h"
#include "MonteCarlo.h"
#include "EulerMonteCarlo.h"
#include "Greeks.h"

int main() {
    double S = 100;
    double K = 105;
    double r = 0.05;
    double d = 0.0;
    double T = 1.0;
    double vol = 0.2;
    
    int nPaths = 500000;
    double eps = 1e-4;
    
    CallOption call(K, T);
    MonteCarloEngine mc(42);
    
    std::cout << "=== Analytical Prices ===" << std::endl;
    std::cout << "BS Call: " << bs_call(S, K, r, d, T, vol) << std::endl;
    std::cout << "BS Put : " << bs_put(S, K, r, d, T, vol) << std::endl;
    
    std::cout << "\n=== Monte Carlo Greeks (Finite Difference) ===" << std::endl;
    std::cout << "Delta: " << mc.delta_fd(call, S, r, d, vol, nPaths, eps) << std::endl;
    std::cout << "Gamma: " << mc.gamma_fd(call, S, r, d, vol, nPaths, eps) << std::endl;
    std::cout << "Vega : " << mc.vega_fd(call, S, r, d, vol, nPaths, eps) << std::endl;
    std::cout << "Rho : " << mc.rho_fd(call, S, r, d, vol, nPaths, eps) << std::endl;
    std::cout << "Theta: " << mc.theta_fd(call, S, r, d, vol, nPaths, eps) << std::endl;
    
    std::cout << "\n=== Analytical Greeks (Black–Scholes) ===" << std::endl;
    std::cout << "Delta: " << bs_call_delta(S, K, r, d, T, vol) << std::endl;
    std::cout << "Gamma: " << bs_call_gamma(S, K, r, d, T, vol) << std::endl;
    std::cout << "Vega : " << bs_call_vega(S, K, r, d, T, vol) << std::endl;
    std::cout << "Rho : " << bs_call_rho(S, K, r, d, T, vol) << std::endl;
    std::cout << "Theta: " << bs_call_theta(S, K, r, d, T, vol) << std::endl;
    
    std::cout << "\n=== Monte Carlo Greeks (Pathwise & LR) ===" << std::endl;
    std::cout << "Delta (Pathwise): " << mc.delta_pathwise(S, K, r, d, vol, T, nPaths) << std::endl;
    std::cout << "Delta (LR): " << mc.delta_likelihood(S, K, r, d, vol, T, nPaths) << std::endl;
    
    std::ofstream file("greeks.csv");
    file << "S,Delta_MC,Gamma_MC,Vega_MC,Rho_MC,Theta_MC,Delta_BS,Gamma_BS,Vega_BS,Rho_BS,Theta_BS\n";
    
    std::cout << "\nGenerating Greeks data for different spot prices..." << std::endl;
    
    int sweep_nPaths = 1000000;
    
    for (double s = 70; s <= 130; s += 5.0) {
        std::cout << "Processing S = " << s << std::endl;
        
        MonteCarloEngine mc_sweep(42);
        
        double delta_mc = mc_sweep.delta_fd(call, s, r, d, vol, sweep_nPaths, eps);
        double gamma_mc = mc_sweep.gamma_fd(call, s, r, d, vol, sweep_nPaths, eps);
        double vega_mc = mc_sweep.vega_fd(call, s, r, d, vol, sweep_nPaths, eps);
        double rho_mc = mc_sweep.rho_fd(call, s, r, d, vol, sweep_nPaths, eps);
        double theta_mc = mc_sweep.theta_fd(call, s, r, d, vol, sweep_nPaths, eps);
        
        double delta_bs = bs_call_delta(s, K, r, d, T, vol);
        double gamma_bs = bs_call_gamma(s, K, r, d, T, vol);
        double vega_bs = bs_call_vega(s, K, r, d, T, vol);
        double rho_bs = bs_call_rho(s, K, r, d, T, vol);
        double theta_bs = bs_call_theta(s, K, r, d, T, vol);
        
        file << s << "," 
             << delta_mc << "," << gamma_mc << "," << vega_mc << "," << rho_mc << "," << theta_mc << ","
             << delta_bs << "," << gamma_bs << "," << vega_bs << "," << rho_bs << "," << theta_bs << "\n";
    }
    file.close();
    
    std::cout << "\nSaved results to greeks.csv (for graphing)" << std::endl;
    std::cout << "CSV includes both Monte Carlo and analytical Greeks for comparison" << std::endl;
    
    return 0;
}
