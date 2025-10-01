#include <iostream>
#include <fstream>
#include <iomanip>
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
    
    // ============ MONTE CARLO OPTIMIZATION DEMONSTRATION ============
    // This section demonstrates the advanced Monte Carlo optimization techniques:
    // 1. Antithetic Sampling: Reduces variance by using paired random numbers (Z, -Z)
    // 2. Control Variates: Uses known expected values to reduce estimation error
    // 3. Parallel Processing: Leverages OpenMP for multi-threaded simulations
    // 4. Combined Optimization: Merges antithetic sampling with parallel processing
    std::cout << "\n=== MONTE CARLO OPTIMIZATION DEMONSTRATION ===" << std::endl;
    
    int opt_nPaths = 100000;
    std::cout << "\nTesting with " << opt_nPaths << " paths..." << std::endl;
    
    // Demonstrate individual optimization techniques
    std::cout << "\n--- Individual Optimization Techniques ---" << std::endl;
    
    MonteCarloEngine mc_opt(42);
    
    // 1. Standard Monte Carlo - baseline implementation
    std::cout << "\n1. Standard Monte Carlo:" << std::endl;
    double standard_price = mc_opt.simulate(call, S, r, d, vol, opt_nPaths);
    double standard_delta = mc_opt.delta_fd(call, S, r, d, vol, opt_nPaths, eps);
    
    // 2. Antithetic Sampling - variance reduction using paired random numbers (Z, -Z)
    std::cout << "\n2. Antithetic Sampling:" << std::endl;
    double antithetic_price = mc_opt.simulate_antithetic(call, S, r, d, vol, opt_nPaths);
    double antithetic_delta = mc_opt.delta_fd_antithetic(call, S, r, d, vol, opt_nPaths, eps);
    
    // 3. Control Variates - error reduction using known expected values
    std::cout << "\n3. Control Variates:" << std::endl;
    double cv_price = mc_opt.simulate_control_variate(call, S, r, d, vol, opt_nPaths);
    double cv_delta = mc_opt.delta_fd_control_variate(call, S, r, d, vol, opt_nPaths, eps);
    
    // 4. Parallel Processing - multi-threaded simulation using OpenMP
    std::cout << "\n4. Parallel Processing:" << std::endl;
    double parallel_price = mc_opt.simulate_parallel(call, S, r, d, vol, opt_nPaths);
    double parallel_delta = mc_opt.delta_fd_parallel(call, S, r, d, vol, opt_nPaths, eps);
    
    // 5. Combined Optimization - antithetic sampling + parallel processing
    std::cout << "\n5. Combined Optimization (Antithetic + Parallel):" << std::endl;
    double optimized_price = mc_opt.simulate_optimized(call, S, r, d, vol, opt_nPaths);
    double optimized_delta = mc_opt.delta_fd_optimized(call, S, r, d, vol, opt_nPaths, eps);
    
    // Display results comparison
    std::cout << "\n--- Results Comparison ---" << std::endl;
    std::cout << "Analytical BS Call Price: " << bs_call(S, K, r, d, T, vol) << std::endl;
    std::cout << "Analytical BS Delta: " << bs_call_delta(S, K, r, d, T, vol) << std::endl;
    std::cout << std::endl;
    
    std::cout << std::setw(25) << "Method" << std::setw(15) << "Price" << std::setw(15) << "Delta" << std::endl;
    std::cout << std::string(55, '-') << std::endl;
    std::cout << std::setw(25) << "Standard MC" << std::setw(15) << std::fixed << std::setprecision(6) << standard_price << std::setw(15) << standard_delta << std::endl;
    std::cout << std::setw(25) << "Antithetic" << std::setw(15) << antithetic_price << std::setw(15) << antithetic_delta << std::endl;
    std::cout << std::setw(25) << "Control Variate" << std::setw(15) << cv_price << std::setw(15) << cv_delta << std::endl;
    std::cout << std::setw(25) << "Parallel" << std::setw(15) << parallel_price << std::setw(15) << parallel_delta << std::endl;
    std::cout << std::setw(25) << "Optimized" << std::setw(15) << optimized_price << std::setw(15) << optimized_delta << std::endl;
    
    // Performance benchmarking - comprehensive timing and accuracy comparison
    std::cout << "\n--- Performance Benchmarking ---" << std::endl;
    int benchmark_paths = 500000;
    std::cout << "Running performance benchmark with " << benchmark_paths << " paths..." << std::endl;
    
    MonteCarloEngine mc_bench(42);
    mc_bench.compare_methods(call, S, r, d, vol, benchmark_paths);
    
    std::cout << "\n=== OPTIMIZATION SUMMARY ===" << std::endl;
    std::cout << "1. Antithetic Sampling: Reduces variance by using paired random numbers" << std::endl;
    std::cout << "2. Control Variates: Uses known expected values to reduce estimation error" << std::endl;
    std::cout << "3. Parallel Processing: Leverages multiple CPU cores for faster computation" << std::endl;
    std::cout << "4. Combined Optimization: Uses antithetic sampling with parallel processing" << std::endl;
    std::cout << "\nThe optimized methods provide better accuracy and/or faster computation!" << std::endl;
    
    return 0;
}
