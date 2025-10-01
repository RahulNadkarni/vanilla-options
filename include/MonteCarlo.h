#pragma once 
#include "Options.h"
#include <random> 
#include <vector>
#include <chrono>

/**
 * MonteCarloEngine - Advanced Monte Carlo simulation engine for option pricing
 * 
 * This class implements multiple Monte Carlo optimization techniques:
 * 1. Antithetic Sampling: Reduces variance by using paired random numbers (Z, -Z)
 * 2. Control Variates: Uses known expected values to reduce estimation error
 * 3. Parallel Processing: Leverages OpenMP for multi-threaded simulations
 * 4. Combined Optimization: Merges antithetic sampling with parallel processing
 * 
 * Performance improvements:
 * - Antithetic sampling: ~50-80% variance reduction
 * - Control variates: Significant error reduction when highly correlated
 * - Parallel processing: Near-linear speedup on multi-core systems
 * - Combined methods: Maximum performance with maintained accuracy
 */
class MonteCarloEngine { 
    private:
        std::mt19937 rng; 
        std::normal_distribution<> norm; 

        // Helper method to generate antithetic random numbers
        std::pair<double, double> generate_antithetic_pair();

    public: 
        MonteCarloEngine(unsigned seed = 42); 

        // Original methods
        double simulate(const Option& opt, double S0, double r, double d, double vol, int nPaths); 
        double delta_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4); 
        double gamma_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4); 
        double vega_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4); 
        double rho_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps);
        double theta_fd(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps);
        double gamma_pathwise(double S0, double K, double r, double d, double vol, double T, int nPaths);
        double delta_pathwise(double S0, double K, double r, double d, double vol, double T, int nPaths); 
        double delta_likelihood(double s0, double K, double r, double d, double vol, double T, int nPaths); 

        // Optimized methods with antithetic sampling
        double simulate_antithetic(const Option& opt, double S0, double r, double d, double vol, int nPaths);
        double delta_fd_antithetic(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4);
        double vega_fd_antithetic(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4);
        double gamma_fd_antithetic(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4);

        // Control variates methods
        double simulate_control_variate(const Option& opt, double S0, double r, double d, double vol, int nPaths);
        double delta_fd_control_variate(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4);

        // Parallel methods
        double simulate_parallel(const Option& opt, double S0, double r, double d, double vol, int nPaths, int num_threads = 0);
        double delta_fd_parallel(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4, int num_threads = 0);

        // Combined optimization methods
        double simulate_optimized(const Option& opt, double S0, double r, double d, double vol, int nPaths, int num_threads = 0);
        double delta_fd_optimized(const Option& opt, double S0, double r, double d, double vol, int nPaths, double eps=1e-4, int num_threads = 0);

        // Benchmarking methods
        struct BenchmarkResult {
            double value;
            double standard_error;
            double execution_time_ms;
            std::string method_name;
        };
        
        BenchmarkResult benchmark_simulate(const Option& opt, double S0, double r, double d, double vol, int nPaths, const std::string& method_name);
        void compare_methods(const Option& opt, double S0, double r, double d, double vol, int nPaths);
}; 