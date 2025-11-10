#include "IVSolver.h"
#include "BlackScholes.h"
#include <algorithm>
#include <cmath>

double IVSolver::solve_implied_volatility(double market_price,
                                          double S,
                                          double K,
                                          double r,
                                          double T,
                                          double d,
                                          bool is_call,
                                          double tolerance,
                                          int max_iterations) {
    double lower = 1e-4;
    double upper = 5.0;
    double f_lower = price_difference(lower, market_price, S, K, r, T, d, is_call);
    double f_upper = price_difference(upper, market_price, S, K, r, T, d, is_call);

    for (int expand = 0; expand < 8 && f_lower * f_upper > 0.0; ++expand) {
        lower = std::max(lower * 0.5, 1e-6);
        upper *= 2.0;
        f_lower = price_difference(lower, market_price, S, K, r, T, d, is_call);
        f_upper = price_difference(upper, market_price, S, K, r, T, d, is_call);
    }

    if (f_lower * f_upper > 0.0) {
        return std::max(lower, tolerance);
    }

    double mid = 0.5 * (lower + upper);
    for (int i = 0; i < max_iterations; ++i) {
        mid = 0.5 * (lower + upper);
        double f_mid = price_difference(mid, market_price, S, K, r, T, d, is_call);

        if (std::abs(f_mid) < tolerance) {
            return std::max(mid, tolerance);
        }

        if (f_lower * f_mid < 0.0) {
            upper = mid;
            f_upper = f_mid;
        } else {
            lower = mid;
            f_lower = f_mid;
        }

        if (std::abs(upper - lower) < tolerance) {
            return std::max(mid, tolerance);
        }
    }

    return std::max(mid, tolerance);
}

double IVSolver::price_difference(double vol,
                                  double market_price,
                                  double S,
                                  double K,
                                  double r,
                                  double T,
                                  double d,
                                  bool is_call) {
    double price = is_call
        ? bs_call(S, K, r, d, T, vol)
        : bs_put(S, K, r, d, T, vol);
    return price - market_price;
}