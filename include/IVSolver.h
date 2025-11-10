#pragma once

class IVSolver {
public:
    static double solve_implied_volatility(double market_price,
                                           double S,
                                           double K,
                                           double r,
                                           double T,
                                           double d,
                                           bool is_call,
                                           double tolerance = 1e-6,
                                           int max_iterations = 200);

private:
    static double price_difference(double vol,
                                   double market_price,
                                   double S,
                                   double K,
                                   double r,
                                   double T,
                                   double d,
                                   bool is_call);
};

