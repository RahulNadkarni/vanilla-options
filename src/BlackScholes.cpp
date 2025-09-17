#include "BlackScholes.h"
#include "Utils.h"
#include <cmath>

double forward_price(double S, double K, double r, double d, double T) {
    return S * std::exp(-d * T) - K * std::exp(-r * T); 
}

double bs_call(double S, double K, double r, double d, double T, double vol) { 
    double d1 = (std::log(S/K) + (r - d + 0.5 * vol * vol) *  T) / (vol * std::sqrt(T)); 
    double d2 = d1 - vol * std::sqrt(T); 
    return S * std::exp(-d * T) * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2); 
}

double bs_put(double S, double K, double r, double d, double T, double vol){
    double d1 = (std::log(S/K) + (r - d + 0.5 * vol * vol) *  T) / (vol * std::sqrt(T)); 
    double d2 = d1 - vol * std::sqrt(T); 
    return K * std::exp(-d * T) * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2); 
}


double bs_digital_call(double S, double K, double r, double d, double T, double vol) {
    double d2 = (std::log(S/K) + (r - d - 0.5 * vol * vol) * T) / (vol * std::sqrt(T)); 
    return std::exp(-r * T) * norm_cdf(d2); 
}

double bs_digital_put(double S, double K, double r, double d, double T, double vol){
    double d2 = (std::log(S/K) + (r - d - 0.5 * vol * vol) * T) / (vol * std::sqrt(T)); 
    return std::exp(-r * T) * norm_cdf(-d2); 
}

double zero_coupon_bond(double r, double T) { 
    return std::exp(-r * T); 
}
