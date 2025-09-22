#include "Greeks.h"
#include "Utils.h"
#include <cmath>

static void d1d2(double S, double K, double r, double d, double T, double vol, double &d1, double &d2){
    double vsqrt = vol * std::sqrt(T); 
    d1 = (std::log(S/K) + (r - d + 0.5 * vol * vol) * T) / vsqrt; 
    d2 = d1 - vsqrt; 
}

double bs_call_delta(double S, double K, double r, double d, double T, double vol){
    double d1, d2;
    d1d2(S, K, r, d, T, vol, d1, d2); 
    return std::exp(-d*T) * norm_cdf(d1); 
}

double bs_call_gamma(double S, double K, double r, double d, double T, double vol){
    double d1, d2; 
    d1d2(S, K, r, d, T, vol, d1, d2); 
    return std::exp(-d*T) * norm_pdf(d1) / (S * vol * std::sqrt(T));
}

double bs_call_vega(double S, double K, double r, double d, double T, double vol){
    double d1, d2; 
    d1d2(S, K, r, d, T, vol, d1, d2); 
    return S * std::exp(-d*T) * norm_pdf(d1) * std::sqrt(T); 
}

double bs_call_rho(double S, double K, double r, double d, double T, double vol){
    double d1, d2; 
    d1d2(S, K, r, d, T, vol, d1, d2); 
    return K * T * std::exp(-r*T) * norm_cdf(d2); 
}

double bs_call_theta(double S, double K, double r, double d, double T, double vol){
    double d1, d2; 
    d1d2(S, K, r, d, T, vol, d1, d2); 

    double term1 = -(S * std::exp(-d*T) * norm_pdf(d1) * vol) / (2.0 * std::sqrt(T)); 
    double term2 = -r * K * std::exp(-r*T) * norm_cdf(d2); 
    double term3 = d * S * std::exp(-d*T) * norm_cdf(d1); 

    return term1 + term2 + term3; 
}
