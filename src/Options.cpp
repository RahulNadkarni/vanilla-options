#include "./include/Options.h"

CallOption::CallOption(double K_, double T_) {
    K = K_; 
    expiry = T_; 
}
double CallOption::payoff(double ST) const {
    return std::max(ST - K, 0.0); 
}

PutOption::PutOption(double K_, double T_) {
    K = K_; 
    expiry = T_; 
}
double PutOption::payoff(double ST) const { 
    return std::max(K - ST, 0.0); 
}