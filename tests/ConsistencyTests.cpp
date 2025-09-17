#include <iostream> 
#include <cmath> 
#include <cassert> 
#include <BlackScholes.h>

bool approx_equal(double a, double b, double tol=1e-6) { 
    return std::fabs(a-b) < tol; 
}

void test_put_call_parity(double S, double K, double r, double d, double T, double vol) { 
    double call = bs_call(S,K,r,d,T,vol); 
    double put = bs_put(S,K,r,d,T,vol); 
    double fwd = forward_price(S,K,r,d,T); 

    std::cout << "[Check] Put-Call Parity: C-P = F ..." << std::endl; 
    std::cout << "  LHS: " << (call - put) << "  RHS:  " << fwd<<std::endl; 
    assert(approx_equal(call - put, fwd, 1e-6)); 
}

void test_monotonicity_in_strike(double S, double r, double d, double T, double vol) {
    std::cout << "[Check] Call decreasing in strike..." << std::endl;
    double c1 = bs_call(S,90,r,d,T,vol); 
    double c2 = bs_call(S,100,r,d,T,vol); 
    double c3 = bs_call(S,110, r, d, T, vol); 
    assert(c1 >= c2 && c2 >= c3); 
}

void test_call_bounds(double S, double K, double r, double d, double T, double vol) {
    std::cout << "[Check] Call bounds... " << std::endl; 

    double call = bs_call(S,K,r,d,T,vol); 
    double upper = S; 
    double lower = std::max(0.0, S - K * std::exp(-r*T)); 
    assert(call <= upper + 1e-6); 
    assert(call >= lower - 1e-6); 
}

void test_volatility_monotonicity(double S, double K, double r, double d, double T){
    std::cout << "[Check] Call increasing in volatility..." << std::endl; 
    double c1 = bs_call(S,K,r,d,T, 0.1); 
    double c2 = bs_call(S,K,r,d,T, 0.3); 
    double c3 = bs_call(S,K,r,d,T,0.6);
    assert(c1 <= c2 && c2 <= c3); 
}
void test_time_monotonicity(double S, double K, double r, double T, double vol) {
    std::cout << "[Check] Call increasing in time (d=0)..." << std::endl; 
    double c1 = bs_call(S,K,r,0.0,0.5, vol); 
    double c2 = bs_call(S,K,r,0.0, 1.0, vol); 
    double c3 = bs_call(S,K,r,0.0, 2.0, vol); 
    
    assert(c1 <= c2 && c2 <= c3); 
}

void test_convexity_in_strike(double S, double r, double d, double T, double vol) {
    std::cout <<"[Check] Call convex in strike..." <<std::endl; 
    double c1 = bs_call(S,90,r,d,T,vol); 
    double c2 = bs_call(S,100,r,d,T,vol);
    double c3 = bs_call(S,110,r,d,T,vol);
    assert(c2 <= .5 * (c1 +c3) + 1e-6); 
}

void test_digital_relation(double S, double K, double r, double d, double T, double vol) {
    std::cout << "[Check] Digital-call + Digital-put = ZCB..." << std::endl;
    double dc = bs_digital_call(S,K,r,d,T,vol);
    double dp = bs_digital_put(S,K,r,d,T,vol);
    double bond = zero_coupon_bond(r,T);
    assert(approx_equal(dc + dp, bond, 1e-6));
}

int main() {
    double S=100, K=100, r=0.05, d=0.02, T=1.0, vol=0.2;

    test_put_call_parity(S,K,r,d,T,vol);
    test_monotonicity_in_strike(S,r,d,T,vol);
    test_call_bounds(S,K,r,d,T,vol);
    test_volatility_monotonicity(S,K,r,d,T);
    test_time_monotonicity(S,K,r,T,vol);
    test_convexity_in_strike(S,r,d,T,vol);
    test_digital_relation(S,K,r,d,T,vol);

    std::cout << "All consistency checks passed" << std::endl;
    return 0;
}