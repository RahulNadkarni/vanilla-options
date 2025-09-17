#include <iostream> 
#include "BlackScholes.h"
#include "Options.h"
#include "MonteCarlo.h"
#include "EulerMonteCarlo.h"
int main() { 
    double S = 100; 
    double K = 105; 
    double r = 0.05; 
    double d = 0.0; 
    double T = 1.0; 
    double vol = 0.2; 

    std::cout << "BS Call : " << bs_call(S,K,r,d,T,vol) << std::endl; 
    std::cout << "BS Put : " << bs_put(S,K,r,d,T,vol) << std::endl; 

    CallOption call(K , T); 
    MonteCarloEngine mc; 
    mc.simulate(call, S, r, d, vol, 100000); 

    EulerMonteCarloEngine emc; 
    emc.simulate(call, S, r, d, vol, 100000, 100); 

    return 0; 
}