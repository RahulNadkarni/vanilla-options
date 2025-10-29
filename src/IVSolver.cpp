#include "IVSolver.h" 
#include "BlackScholes.h"
#include <cmath> 
#include <algorithm> 


double IVSolver::solveimpliedvolatility(
double marketprice, double S, double K, double R, double T, double d, 
bool isCall, double tolerance, int maxiterations){

    double lower = 0.01; 
    double upper = 5.0; 

    double vol = (lower + upper) / 2.0; 

    for (int i =0; i < maxiterations; i++){
        double pricediff = price_difference(vol, marketprice, S, K, R, T, d, isCall); 

        if (std::abs(pricediff) < tolerance){
            return vol; 
        } 

        if (pricediff > 0.0){
            upper = vol; 
        } else { 
            lower = vol; 
        }

        vol = (lower + upper) / 2.0;
    }
    
    return vol; 
}

double IVSolver::pricedifference(double vol, double marketprice, double S, double K, double R, double T, double d, bool isCall){
    duoble modelprice; 

    if (isCall){
        modelprice = bs_call(S, K, R, d, T, vol); 
    } else { 
        modelprice = bs_put(S, K, R, d, T, vol); 
    }

    return modelprice - marketprice;
}

double IVSolver::brent(std::function<double(double)> func, double lower, double upper, double tolerance, int maxiterations){
    double a = lower; 
    double b = upper; 
    double fa = func(a); 
    double fb = func(b); 
    if (fa * fb > 0.0){
        return (a+b)/2.0; 
    }

    double c= a; 
    double fc = fa; 
    bool mflag = true; 

    for (int i = 0; i < maxiterations; i++){
        if(std::abs(fb) < tolerance){
            return b; 
        } 

        if (std::abs(fa) < tolerance){ 
            return a; 
        }

        if (fa != fc && fb != fc) { 
            double s = (a * fb * fc) / ((fa-fb) * (fa-fc)) + 
            (b * fa * fc) / ((fb-fa) * (fb-fc)) + 
            (c * fa * fb) / ((fc-fa) * (fc-fb)); 
        } else { 
            double s = b - fb * (b-a) / (fb - fa); 
        }

        double checklower = (3 * a + b) / 4; 
        double checkupper = b; 

        if (s < checklower || s > checkupper){
            s = (a + b) / 2; 
            mflag = true; 
        } else { 
            mflag = false; 
        }

        double fs = func(s); 
        c = b; 
        fc = fb; 

        if (fa * fs < 0) { 
            b  = s; 
            fb = fs; 
        } else{ 
            a = s; 
            fa = fs; 
        }
        if (std::abs(fa) < std::abs(fb)){
            std::swap(a, b); 
            std::swap(fa, fb); 
        }
    }
    return b; 
}