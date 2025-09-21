#pragma once 

double bs_call_delta(double S, double K, double r, double d, double T, double vol); 
double bs_call_gamma(double S, double K, double r, double d, double T, double vol); 
double bs_call_vega(double S, double K, double r, double d, double T, double vol); 
double bs_call_rho(double S, double K, double r, double d, double T, double vol); 
double bs_call_theta(double S, double K, double r, double d, double T, double vol); 

