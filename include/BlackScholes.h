#pragma once 

double forward_price(double S, double K, double r, double d, double T); 
double bs_call(double S, double K, double r, double d, double T, double vol); 
double bs_put(double S, double K, double r, double d, double T, double vol); 
double bs_digital_call(double S, double K, double r, double d, double T, double vol); 
double bs_digital_put(double S, double K, double r, double d, double T, double vol); 
double zero_coupon_bond(double r, double T); 