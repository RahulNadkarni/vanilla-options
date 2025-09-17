#pragma once 
#include <algorithm>

class Option {
    public: 
        double expiry; 
        virtual double payoff(double ST) const = 0; 
        virtual ~Option() {}
}; 

class CallOption : public Option {
    public:
        double K; 
        CallOption(double K_, double T_); 
        double payoff(double ST) const override; 
}; 

class PutOption : public Option {
    public: 
        double K; 
        PutOption(double K_, double T_); 
        double payoff(double ST) const override; 
}; 
