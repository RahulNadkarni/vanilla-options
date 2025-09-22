#include "Options.h"
#include <algorithm>

CallOption::CallOption(double K_, double T_) : Option(K_, T_) {}

double CallOption::payoff(double ST) const {
    return std::max(ST - K, 0.0);
}

std::unique_ptr<Option> CallOption::clone() const {
    return std::make_unique<CallOption>(*this);
}

PutOption::PutOption(double K_, double T_) : Option(K_, T_) {}

double PutOption::payoff(double ST) const {
    return std::max(K - ST, 0.0);
}

std::unique_ptr<Option> PutOption::clone() const {
    return std::make_unique<PutOption>(*this);
}