#ifndef OPTIONS_H
#define OPTIONS_H

#include <memory>
#include <algorithm>

class Option {
public:
    double K;  
    double expiry;  

    Option(double K_, double T_) : K(K_), expiry(T_) {}
    virtual ~Option() = default;

    virtual double payoff(double ST) const = 0;
    virtual std::unique_ptr<Option> clone() const = 0;
};

class CallOption : public Option {
public:
    CallOption(double K_, double T_);
    double payoff(double ST) const override;
    std::unique_ptr<Option> clone() const override;
};

class PutOption : public Option {
public:
    PutOption(double K_, double T_);
    double payoff(double ST) const override;
    std::unique_ptr<Option> clone() const override;
};

#endif // OPTIONS_H