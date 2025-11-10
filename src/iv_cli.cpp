#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "IVSolver.h"

namespace {
double parse_double(const std::string& value) {
    size_t pos = 0;
    double result = std::stod(value, &pos);
    if (pos != value.size()) {
        throw std::invalid_argument("invalid number");
    }
    return result;
}
}

int main(int argc, char** argv) {
    if (argc != 8) {
        std::cerr << "usage: iv_solver <market_price> <spot> <strike> <rate> <tenor> <dividend> <call|put>\n";
        return 1;
    }

    try {
        double market_price = parse_double(argv[1]);
        double spot = parse_double(argv[2]);
        double strike = parse_double(argv[3]);
        double rate = parse_double(argv[4]);
        double tenor = parse_double(argv[5]);
        double dividend = parse_double(argv[6]);
        std::string type = argv[7];

        bool is_call;
        if (type == "call") {
            is_call = true;
        } else if (type == "put") {
            is_call = false;
        } else {
            std::cerr << "option type must be call or put\n";
            return 1;
        }

        double iv = IVSolver::solve_implied_volatility(market_price, spot, strike, rate, tenor, dividend, is_call);
        std::cout << "{\"impliedVol\":" << std::fixed << std::setprecision(8) << iv << "}\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        return 1;
    }
}

