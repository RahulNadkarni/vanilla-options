#pragma once

#include <vector>

struct IVPoint {
    double moneyness;
    double tenor;
    double market_price;
    bool is_call;
    double implied_vol;
};

class IVSurface {
public:
    void add_iv_point(double S,
                      double K,
                      double r,
                      double T,
                      double d,
                      double market_price,
                      bool is_call);

    void build_surface(int moneyness_points, int tenor_points);

    double interpolate_iv(double moneyness, double tenor) const;

    void clear();

    const std::vector<double>& get_moneyness_grid() const;
    const std::vector<double>& get_tenor_grid() const;
    const std::vector<std::vector<double>>& get_surface() const;
    const std::vector<IVPoint>& get_points() const;

private:
    std::vector<IVPoint> iv_points;
    std::vector<double> moneyness_grid;
    std::vector<double> tenor_grid;
    std::vector<std::vector<double>> iv_surface;
};

