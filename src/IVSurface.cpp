#include "IVSurface.h"
#include "IVSolver.h"
#include <algorithm>
#include <cmath>

void IVSurface::add_iv_point(double S, double K, double r, double T, double d,
                             double market_price, bool is_call) {
    IVPoint point;
    point.moneyness = S / K;
    point.tenor = T;
    point.market_price = market_price;
    point.is_call = is_call;

    point.implied_vol = IVSolver::solve_implied_volatility(
        market_price, S, K, r, T, d, is_call
    );

    iv_points.push_back(point);
}

void IVSurface::build_surface(int moneyness_points, int tenor_points) {
    if (iv_points.empty() || moneyness_points <= 1 || tenor_points <= 1) {
        iv_surface.clear();
        moneyness_grid.clear();
        tenor_grid.clear();
        return;
    }

    double min_m = iv_points.front().moneyness;
    double max_m = iv_points.front().moneyness;
    double min_t = iv_points.front().tenor;
    double max_t = iv_points.front().tenor;

    for (const auto& p : iv_points) {
        min_m = std::min(min_m, p.moneyness);
        max_m = std::max(max_m, p.moneyness);
        min_t = std::min(min_t, p.tenor);
        max_t = std::max(max_t, p.tenor);
    }

    moneyness_grid.resize(moneyness_points);
    tenor_grid.resize(tenor_points);
    for (int i = 0; i < moneyness_points; ++i) {
        moneyness_grid[i] = min_m + (max_m - min_m) * (static_cast<double>(i) / (moneyness_points - 1));
    }
    for (int j = 0; j < tenor_points; ++j) {
        tenor_grid[j] = min_t + (max_t - min_t) * (static_cast<double>(j) / (tenor_points - 1));
    }

    iv_surface.assign(moneyness_points, std::vector<double>(tenor_points, 0.0));
    for (int i = 0; i < moneyness_points; ++i) {
        for (int j = 0; j < tenor_points; ++j) {
            iv_surface[i][j] = interpolate_iv(moneyness_grid[i], tenor_grid[j]);
        }
    }
}

double IVSurface::interpolate_iv(double moneyness, double tenor) const {
    if (iv_points.empty()) return 0.0;

    double total_weight = 0.0;
    double weighted_sum = 0.0;

    for (const auto& p : iv_points) {
        double dm = p.moneyness - moneyness;
        double dt = p.tenor - tenor;
        double dist2 = dm * dm + dt * dt;

        if (dist2 < 1e-12) {
            return p.implied_vol; // near-exact match
        }

        double w = 1.0 / dist2; // inverse distance squared
        total_weight += w;
        weighted_sum += w * p.implied_vol;
    }

    if (total_weight <= 0.0) return 0.0;
    return weighted_sum / total_weight;
}

void IVSurface::clear() {
    iv_points.clear();
    moneyness_grid.clear();
    tenor_grid.clear();
    iv_surface.clear();
}


