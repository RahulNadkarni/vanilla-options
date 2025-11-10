#include "IVExporter.h"
#include <fstream>
#include <iomanip>

bool IVExporter::export_surface(const IVSurface& surface, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    const auto& grid_m = surface.get_moneyness_grid();
    const auto& grid_t = surface.get_tenor_grid();
    const auto& data = surface.get_surface();

    file << "moneyness";
    for (double tenor : grid_t) {
        file << ",T_" << std::fixed << std::setprecision(3) << tenor;
    }
    file << "\n";

    for (size_t i = 0; i < grid_m.size(); ++i) {
        file << std::fixed << std::setprecision(3) << grid_m[i];
        for (size_t j = 0; j < grid_t.size(); ++j) {
            double value = (i < data.size() && j < data[i].size()) ? data[i][j] : 0.0;
            file << "," << std::fixed << std::setprecision(6) << value;
        }
        file << "\n";
    }

    return true;
}

bool IVExporter::export_points(const std::vector<IVPoint>& points, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    file << "moneyness,tenor,implied_vol,market_price,is_call\n";
    for (const auto& point : points) {
        file << std::fixed << std::setprecision(6)
             << point.moneyness << ","
             << point.tenor << ","
             << point.implied_vol << ","
             << point.market_price << ","
             << (point.is_call ? 1 : 0) << "\n";
    }

    return true;
}
