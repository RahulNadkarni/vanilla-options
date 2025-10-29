#include "IVExporter.h"
#include <fstream>
#include <iomanip>

bool IVExporter::export_to_csv(const IVSurface& surface, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) return false;

    const auto& grid_m = surface.get_moneyness_grid();
    const auto& grid_t = surface.get_tenor_grid();
    const auto& data = surface.get_surface();

    file << "Moneyness";
    for (size_t j = 0; j < grid_t.size(); ++j) {
        file << ",T_" << std::fixed << std::setprecision(3) << grid_t[j];
    }
    file << "\n";

    for (size_t i = 0; i < grid_m.size(); ++i) {
        file << std::fixed << std::setprecision(3) << grid_m[i];
        for (size_t j = 0; j < grid_t.size(); ++j) {
            double v = (i < data.size() && j < data[i].size()) ? data[i][j] : 0.0;
            file << "," << std::fixed << std::setprecision(6) << v;
        }
        file << "\n";
    }

    return true;
}

bool IVExporter::export_points_to_csv(const std::vector<IVPoint>& points, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) return false;

    file << "moneyness,tenor,implied_vol,market_price,is_call\n";
    for (const auto& p : points) {
        file << std::fixed << std::setprecision(6)
             << p.moneyness << ","
             << p.tenor << ","
             << p.implied_vol << ","
             << p.market_price << ","
             << (p.is_call ? 1 : 0) << "\n";
    }

    return true;
}

bool IVExporter::export_surface_with_headers(const IVSurface& surface, const std::string& filename) {
    return export_to_csv(surface, filename);
}


