#pragma once

#include <string>
#include <vector>

#include "IVSurface.h"

class IVExporter {
public:
    static bool export_surface(const IVSurface& surface, const std::string& filename);
    static bool export_points(const std::vector<IVPoint>& points, const std::string& filename);
};

