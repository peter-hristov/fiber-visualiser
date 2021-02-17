#pragma once

#include "./MeshTriangle.h"
#include <utility>
#include <vector>

namespace tv9k {
namespace utility {

// Changes per isovalue
class SurfaceMesh
{
  public:
    SurfaceMesh(){};

    std::vector<tv9k::utility::MeshTriangle> triangles;

    // Hold the object ID mask for a given caches isovalue
    std::vector<std::vector<std::vector<int>>> visited;

    void clear();

    // For a given object id compute, take all points on all triangles that belong
    // to it The get the bounding box of that as the min x, y and z and max x, y
    // and z of those points Then return the point that halfway between the min
    // and max point If no objects have this ID then return the origin (0, 0, 0)
    std::vector<GLfloat> getCenter(const int);

    std::pair<std::vector<GLfloat>, std::vector<GLfloat>> getMinMaxPointsObject(const int);
};
}
}
