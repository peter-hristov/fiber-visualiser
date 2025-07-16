#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "./DisjointSet.h"

class TetMesh
{
  public:
    TetMesh() {}

    // Min/max range coordinates (for a bounding box)
    float minF, maxF;
    float minG, maxG;

    // Min/max domain coordinates (for a bounding box)
    float minX, maxX;
    float minY, maxY;
    float minZ, maxZ;

    // (Optional) indicative names for the two scalar fields and their units
    std::string longnameF, longnameG, units;

    // Tets, edges, triangles, etc
    std::vector<std::vector<size_t>> tetrahedra;
    std::map<std::pair<int, int>, int> edges;
    std::vector<std::set<int>> indexToTriangle;
    std::unordered_map<std::set<int>, int, MyHash<std::set<int>>> triangleToIndex;

    // What does this do again?
    std::vector<std::vector<int>> adjacentTrianglesIndex;

    std::vector<float> vertexCoordinatesF;
    std::vector<float> vertexCoordinatesG;
    std::vector<std::vector<float>> vertexDomainCoordinates;

    // Ideally we do want this with a proper data structure with a STAR, then there's no write conflits
    // Make sure to always keep the edge (u, v) such that u < v in index value,
    std::map<std::pair<int, int>, std::set<int>> upperLink;
    std::map<std::pair<int, int>, std::set<int>> lowerLink;

    // Tell us which triangles (as sets of IDs) are connected (part of a tetrahedron)
    std::set<std::pair<std::set<int>, std::set<int>>> connectedTriangles;

    void sortVertices();
    void printMesh();
    void perturbRangeValues(const float &perturbationEpsilon);
    void computeMinMaxRangeDomainCoordinates();
};

