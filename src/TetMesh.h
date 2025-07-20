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

    // Domain and range coordinates
    std::vector<float> vertexCoordinatesF;
    std::vector<float> vertexCoordinatesG;
    std::vector<std::vector<float>> vertexDomainCoordinates;

    // Bounding box min/max for the domain and range coordinates of all vertices
    float minF, maxF, minG, maxG;
    float minX, maxX, minY, maxY, minZ, maxZ;

    // (Optional) indicative names for the two scalar fields and their units
    std::string longnameF, longnameG, units;

    // Combinatorial structure
    std::vector<std::array<int, 4>> tetrahedra;

    std::vector<std::array<int, 2>> edges;
    std::map<std::array<int, 2>, int> edgeSingularTypes;
    // Ideally we do want this with a proper data structure with a STAR, then there's no write conflits
    // Make sure to always keep the edge (u, v) such that u < v in index value,
    std::map<std::array<int, 2>, std::set<int>> upperLink;
    std::map<std::array<int, 2>, std::set<int>> lowerLink;
    std::map<std::array<int, 2>, std::set<int>> upperStarTriangles;
    std::map<std::array<int, 2>, std::set<int>> lowerStarTriangles;

    std::vector<std::set<int>> triangles;
    std::unordered_map<std::set<int>, int, MyHash<std::set<int>>> triangleIndices;
    // Which triangles belong to the same tet as other triangles
    std::vector<std::vector<int>> tetIncidentTriangles;


    void computeBoundingBoxes();
    void computeDomainBoundingBox();
    void computeRangeBoundingBox();

    void sortVertices();
    void printMesh();
    void perturbRangeValues(const float &perturbationEpsilon);


    // We only add upperLink/lowerLink to data, the rest of data is unchagned
    void computeUpperLowerLinkAndStar();
    void computeSingularEdgeTypes();

    // From the tet soup, get the edges and triangles and tet-incidence of the triangles
    void computeCombinatorialStructure();

    // Give the edge (aIndex, bIndex), is the vertex vIndex from its link in the upper and lower link of the edge
    // We assume that aIndex < bIndex for consistent orientation.
    bool isUpperLinkEdgeVertex(int, int, int);
};
