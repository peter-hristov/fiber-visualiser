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

    // Bounding box min/max for the domain and range coordinates of all vertices
    float minF, maxF;
    float minG, maxG;
    float minX, maxX;
    float minY, maxY;
    float minZ, maxZ;

    // (Optional) indicative names for the two scalar fields and their units
    std::string longnameF, longnameG, units;

    // Tets, edges, triangles, etc
    std::vector<std::array<size_t, 4>> tetrahedra;

    std::vector<std::array<int, 2>> edges;
    std::map<std::array<int, 2>, int> edgeSingularTypes;


    std::vector<std::set<int>> indexToTriangle;
    std::unordered_map<std::set<int>, int, MyHash<std::set<int>>> triangleToIndex;


    // What does this do again?
    std::vector<std::vector<int>> adjacentTrianglesIndex;

    std::vector<float> vertexCoordinatesF;
    std::vector<float> vertexCoordinatesG;
    std::vector<std::vector<float>> vertexDomainCoordinates;

    // Ideally we do want this with a proper data structure with a STAR, then there's no write conflits
    // Make sure to always keep the edge (u, v) such that u < v in index value,
    std::map<std::array<int, 2>, std::set<int>> upperLink;
    std::map<std::array<int, 2>, std::set<int>> lowerLink;

    std::map<std::array<int, 2>, std::set<int>> upperStarTriangles;
    std::map<std::array<int, 2>, std::set<int>> lowerStarTriangles;

    // Tell us which triangles (as sets of IDs) are connected (part of a tetrahedron)
    std::set<std::pair<std::set<int>, std::set<int>>> connectedTriangles;


    void computeBoundingBoxes();
    void computeDomainBoundingBox();
    void computeRangeBoundingBox();

    void sortVertices();
    void printMesh();
    void perturbRangeValues(const float &perturbationEpsilon);


    // We only add upperLink/lowerLink to data, the rest of data is unchagned
    void computeUpperLowerLinkVertices();
    void computeSingularEdgeTypes();
    // We only add upperLink/lowerLink to data, the rest of data is unchagned
    void computeTriangleAdjacency();

    // Give the edge (aIndex, bIndex), is the vertex vIndex from its link in the upper and lower link of the edge
    // We assume that aIndex < bIndex for consistent orientation.
    bool isUpperLinkEdgeVertex(int, int, int);
};
