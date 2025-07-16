#pragma once

#include "src/Arrangement.h"
#include "src/CGALTypedefs.h"
#include "src/ReebSpace.h"
#include "src/TetMesh.h"
#include <GL/gl.h>
#include <unordered_map>
#ifdef __APPLE__
#include <OpenGL/glu.h>
#include <OpenGL/gl.h>
#else
#include <GL/glu.h>
#include <GL/gl.h>
#endif

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "./FaceFiber.h"
#include "./DisjointSet.h"

#include <qpoint.h>
#include <qvector.h>

#include "./CGALTypedefs.h"
#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./ReebSpace.h"

class Data
{
  public:

    // Ideally, I would like a move constructor, but there's some issues in computing the search structure for the arrangement, something is not moved.
    // For not, just pass by reference

    TetMesh &tetMesh;
    Arrangement &arrangement;
    ReebSpace &reebSpace;

    Data(TetMesh& tm, Arrangement& a, ReebSpace& rs)
        : tetMesh(tm),
        arrangement(a),
        reebSpace(rs)
    {}




    //
    // Other
    //
    std::vector<float> currentFiberPoint;
    void printSheetHistogram();

    // Fibers class?
    std::vector<bool> tetsWithFibers;

    std::string fibersFile = "./fibers.vtp";
    std::vector<FaceFiberPoint> faceFibers;
    void saveFibers();

    void generatefFaceFibersForSheet(const int sheetId, const int numberOfFiberPoints);
    void generatefFaceFibersForSheets(const int sheetOutputCount, const int numberOfFiberPoints, const std::string);

    void computeTetExitPoints(const float, const float, const std::vector<float> = {1,1,1});
    void computeTetExitPointsNew(const float, const float, const std::vector<float> = {1,1,1});
    void computeTetExitPointsNewNew(const float, const float, const bool, const int reebSheetIdOnly = -1, const std::vector<float> = {1,1,1});

    QVector<QPointF> mousePoints;

    // Colour map for Reeb space sheets and fibers
    const std::vector<std::vector<float>> fiberColours = {
        {1.0f, 0.0f, 0.0f, 0.392f},    // Vivid Red
        {0.0f, 1.0f, 0.0f, 0.392f},    // Bright Green
        {0.0f, 0.0f, 1.0f, 0.392f},    // Pure Blue
        {1.0f, 1.0f, 0.0f, 0.392f},    // Bright Yellow
        {0.0f, 1.0f, 1.0f, 0.392f},    // Cyan
        {1.0f, 0.0f, 1.0f, 0.392f},    // Magenta
        {0.58f, 0.0f, 1.0f, 0.392f},   // Deep Purple
        {0.0f, 0.45f, 0.7f, 0.392f},   // Ocean Blue
        {1.0f, 0.5f, 0.0f, 0.392f},    // Orange
        {0.0f, 0.6f, 0.5f, 0.392f},    // Teal
        {1.0f, 0.84f, 0.0f, 0.392f},   // Gold
        {0.85f, 0.4f, 0.55f, 0.392f},  // Mauve
        {0.4f, 0.8f, 1.0f, 0.392f},    // Sky Blue
        {0.2f, 0.8f, 0.2f, 0.392f},    // Leaf Green
        {0.9f, 0.3f, 0.3f, 0.392f},    // Coral Red
        {0.6f, 0.6f, 0.0f, 0.392f}     // Olive
    };
};
