#pragma once

#include "src/CGALTypedefs.h"
#include <GL/gl.h>
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

class Data
{
  public:
    Data() {}

    // Min/max range coordinates
    GLfloat minF, maxF;
    GLfloat minG, maxG;

    // Min/max domain coordinates
    GLfloat minX, maxX;
    GLfloat minY, maxY;
    GLfloat minZ, maxZ;

    // Compute the min/max F, G and X, Y, Z coordinates
    void computeMinMaxRangeDomainCoordinates();

    // Which of the tets in the tetrahedra array contain a fiber
    std::vector<bool> tetsWithFibers;

    std::string longnameF, longnameG, units;

    QVector<QPointF> mousePoints;

    std::vector<std::vector<size_t>> tetrahedra;
    std::vector<std::vector<GLfloat>> vertexDomainCoordinates;

    std::vector<GLfloat> vertexCoordinatesF;
    std::vector<GLfloat> vertexCoordinatesG;


    std::vector<FaceFiberPoint> faceFibers;

    void computeTetExitPoints(const float, const float, const std::vector<float> = {1,1,1});


    void readData(std::string);



    // Reeb stuff relted content
    Arrangement_2 arr;

    // Ideally we do want this with a proper data structure with a STAR, then there's no write conflits
    // Make sure to always keep the edge (u, v) such that u < v in index value,
    std::map<std::pair<int, int>, std::set<int>> upperLink;
    std::map<std::pair<int, int>, std::set<int>> lowerLink;

    // Map from Int -> Point_2
    std::vector<Point_2> arrangementPoints;

    // Map from Point_2 -> Int
    std::map<Point_2, int> arrangementPointsIdices;
    std::map<Arrangement_2::Face_const_handle, int> arrangementFacesIdices;
    std::set<std::pair<std::set<int>, std::set<int>>> connectedTriangles;
    std::vector<int> arrangementFiberComponents;
    DisjointSet<std::pair<int, int>> reebSpace;
};




