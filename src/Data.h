#pragma once

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

#include <qpoint.h>
#include <qvector.h>

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

    void computeTetExitPoints(const float, const float);


    void readData(std::string);
};
