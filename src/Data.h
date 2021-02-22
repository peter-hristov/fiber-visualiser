#pragma once

#include <GL/gl.h>
#ifdef __APPLE__
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <qpoint.h>
#include <qvector.h>
#include "GlobalConfig.h"

class Data
{
  public:
    Data() {}

    GLfloat minF, maxF;
    GLfloat minG, maxG;

    std::string longName, units;

    QVector<QPointF> mousePoints;

    std::vector<std::vector<size_t>> tetrahedra;
    std::vector<std::vector<GLfloat>> vertexDomainCoordinates;
    std::vector<std::vector<GLfloat>> vertexRangeCoordinates;

    struct FaceFiber{
        float alpha;
        float betta;
        std::vector<size_t> vertices;
    };


    std::vector<FaceFiber> faceFibers;


    // Original dimensions of the data before cropping and downsampling - used for
    // reference
    int originalXdim, originalYdim, originalZdim, originalTdim;

    // New dimensions after cropping and downsampling - used in the application
    int xdim, ydim, zdim, tdim;

    void readNcData(tv9k::InputInformation);

    size_t trippleToIndex(const size_t i, const size_t j, const size_t k);
    void addTetsForCube(const size_t i, const size_t j, const size_t k);
};
