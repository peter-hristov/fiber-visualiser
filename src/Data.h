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

    GLfloat minX, maxX;
    GLfloat minY, maxY;
    GLfloat minZ, maxZ;

    std::string longnameF, longnameG, units;

    QVector<QPointF> mousePoints;

    std::vector<std::vector<size_t>> tetrahedra;
    std::vector<std::vector<GLfloat>> vertexDomainCoordinates;

    std::vector<GLfloat> vertexCoordinatesF;
    std::vector<GLfloat> vertexCoordinatesG;

    // The location of a fiber on a single tet
    struct FaceFiber{
        float alpha;
        float betta;
        std::vector<size_t> vertices;
    };

    // The location 
    struct MeshTriangle{
        std::vector<float> vertixA;
        std::vector<float> vertixB;
        std::vector<float> vertixC;
    };


    std::vector<FaceFiber> faceFibers;
    std::vector<MeshTriangle> meshTriangles;


    // Original dimensions of the data before cropping and downsampling - used for
    // reference
    int originalXdim, originalYdim, originalZdim, originalTdim;

    // New dimensions after cropping and downsampling - used in the application
    int xdim, ydim, zdim, tdim;

    void readNcData(tv9k::InputInformation);

    size_t trippleToIndex(const size_t i, const size_t j, const size_t k);
    void addTetsForCube(const size_t i, const size_t j, const size_t k);
};
