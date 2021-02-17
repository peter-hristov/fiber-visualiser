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

#include <QImage>

#include "./ScatterPlot.h"
#include "./utility/ScalarField.h"
#include "./utility/SurfaceMesh.h"
#include "./utility/Utility.h"

#include "GlobalConfig.h"

// @TODO This shouldn't be defined here.
enum SurfaceType {none, isosurface, fibersurface, combinedSurface};

class Data
{
  public:
    Data() {}
    GLfloat min, max;
    std::string longName, units;


    QVector<QPointF> mousePoints;


    // Original dimensions of the data before cropping and downsampling - used for
    // reference
    int originalXdim, originalYdim, originalZdim, originalTdim;

    // New dimensions after cropping and downsampling - used in the application
    int xdim, ydim, zdim, tdim;



    void readNcData(tv9k::InputInformation);

    // These should be somewhere else

    // @TODO Move these to a separate namespace for input data
    bool render3DLabels = false;
    bool drawLines = false;
    bool continuousScatterPlot = false;
    bool flatNormals = false;
    bool precomputeMergeTrees = false;
    bool cacheJoinTree = false;
    int projectionType = 2;
    bool dynamicPolygon = 2;

};
