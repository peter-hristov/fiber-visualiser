#pragma once

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

#include <qpoint.h>
#include <qvector.h>

#include "./FiberPoint.h"
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
    Arrangement &singularArrangement;
    ReebSpace &reebSpace;

    Data(TetMesh& tm, Arrangement& a, Arrangement& sa, ReebSpace& rs)
        : tetMesh(tm),
        arrangement(a),
        singularArrangement(sa),
        reebSpace(rs)
    {}

    std::string outputFibersFile = "./fibers.vtp";
};
