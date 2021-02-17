#pragma once

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <vector>

#include "../Data.h"
#include "./MeshTriangle.h"
#include "./SurfaceMesh.h"

namespace tv9k {
namespace utility {
namespace MarchingCubes {
void
computeTriangles(GLfloat,
                 const std::vector<std::vector<std::vector<GLfloat>>>&,
                 const std::vector<std::vector<std::vector<GLfloat>>>&,
                 const std::vector<std::vector<std::vector<GLfloat>>>&,
                 int,
                 tv9k::utility::SurfaceMesh&,
                 // std::vector<std::tuple<float, float, int>>&,
                 const int,
                 Data*);

void
prepareCube(int,
            int,
            int,
            GLfloat,
            const std::vector<std::vector<std::vector<GLfloat>>>&,
            const std::vector<std::vector<std::vector<GLfloat>>>&,
            const std::vector<std::vector<std::vector<GLfloat>>>&,
            int,
            tv9k::utility::SurfaceMesh&,
            // std::vector<std::tuple<float, float, int>>&,
            const int,
            Data*);

void
processCube(int,
            const GLfloat[8],
            const GLfloat[8],
            const GLfloat[8],
            GLfloat,
            int,
            int,
            int,
            int,
            tv9k::utility::SurfaceMesh&,
            // std::vector<std::tuple<float, float, int>>&,
            const int,
            const std::vector<std::vector<std::vector<GLfloat>>>&,
            Data*);

};
}
}
