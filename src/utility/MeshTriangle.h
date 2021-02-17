#pragma once

#include <vector>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include <QPointF>
#include <QVector>

namespace tv9k {
namespace utility {
class MeshTriangle
{
  public:
    GLint triangleId{ -1 };
    std::vector<GLfloat> color;
    std::vector<std::vector<GLfloat>> vertices;
    std::vector<std::vector<GLfloat>> normals;
    std::vector<std::vector<GLfloat>> projectedVertices;

    QVector<QPointF> projectedTriangle;

    inline void clear()
    {
        this->vertices.clear();
        this->normals.clear();
        this->color.clear();
        this->projectedVertices.clear();
        triangleId = -1;
    }
};
}
}
