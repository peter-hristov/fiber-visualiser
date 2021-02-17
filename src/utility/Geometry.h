#pragma once

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "../Data.h"
#include <QPointF>
#include <QVector>
#include <vector>

namespace tv9k {
namespace geometry {
const GLfloat cubeVertices[8][3] = { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, 0 }, { 0, 1, 1 },
                                     { 1, 0, 0 }, { 1, 0, 1 }, { 1, 1, 0 }, { 1, 1, 1 } };

std::vector<GLfloat>
computeCentralDifferencingNormal(std::vector<GLfloat> currentVertex,
                                 int,
                                 int,
                                 GLfloat,
                                 int,
                                 int,
                                 int,
                                 const std::vector<std::vector<std::vector<float>>>&,
                                 bool);

float dotProduct(QPointF, QPointF);
QPointF getNormal(QPointF);
float getDistancePointPoint(QPointF, QPointF);
float getDistancePointLine(QPointF, QPointF, QPointF);
float getCross(QPointF, QPointF);
std::pair<float, size_t> getDistancePointPolygon(QVector<QPointF>, QPointF);

GLfloat
bilinearInterpolation(const float x, const float y, const QVector<QVector<float>> distanceField);

QPointF
scaleProjectedPoint(const Data* data, const float resolution, GLfloat pointU, const GLfloat pointV);

std::vector<std::vector<std::vector<GLfloat>>>
compute3DDistanceField(Data* data,
                       const QVector<QPointF> polyPoints,
                       const int resolution,
                       const std::string interpolationType,
                       const int mult);

}
}
