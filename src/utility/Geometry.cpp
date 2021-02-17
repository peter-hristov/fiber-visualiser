#include <QtMath>
#include <assert.h>

#include "Geometry.h"

using namespace std;

// GLfloat Geometry::cubeVertices[8][3] = {{0,  0,  0 }, { 0,  0,  1 }, { 0,  1,
// 0 }, { 0,  1,  1 }, { 1,  0,  0 }, { 1,  0,  1 }, { 1,  1,  0 }, { 1,  1,  1 }
// };

float
tv9k::geometry::dotProduct(QPointF a, QPointF b)
{
    return a.x() * b.x() + a.y() * b.y();
}

QPointF
tv9k::geometry::getNormal(QPointF v)
{
    QPointF normal(-v.y(), v.x());
    float length = qSqrt(normal.x() * normal.x() + normal.y() * normal.y());
    return QPointF(normal.x() / length, normal.y() / length);
}

float
tv9k::geometry::getDistancePointPoint(QPointF p1, QPointF p2)
{
    return qSqrt(dotProduct(p2 - p1, p2 - p1));
}

float
tv9k::geometry::getDistancePointLine(QPointF p1, QPointF p2, QPointF x)
{
    return dotProduct(x - p1, getNormal(p1 - p2));
}

float
tv9k::geometry::getCross(QPointF u, QPointF v)
{
    return u.x() * v.y() - u.y() * v.x();
}

std::pair<float, size_t>
tv9k::geometry::getDistancePointPolygon(QVector<QPointF> points, QPointF x)
{
    int n = points.size();

    // If we have less than 2 points, don't compute
    if (points.size() < 2) {
        return { -1, 0 };
        // If we have 2 points, it's a half plane
    } else if (points.size() == 2) {
        return { -1 * getDistancePointLine(points[0], points[1], x), 0 };
    }

    // For a polygon make it circular, by inserting the first point last
    points.push_back(QPointF(points[0].x(), points[0].y()));
    n++;

    float distance[n - 1];

    // Determine whether the polygon is clockwise -1.0 or anticlockwise 1.0
    float cross = getCross(points[1] - points[0], points[2] - points[0]) < 0 ? -1.0 : 1.0;

    // @TODO This is a horrible implementation, rewrite all of it
    if (cross < 0) {
        // std::cout << "Printing : " << std::endl;
        for (int i = n - 1; i > 0; i--) {
            float angle1 = dotProduct(points[i - 1] - points[i], x - points[i]);
            float angle2 = dotProduct(points[i] - points[i - 1], x - points[i - 1]);
            float sign = getDistancePointLine(points[i], points[i - 1], x) < 0 ? -1.0 : 1.0;

            // @TODO Optimise this using Pythogoras
            if (angle1 < 0) {
                distance[i - 1] = sign * qSqrt(dotProduct(points[i] - x, points[i] - x));
                // std::cout << "Distance is 1 from (" << x.x() << ", " << x.y() << ")
                // to (" << points[i].x() << ", " << points[i].y() <<  ") is " <<
                // distance[i] << std::endl;
            } else if (angle2 < 0) {
                distance[i - 1] = sign * qSqrt(dotProduct(points[i - 1] - x, points[i - 1] - x));
                // std::cout << "Distance is 2 from (" << x.x() << ", " << x.y() << ")
                // to (" << points[i + 1].x() << ", " << points[i + 1].y() <<  ") is " <<
                // distance[i] << std::endl;
            } else {
                distance[i - 1] = getDistancePointLine(points[i], points[i - 1], x);
                // std::cout << "Distance is 3 from is " << distance[i] << std::endl;
            }
        }
    } else {
        for (int i = 0; i < n - 1; i++) {
            float angle1 = dotProduct(points[i + 1] - points[i], x - points[i]);
            float angle2 = dotProduct(points[i] - points[i + 1], x - points[i + 1]);
            float sign = getDistancePointLine(points[i], points[i + 1], x) < 0 ? -1.0 : 1.0;

            // @TODO Optimise this using Pythogoras
            if (angle1 < 0) {
                distance[i] = sign * qSqrt(dotProduct(points[i] - x, points[i] - x));
                // std::cout << "Distance is 1 from (" << x.x() << ", " << x.y() << ")
                // to (" << points[i].x() << ", " << points[i].y() <<  ") is " <<
                // distance[i] << std::endl;
            } else if (angle2 < 0) {
                distance[i] = sign * qSqrt(dotProduct(points[i + 1] - x, points[i + 1] - x));
                // std::cout << "Distance is 2 from (" << x.x() << ", " << x.y() << ")
                // to (" << points[i + 1].x() << ", " << points[i + 1].y() <<  ") is " <<
                // distance[i] << std::endl;
            } else {
                distance[i] = getDistancePointLine(points[i], points[i + 1], x);
                // std::cout << "Distance is 3 from is " << distance[i] << std::endl;
            }
        }
    }

    // std::cout << std::endl << std::endl;

    bool foundPositive = false;

    float maxNegative = -1.0 * std::numeric_limits<float>::max();
    float minPositive = std::numeric_limits<float>::max();

    int maxNegativeIndex = -1;
    int minPositiveIndex = -1;

    for (int i = 0; i < n - 1; i++) {

        // cout << "The distance is " << distance[i] << endl;

        if (distance[i] < 0) {
            if (maxNegative < distance[i]) {
                maxNegative = distance[i];
                maxNegativeIndex = i;
            }
            // maxNegative = std::max(maxNegative, distance[i]);
        }

        // At least one outside
        else if (distance[i] > 0) {
            if (minPositive > distance[i]) {
                minPositive = distance[i];
                minPositiveIndex = i;
            }
            // minPositive = std::min(minPositive, distance[i]);
            foundPositive = true;
        }

        // On the line
        else {
            return { 0.0, i };
        }
    }

    if (foundPositive == true) {
        return { minPositive, minPositiveIndex };
    }

    return { maxNegative, maxNegativeIndex };
}

std::vector<GLfloat>
tv9k::geometry::computeCentralDifferencingNormal(std::vector<GLfloat> currentVertex,
                                                 int u,
                                                 int v,
                                                 GLfloat l,
                                                 int xdim,
                                                 int ydim,
                                                 int zdim,
                                                 const std::vector<std::vector<std::vector<float>>>& vals,
                                                 bool flipNormal)
{
    int x1 = currentVertex[0] + cubeVertices[u][0];
    int y1 = currentVertex[1] + cubeVertices[u][1];
    int z1 = currentVertex[2] + cubeVertices[u][2];

    int x2 = currentVertex[0] + cubeVertices[v][0];
    int y2 = currentVertex[1] + cubeVertices[v][1];
    int z2 = currentVertex[2] + cubeVertices[v][2];

    // cout << "Here are the coordinates (" << x1 << " " << y1 << " " << z1 << ")
    // , (" << x2 << " " << y2 << " " << z2 << ")" << endl;

    float df_dx1 = -1.0;
    float df_dy1 = -1.0;
    float df_dz1 = -1.0;

    float df_dx2 = -1.0;
    float df_dy2 = -1.0;
    float df_dz2 = -1.0;

    if (0 < x1 && x1 < xdim - 1) {
        df_dx1 = (vals[x1 + 1][y1][z1] - vals[x1 - 1][y1][z1]) / 2;
    } else if (x1 == 0) {
        df_dx1 = vals[x1 + 1][y1][z1] - vals[x1][y1][z1];
    } else if (x1 == xdim - 1) {
        df_dx1 = vals[x1][y1][z1] - vals[x1 - 1][y1][z1];
    } else {
        assert(false);
    }

    if (0 < y1 && y1 < ydim - 1) {
        df_dy1 = (vals[x1][y1 + 1][z1] - vals[x1][y1 - 1][z1]) / 2;
    } else if (0 == y1) {
        df_dy1 = vals[x1][y1 + 1][z1] - vals[x1][y1][z1];
    } else if (y1 = ydim - 1) {
        df_dy1 = vals[x1][y1][z1] - vals[x1][y1 - 1][z1];
    } else {
        assert(false);
    }

    if (0 < z1 && z1 < zdim - 1) {
        df_dz1 = (vals[x1][y1][z1 + 1] - vals[x1][y1][z1 - 1]) / 2;
    } else if (0 == z1) {
        df_dz1 = vals[x1][y1][z1 + 1] - vals[x1][y1][z1];
    } else if (z1 == zdim - 1) {
        df_dz1 = vals[x1][y1][z1] - vals[x1][y1][z1 - 1];
    } else {
        assert(false);
    }

    if (0 < x2 && x2 < xdim - 1) {
        df_dx2 = (vals[x2 + 1][y2][z2] - vals[x2 - 1][y2][z2]) / 2;
    } else if (0 == x2) {
        df_dx2 = vals[x2 + 1][y2][z2] - vals[x2][y2][z2];
    } else if (x2 == xdim - 1) {
        df_dx2 = vals[x2][y2][z2] - vals[x2 - 1][y2][z2];
    } else {
        assert(false);
    }

    if (0 < y2 && y2 < ydim - 1) {
        df_dy2 = (vals[x2][y2 + 1][z2] - vals[x2][y2 - 1][z2]) / 2;
    } else if (0 == y2) {
        df_dy2 = vals[x2][y2 + 1][z2] - vals[x2][y2][z2];
    } else if (y2 == ydim - 1) {
        df_dy2 = vals[x2][y2][z2] - vals[x2][y2 - 1][z2];
    } else {
        assert(false);
    }

    if (0 < z2 && z2 < (zdim - 1)) {
        df_dz2 = (vals[x2][y2][z2 + 1] - vals[x2][y2][z2 - 1]) / 2;
    } else if (0 == z2) {
        df_dz2 = vals[x2][y2][z2 + 1] - vals[x2][y2][z2];
    } else if (z2 == zdim - 1) {
        df_dz2 = vals[x2][y2][z2] - vals[x2][y2][z2 - 1];
    } else {
        assert(false);
    }

    float normalX = l * df_dx1 + (1 - l) * df_dx2;
    float normalY = l * df_dy1 + (1 - l) * df_dy2;
    float normalZ = l * df_dz1 + (1 - l) * df_dz2;

    float normalLength = sqrt(normalX * normalX + normalY * normalY + normalZ * normalZ);

    vector<GLfloat> normal = { float(1.0 * normalX / normalLength),
                               float(1.0 * normalY / normalLength),
                               float(1.0 * normalZ / normalLength) };

    if (flipNormal) {
        normal[0] *= -1.0;
        normal[1] *= -1.0;
        normal[2] *= -1.0;
    }

    return normal;
}

GLfloat
tv9k::geometry::bilinearInterpolation(const float x, const float y, const QVector<QVector<float>> distanceField)
{
    float leftX = floor(x);
    float rightX = ceil(x);

    float bottomY = floor(y);
    float topY = ceil(y);

    float values[2][2] = { { distanceField[static_cast<int>(leftX)][static_cast<int>(bottomY)],
                             distanceField[static_cast<int>(rightX)][static_cast<int>(bottomY)] },
                           { distanceField[static_cast<int>(leftX)][static_cast<int>(bottomY)],
                             distanceField[static_cast<int>(rightX)][static_cast<int>(bottomY)] } };

    float xxRatio = (x - leftX) / (rightX - leftX);
    float yyRatio = (y - bottomY) / (bottomY - bottomY);

    float top = xxRatio * values[0][0] + (1 - xxRatio) * values[0][1];
    float bottom = xxRatio * values[1][0] + (1 - xxRatio) * values[1][1];

    float middle = yyRatio * top + (1 - yyRatio) * bottom;

    return middle;
}

QPointF
tv9k::geometry::scaleProjectedPoint(const Data* data, const float resolution, GLfloat pointU, const GLfloat pointV)
{
    // Rescale to squre
    float rescaledU = 10.0;
      //(static_cast<int>(resolution) / (data->uField->max - data->uField->min)) * (pointU - data->uField->min);
    float rescaledV = 10.0;
      //(static_cast<int>(resolution) / (data->vField->max - data->vField->min)) * (pointV - data->vField->min);

    return QPointF(rescaledU, rescaledV);
}

std::vector<std::vector<std::vector<GLfloat>>>
tv9k::geometry::compute3DDistanceField(Data* data,
                                       const QVector<QPointF> polyPoints,
                                       const int resolution,
                                       const std::string interpolationType,
                                       const int mult)
{
    QVector<QVector<float>> distanceField;
    std::vector<std::vector<std::vector<GLfloat>>> signedDistanceField(
      data->xdim, std::vector<std::vector<GLfloat>>(data->ydim, std::vector<GLfloat>(data->zdim)));
    return signedDistanceField;
}
