#include <QApplication>
#include <QDesktopWidget>
#include <QGraphicsScene>
#include <QLabel>
#include <QMainWindow>
#include <QMessageBox>
#include <QObject>
#include <QPainter>
#include <QVector>
#include <QtGui>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <qnamespace.h>
#include <utility>

#include "./TracerVisualiserWindow.h"
#include "./utility/Geometry.h"
#include "./utility/Utility.h"

#include "GlobalConfig.h"
#include "PlotWidget.h"
#include "src/Data.h"

using namespace std;

const bool DRAW_GRIDLINES = true;

PlotWidget::PlotWidget(QWidget* parent, Data* _data, string _interpolationType, tv9k::InputInformation input)
  : QWidget(parent)
  , data(_data)
  , resolution(input.scatterplotResolution)
  , interpolationType(_interpolationType)
{
    setMouseTracking(true);
    this->verticalLineNumbers = input.verticalLineNumbers;
    this->horizontalLineNumbers = input.horizontalLineNumbers;
}

void
PlotWidget::resetPoints(int mode)
{
    assert(mode == 0 || mode == 1);
    points[mode].clear();
    this->update();
}
void
PlotWidget::resetTriangles(int mode)
{
    assert(mode == 0 || mode == 1);
    triangles[mode].clear();
    this->update();
}

void
PlotWidget::addPoint(float x, float y, const int c, int mode)
{
    assert(mode == 0 || mode == 1);

    // Rescale [min, max] to [0, 100]
    float xx = (resolution / (data->uField->max - data->uField->min)) * (x - data->uField->min);
    float yy = (resolution / (data->vField->max - data->vField->min)) * (y - data->vField->min);

    // this->points.push_back(QPointF(10, 10));
    auto p = std::make_pair(QPointF((int)xx, (int)yy), c);

    this->points[mode].push_back(p);
    // this->update();
}

void
PlotWidget::addTriangle(float x1, float y1, float x2, float y2, float x3, float y3, int id, int mode)
{
    assert(mode == 0 || mode == 1);

    float xx1 = (resolution / (data->uField->max - data->uField->min)) * (x1 - data->uField->min);
    float yy1 = (resolution / (data->vField->max - data->vField->min)) * (y1 - data->vField->min);

    float xx2 = (resolution / (data->uField->max - data->uField->min)) * (x2 - data->uField->min);
    float yy2 = (resolution / (data->vField->max - data->vField->min)) * (y2 - data->vField->min);

    float xx3 = (resolution / (data->uField->max - data->uField->min)) * (x3 - data->uField->min);
    float yy3 = (resolution / (data->vField->max - data->vField->min)) * (y3 - data->vField->min);

    QVector<QPointF> triangle = { QPointF(xx1, yy1), QPointF(xx2, yy2), QPointF(xx3, yy3) };

    triangles[mode].push_back(make_pair(triangle, id));
    // this->update();
}

void
PlotWidget::keyPressEvent(QKeyEvent* event)
{
    // if (event->key() == Qt::Key_W)
    //{
    // cout << "I'M PRESSING!!!!!!";
    //}
}

void
PlotWidget::mouseDoubleClickEvent(QMouseEvent* event)
{}

void
PlotWidget::mouseMoveEvent(QMouseEvent* event)
{
#ifdef __APPLE__
    const QPointF clickPoint = { event->localPos().x() * 2, event->localPos().y() * 2 };
#else
    const QPointF clickPoint = event->localPos();
#endif
    int m = 0;

    auto &mousePoints = this->data->mousePoints[this->data->uField->name][this->data->vField->name];

    // Find out whether we're hovering on a vertex of the polygon and highlight it
    if (event->button() == Qt::NoButton) {
        // Have you clicked on a mouse point? Drag it.
        if (mousePoints.size() > 0) {
            for (int i = 0; i < mousePoints.size(); i++) {
                if (tv9k::geometry::getDistancePointPoint(clickPoint, mousePoints[i]) <
                    tv9k::geometry::getDistancePointPoint(clickPoint, mousePoints[m])) {
                    m = i;
                }
            }

            if (tv9k::geometry::getDistancePointPoint(clickPoint, mousePoints[m]) < sphereRadius * 4) {
                closePoint = m;
            } else {
                closePoint = -1;
            }
        } else {
            closePoint = -1;
        }
    }

    if (event->buttons() == Qt::LeftButton) {
        // Draggin a vertex
        if (MouseDragMode::Vertex == dragMode) {
            assert(movePoint >= 0 && movePoint < mousePoints.size());
            mousePoints[movePoint] = clickPoint;
        }

        // Draggin the polygon
        else if (MouseDragMode::Polygon == dragMode) {
            // Translate all the vertices of the polygon
            for (int i = 0; i < mousePoints.size(); i++) {
                mousePoints[i] -= (initialMovePoint - clickPoint);
            }
            this->initialMovePoint = clickPoint;
        }
    }

    this->update();
}

void
PlotWidget::mouseReleaseEvent(QMouseEvent* event)
{
    // If we've just dragged the polygon and it's finalised render fiber surface
    if (MouseDragMode::Polygon == dragMode && true == this->polygonLocked) {
        redoFS();
    }

    this->dragMode = MouseDragMode::Nothing;
}

void
PlotWidget::mousePressEvent(QMouseEvent* event)
{

#ifdef __APPLE__
    const QPointF clickPoint = { event->localPos().x() * 2, event->localPos().y() * 2 };
#else
    const QPointF clickPoint = event->localPos();
#endif

    auto &mousePoints = this->data->mousePoints[this->data->uField->name][this->data->vField->name];

    // 0 is outside the fscp, 1 in on a vertex of the fscp, 2 is inside the fscp
    enum ClickLocation
    {
        outside,
        vertex,
        inside
    } clickLocation = ClickLocation::outside;

    size_t closestEdge = -1;
    GLfloat distanceFromPolygon = -1;

    if (-1 != closePoint) {
        clickLocation = ClickLocation::vertex;
    } else {
        // If we at least have a triangle
        if (mousePoints.size() >= 3) {
            std::tie(distanceFromPolygon, closestEdge) =
              tv9k::geometry::getDistancePointPolygon(mousePoints, QPointF(clickPoint.x(), clickPoint.y()));

            if (distanceFromPolygon < 0) {
                clickLocation = ClickLocation::inside;
            }
        }
    }

    if (event->button() == Qt::RightButton) {
        if (clickLocation == ClickLocation::vertex) {
            mousePoints.erase(mousePoints.begin() + closePoint);
            this->polygonLocked = false;
            this->update();
        } else if (clickLocation == ClickLocation::inside) {
            this->polygonLocked = true;
            this->redoFS();
        }
    } else if (event->button() == Qt::LeftButton) {
        if (clickLocation == ClickLocation::vertex) {
            // Say we're dragging a vertex
            this->dragMode = MouseDragMode::Vertex;

            // Set the initial location
            this->movePoint = closePoint;

            this->polygonLocked = false;
            this->update();
        } else if (clickLocation == ClickLocation::inside) {
            // Say we're dragging the whole thing
            this->dragMode = MouseDragMode::Polygon;

            // Set the initial move point
            this->initialMovePoint = clickPoint;
        } else {
            mousePoints.insert(mousePoints.begin() + closestEdge + 1, clickPoint);
            this->polygonLocked = false;
            this->update();
        }
    }

    this->update();
}

void
PlotWidget::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.setWindow(QRect(0, 0, resolution, resolution));
    // p.setViewport(QRect(0, 0, resolution, resolution));

    p.save();
    p.setTransform(QTransform(1., 0., 0., -1., 0., resolution));
    p.setPen(Qt::gray);

    // Draw projetions
    if (data->projectionType == 1) {
        drawInteriorPoints(p);
    }

    else if (data->projectionType == 2) {
        drawIsosurfaceTriangles(p);
    }

    // Draw Fiber Surface Control Polygon
    drawAndRecomputeFS(p);

    p.restore();

    drawAxisLabels(p);
}

void
PlotWidget::drawInteriorPoints(QPainter& p)
{
    // Nothing to do if there's no isosurface
    if(this->data->isosurfaceMeshes[data->currentTimestep].triangles.size() == 0)
    {
        return;
    }

    // If we have a feature selected only display it
    if (SurfaceType::isosurface == this->data->selectedSurfaceType[this->data->currentTimestep])
    {
        for (int i = 0 ; i < data->xdim ; i++)
        {
            for (int j = 0 ; j < data->ydim ; j++)
            {
                for (int k = 0 ; k < data->zdim ; k++)
                {
                    // Get pont ID
                    const auto pointFeatureId = this->data->isosurfaceMeshes[data->currentTimestep].visited[i][j][k];

                    // Project point
                    if (
                            this->data->isoField->values[this->data->currentTimestep][i][j][k] > this->data->isoField->currentIsovalue && 
                            pointFeatureId == this->data->selectedID[this->data->currentTimestep])
                    {
                        float u = this->data->uField->values[this->data->currentTimestep][i][j][k];
                        float v = this->data->vField->values[this->data->currentTimestep][i][j][k];

                        float xx = (resolution / (data->uField->max - data->uField->min)) * (u - data->uField->min);
                        float yy = (resolution / (data->vField->max - data->vField->min)) * (v - data->vField->min);

                        // ID is the came as the color
                        p.setPen(QPen(Utility::getColorQt(pointFeatureId)));

                        int colorIndex = pointFeatureId;

                        //if (colorIndex == this->data->selectedID[this->data->currentTimestep]) {
                        //}

                        //const QPointF point(
                                //rescaleScalar(this->data->uField->min, this->data->uField->max, u),
                                //rescaleScalar(this->data->vField->min, this->data->vField->max, v));

                        p.setPen(QPen(Utility::getColorQt(colorIndex)));
                        p.drawPoint(QPointF{xx, yy});
                    }
                }
            }
        }
    }
    else
    {
        for (int i = 0 ; i < data->xdim ; i++)
        {
            for (int j = 0 ; j < data->ydim ; j++)
            {
                for (int k = 0 ; k < data->zdim ; k++)
                {
                    // Get pont ID
                    const auto pointFeatureId = this->data->isosurfaceMeshes[data->currentTimestep].visited[i][j][k];

                    // Project point
                    if (
                            this->data->isoField->values[this->data->currentTimestep][i][j][k] > this->data->isoField->currentIsovalue
                       )
                    {
                        float u = this->data->uField->values[this->data->currentTimestep][i][j][k];
                        float v = this->data->vField->values[this->data->currentTimestep][i][j][k];

                        // ID is the came as the color
                        p.setPen(QPen(Utility::getColorQt(pointFeatureId)));

                        int colorIndex = pointFeatureId;

                        //if (colorIndex == this->data->selectedID[this->data->currentTimestep]) {
                        //}

                        const QPointF point(
                                rescaleScalar(this->data->uField->min, this->data->uField->max, u),
                                rescaleScalar(this->data->vField->min, this->data->vField->max, v));

                        p.setPen(QPen(Utility::getColorQt(colorIndex)));
                        p.drawPoint(point);
                    }
                }
            }
        }

    }

    //for (int i = 0; i < this->points[0].size(); i++) {
        //int colorIndex = triangles[0][i].second;

        //if (colorIndex == this->data->selectedID[this->data->currentTimestep]) {
            //p.setPen(QPen(Utility::getColorQt(colorIndex)));
            //p.drawPoint(this->points[0][i].first);
        //}
    //}
}

GLfloat
PlotWidget::rescaleScalar(const GLfloat min, const GLfloat max, const GLfloat value)
{
    return (this->resolution / (max - min)) * (value - min);
}

void
PlotWidget::drawIsosurfaceTriangles(QPainter& p)
{
    // Draw transparent triangles
    for (const auto& triangle : this->data->isosurfaceMeshes[this->data->currentTimestep].triangles) {
        int colorIndex = triangle.triangleId;
        auto color = Utility::getColorQt(colorIndex);

        // If the haven't selected an isosurface object at all or we're an triangle from an object which is not selected
        if (SurfaceType::isosurface != this->data->selectedSurfaceType[this->data->currentTimestep] ||
            colorIndex != this->data->selectedID[this->data->currentTimestep]) {

            QVector<QPointF> projectedTriangle = {
                QPointF(
                  rescaleScalar(this->data->uField->min, this->data->uField->max, triangle.projectedVertices[0][0]),
                  rescaleScalar(this->data->vField->min, this->data->vField->max, triangle.projectedVertices[0][1])),
                QPointF(
                  rescaleScalar(this->data->uField->min, this->data->uField->max, triangle.projectedVertices[1][0]),
                  rescaleScalar(this->data->vField->min, this->data->vField->max, triangle.projectedVertices[1][1])),
                QPointF(
                  rescaleScalar(this->data->uField->min, this->data->uField->max, triangle.projectedVertices[2][0]),
                  rescaleScalar(this->data->vField->min, this->data->vField->max, triangle.projectedVertices[2][1]))
            };

            QPainterPath path;
            path.addPolygon(projectedTriangle);

            color.setAlpha(100);
            p.fillPath(path, QBrush(color));
        }
    }

    // Draw solid triangles
    for (const auto& triangle : this->data->isosurfaceMeshes[this->data->currentTimestep].triangles) {
        int colorIndex = triangle.triangleId;
        auto color = Utility::getColorQt(colorIndex);

        if (SurfaceType::isosurface == this->data->selectedSurfaceType[this->data->currentTimestep] &&
            colorIndex == this->data->selectedID[this->data->currentTimestep]) {

            QVector<QPointF> projectedTriangle = {
                QPointF(
                  rescaleScalar(this->data->uField->min, this->data->uField->max, triangle.projectedVertices[0][0]),
                  rescaleScalar(this->data->vField->min, this->data->vField->max, triangle.projectedVertices[0][1])),
                QPointF(
                  rescaleScalar(this->data->uField->min, this->data->uField->max, triangle.projectedVertices[1][0]),
                  rescaleScalar(this->data->vField->min, this->data->vField->max, triangle.projectedVertices[1][1])),
                QPointF(
                  rescaleScalar(this->data->uField->min, this->data->uField->max, triangle.projectedVertices[2][0]),
                  rescaleScalar(this->data->vField->min, this->data->vField->max, triangle.projectedVertices[2][1]))
            };

            QPainterPath path;
            path.addPolygon(projectedTriangle);

            color.setAlpha(200);
            p.fillPath(path, QBrush(color));
        }
    }
}

void
PlotWidget::drawAxisLabels(QPainter& p)
{
    auto penGrey = QPen(QColor(0, 0, 0, 50));
    penGrey.setWidth(0.5);

    auto penBlack = QPen(Qt::black);
    penBlack.setWidth(0.5);
    p.setPen(penBlack);

    float boxOffset = 15;

    // X Axis
    p.drawLine(boxOffset, resolution - boxOffset, resolution - boxOffset + 100, resolution - boxOffset);
    // Y Axis
    p.drawLine(boxOffset, resolution - boxOffset, boxOffset, boxOffset - 100);

    // Write out numbers
    QFont font = p.font();
    font.setPixelSize(6);
    // font.setWeight(20);
    p.setFont(font);

    //
    // Write out numbers on axis
    //
    float step = resolution / 15;

    float xMin = data->uField->min;
    float yMin = data->vField->min;
    float xMax = data->uField->max;
    float yMax = data->vField->max;
    float xRange = xMax - xMin;
    float yRange = yMax - yMin;
    float a = 10.0;
    int minSteps = 6;
    float xStepSize = pow(a, round(log(xRange) / log(a)) - 1);
    float yStepSize = pow(a, round(log(yRange) / log(a)) - 1);

    if (xRange / xStepSize < minSteps) {
        xStepSize /= 4.0;
    }
    if (yRange / yStepSize < minSteps) {
        yStepSize /= 4.0;
    }

    // ranges for labels and gridlines
    float xMinPlot = ceil(xMin / xStepSize) * xStepSize;
    float yMinPlot = ceil(yMin / yStepSize) * yStepSize;
    float xMaxPlot = floor(xMax / xStepSize) * xStepSize;
    float yMaxPlot = floor(yMax / yStepSize) * yStepSize;

    //
    // Draw Carthesian Grid and labels
    //

    // X labels and gridlines
    float xCurrent = xMinPlot;
    do {
        int i = int(resolution * (xCurrent - xMin) / xRange);
        if (i < boxOffset || i > (resolution - boxOffset)) {
            // skip this point
        } else {
            // grid line
            if (DRAW_GRIDLINES) {
                p.setPen(penGrey);
                p.drawLine(i, resolution - boxOffset - 20000, i, resolution - boxOffset);
            }
            // label
            p.setPen(penBlack);
            p.drawText(i - 5, resolution - boxOffset - 5, QString::number(xCurrent).mid(0, 6));
        }
        xCurrent += xStepSize;
        xCurrent = round(xCurrent / xStepSize) * xStepSize;
    } while (xCurrent < xMaxPlot);

    // Y labels and gridlines
    float yCurrent = yMinPlot;
    do {
        int i = int(resolution * (1.0 - (yCurrent - yMin) / yRange));
        if (i < boxOffset || i > (resolution - boxOffset)) {
            // skip this point
        } else {
            // grid line
            if (DRAW_GRIDLINES) {
                p.setPen(penGrey);
                p.drawLine(boxOffset, i, boxOffset + 2000, i);
            }
            // label
            p.setPen(penBlack);
            p.drawText(boxOffset + 5, i + 2, QString::number(yCurrent).mid(0, 6));
        }
        yCurrent += yStepSize;
        yCurrent = round(yCurrent / yStepSize) * yStepSize;
    } while (yCurrent < yMaxPlot);

    // pen2 = QPen(QColor(0, 100, 0, 50));
    // pen2.setWidth(0.5);
    // p.setPen(pen2);

    // float yCurrent = -1.0;
    // float yScaling = float(resolution)/(data->vField->max - data->vField->min);
    // do {
    // int i = int(yCurrent*(yScaling);
    // if (i > resolution) {
    // break;
    //}
    // yCurrent += yStepSize;
    //}
    // while (xCurrent < data->uField->max);

    // Y labels and gridlines
    // for (float i = resolution - boxOffset - step; i >= 0; i -= step) {
    //// label
    //// We're inverting the y-axis, hence the "1.0 -"
    // float ratio = 1.0 - (float)(i) / (float(resolution));
    // float a = ratio * (data->vField->max) + (1 - ratio) * (data->vField->min);
    // p.drawText(boxOffset + 5, i + 2, QString::number(a).mid(0, 6));
    //// grid line
    // p.drawLine(boxOffset, i, boxOffset + 2000, i);
    //}

    //
    // Draw custom lines
    //
    auto pen3 = QPen(QColor(0, 100, 0, 100));
    pen3.setWidth(1);
    p.setPen(pen3);

    // Vertical (constant X)
    for (const auto number : this->verticalLineNumbers) {
        if (data->uField->min < number && number < data->uField->max) {
            float ratio = (number - data->uField->min) / ((data->uField->max - data->uField->min));
            float plotPosition = ratio * (resolution);
            p.drawLine(plotPosition, resolution - boxOffset - 20000, plotPosition, resolution - boxOffset);
        }
    }

    // Horizontal (constant Y)
    for (const auto number : this->horizontalLineNumbers) {
        if (data->vField->min < number && number < data->vField->max) {
            // Invert the y axis
            float ratio = 1.0 - (number - data->vField->min) / ((data->vField->max - data->vField->min));
            float plotPosition = ratio * (resolution);
            p.drawLine(boxOffset, plotPosition, boxOffset + 2000, plotPosition);
        }
    }

    // Draw Labels
    p.setPen(penBlack);

    font.setPixelSize(7);
    p.setFont(font);

    // x label
    p.drawText(resolution / 2 - 30,
               resolution - boxOffset + 10,
               QString::fromStdString(data->uField->longName) + " (" + QString::fromStdString(data->uField->units) +
                 ")");

    p.translate(boxOffset - 5, resolution / 2 + 20);
    p.rotate(-90);

    // y label
    p.drawText(
      0, 0, QString::fromStdString(data->vField->longName) + " (" + QString::fromStdString(data->vField->units) + ")");
}

void
PlotWidget::drawAndRecomputeFS(QPainter& p)
{
    auto &mousePoints = this->data->mousePoints[this->data->uField->name][this->data->vField->name];

    // Only recompute the polygon if we have a new point
    polyPoints.erase(polyPoints.begin(), polyPoints.end());

    // Draw Fiber Surface Control Polygon (FSCP)
    for (int i = 0; i < mousePoints.size(); i++) {
        // Transform from the bottom left corner to be (0, 0) to the defautl QT
        // (upper left)
        QPointF a = p.combinedTransform().inverted().map(mousePoints[i]);
        // cout << " The point is was " << mousePoints[i].x() << " , " <<
        // mousePoints[i].y() << endl; cout << " The point is " << a.x() << " , " <<
        // a.y() << endl;
        polyPoints.push_back(a);
    }

    if (dragMode == 2) {
        p.setPen(Qt::green);
    } else {
        p.setPen(Qt::black);
    }
    p.drawConvexPolygon(QPolygonF(polyPoints));

    QPainterPath path;
    path.addPolygon(polyPoints);

    if (true == this->polygonLocked) {
        p.fillPath(path, QColor(100, 100, 100, 100));
    }

    // Draw vertices of the FSCP
    for (int i = 0; i < polyPoints.size(); i++) {
        p.setPen(Qt::black);

        // The current point is the one being dragged.
        if (closePoint == i || dragMode == 2) {
            p.setPen(Qt::green);
        }

        QPointF a = polyPoints[i];
        p.drawEllipse(QPointF(a.x(), a.y()), sphereRadius, sphereRadius);
    }
}

bool
PlotWidget::isPointInsidePolygon(std::vector<GLfloat> trianglePoint)
{
    // If we don't have at least a triangle
    if (this->polyPoints.size() <= 2) {
        return false;
    }

    const QPointF scaledProjectedPoint =
      tv9k::geometry::scaleProjectedPoint(data, this->resolution, trianglePoint[0], trianglePoint[1]);

    float signedDistance;
    std::tie(signedDistance, std::ignore) =
      tv9k::geometry::getDistancePointPolygon(this->polyPoints, scaledProjectedPoint);

    if (signedDistance > 0) {
        return false;
    } else {
        return true;
    }
}

void
PlotWidget::redoFS()
{
    // Update & repaint to update the fscp (rescaled version of the mouse points)
    this->repaint();
    this->update();
}

void
PlotWidget::resetDistanceField()
{}
