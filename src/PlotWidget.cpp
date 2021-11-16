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
#include <qcolor.h>
#include <qnamespace.h>
#include <qpoint.h>
#include <qtransform.h>
#include <utility>

#include "./TracerVisualiserWindow.h"
#include "./utility/Geometry.h"

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
    this->setEnabled(true);
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
}

void
PlotWidget::addTriangle(float x1, float y1, float x2, float y2, float x3, float y3, int id, int mode)
{
}

void
PlotWidget::keyPressEvent(QKeyEvent* event)
{
    //if (event->key() == Qt::Key_I) {
        //this->data->mousePoints[0].setY(this->data->mousePoints[0].y() + 1);
        //this->update();
    //}
    //if (event->key() == Qt::Key_J) {
        //this->data->mousePoints[0].setX(this->data->mousePoints[0].x() - 1);
        //this->update();
    //}
    //if (event->key() == Qt::Key_K) {
        //this->data->mousePoints[0].setY(this->data->mousePoints[0].y() - 1);
        //this->update();
    //}
    //if (event->key() == Qt::Key_L) {
        //this->data->mousePoints[0].setX(this->data->mousePoints[0].x() + 1);
        //this->update();
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

    // Transform the [res, res] grid
    auto mouseTransformedPoint = painterCombinedTransform.inverted().map(clickPoint);

    int m = 0;
    int closestVertexPoint = 0;

    auto &mousePoints = this->data->mousePoints;

    // Find out whether we're hovering on a vertex of the polygon and highlight it
    if (event->button() == Qt::NoButton) {

        if (this->data->vertexCoordinatesF.size() > 0) {

            for (int i = 0; i < this->data->vertexCoordinatesF.size(); i++) {
                QPointF point;
                point.setX((resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[i] - data->minF));
                point.setY((resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[i] - data->minG));

                QPointF bestPoint;
                bestPoint.setX((resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[closestVertexPoint] - data->minF));
                bestPoint.setY((resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[closestVertexPoint] - data->minG));

                if (tv9k::geometry::getDistancePointPoint(mouseTransformedPoint, point) <
                    tv9k::geometry::getDistancePointPoint(mouseTransformedPoint, bestPoint))
                {
                    closestVertexPoint = i;
                }
            }

            QPointF bestPoint;
            bestPoint.setX((resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[closestVertexPoint] - data->minF));
            bestPoint.setY((resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[closestVertexPoint] - data->minG));

            if (tv9k::geometry::getDistancePointPoint(mouseTransformedPoint, bestPoint) < sphereRadius * 4) {
                this->closePointData = closestVertexPoint;
            } else {
                this->closePointData = -1;
            }
        }
        else
        {
            this->closePointData = -1;
        }

        // Have you clicked on a mouse point? Drag it.
        if (mousePoints.size() > 0 && this->closePointData == -1) {
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

        // Dragging a data point
        if (MouseDragMode::DataPoint == dragMode) {
            assert(movePoint >= 0 && movePoint < this->data->vertexCoordinatesF.size());

            float fValue = this->data->minF + (mouseTransformedPoint.x() / resolution) * (this->data->maxF - this->data->minF);
            float gValue = this->data->minG + (mouseTransformedPoint.y() / resolution) * (this->data->maxG - this->data->minG);

            this->data->vertexCoordinatesF[movePoint] = fValue;
            this->data->vertexCoordinatesG[movePoint] = gValue;


            //QPointF bestPoint;
            //bestPoint.setX((resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[closestVertexPoint] - data->minF));
            //bestPoint.setY((resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[closestVertexPoint] - data->minG));





            //mousePoints[movePoint] = clickPoint;
        }

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

    auto &mousePoints = this->data->mousePoints;

    // 0 is outside the fscp, 1 in on a vertex of the fscp, 2 is inside the fscp
    enum ClickLocation
    {
        outside,
        vertex,
        inside,
        dataPoint
    } clickLocation = ClickLocation::outside;

    size_t closestEdge = -1;
    GLfloat distanceFromPolygon = -1;

    if (-1 != this->closePointData) {
        clickLocation = ClickLocation::dataPoint;
    }
    else if (-1 != this->closePoint)
    {
        clickLocation = ClickLocation::vertex;
    }
    else {
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
        if (clickLocation == ClickLocation::dataPoint)
        {
            this->dragMode = MouseDragMode::DataPoint;

            // Set the initial location
            this->movePoint = this->closePointData;
        }
        else if (clickLocation == ClickLocation::vertex) {
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
            mousePoints.clear();
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

    this->painterCombinedTransform = p.combinedTransform();

    // Draw Fiber Surface Control Polygon
    drawAndRecomputeFS(p);


    p.restore();

    drawAxisLabels(p);
}

void
PlotWidget::drawInteriorPoints(QPainter& p)
{
}

GLfloat
PlotWidget::rescaleScalar(const GLfloat min, const GLfloat max, const GLfloat value)
{
    return (this->resolution / (max - min)) * (value - min);
}

void
PlotWidget::drawIsosurfaceTriangles(QPainter& p)
{
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

    float xMin = data->minF;
    float yMin = data->minG;
    float xMax = data->maxF;
    float yMax = data->maxG;
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

    //
    // Draw custom lines
    //
    auto pen3 = QPen(QColor(0, 100, 0, 100));
    pen3.setWidth(1);
    p.setPen(pen3);

    // Vertical (constant X)
    for (const auto number : this->verticalLineNumbers) {
        if (data->minF < number && number < data->maxF) {
            float ratio = (number - data->minF) / ((data->maxF - data->minF));
            float plotPosition = ratio * (resolution);
            p.drawLine(plotPosition, resolution - boxOffset - 20000, plotPosition, resolution - boxOffset);
        }
    }

    // Horizontal (constant Y)
    for (const auto number : this->horizontalLineNumbers) {
        if (data->minG < number && number < data->maxG) {
            // Invert the y axis
            float ratio = 1.0 - (number - data->minG) / ((data->maxG - data->minG));
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
               QString::fromStdString(data->longnameF) + " (" + QString::fromStdString(data->units) +
                 ")");

    p.translate(boxOffset - 5, resolution / 2 + 20);
    p.rotate(-90);

    // y label
    p.drawText(
      0, 0, QString::fromStdString(data->longnameG) + " (" + QString::fromStdString(data->units) + ")");
}

void
PlotWidget::drawAndRecomputeFS(QPainter& p)
{
    auto &mousePoints = this->data->mousePoints;

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
        auto penGrey = QPen(QColor(0, 0, 0, 250));
        p.setPen(penGrey);

        // The current point is the one being dragged.
        if (closePoint == i || dragMode == 2) {
            p.setPen(Qt::green);
        }


        QPointF a = polyPoints[i];

        p.drawEllipse(QPointF(a.x(), a.y()), sphereRadius, sphereRadius);

        
        penGrey.setWidthF(1.0);
        p.setPen(penGrey);

        

        p.drawLine(a.x(), a.y() - resolution, a.x(), a.y() + resolution);
        p.drawLine(a.x() - resolution, a.y(), a.x() + resolution, a.y());
    }



    for(size_t i = 0 ; i <  this->data->vertexCoordinatesF.size() ; i++)
    {

        float x1 = (resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[i] - data->minF);
        float y1 = (resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[i] - data->minG);

        if (this->closePointData == i)
        {
            p.setPen(Qt::blue);
        }
        else
        {
            if (i == 13)
            {
                //p.setPen(Qt::red);
                p.setPen(Qt::black);
            }
            else
            {
                p.setPen(Qt::black);
            }
        }

        p.drawEllipse(QPointF(x1, y1), 3, 3);
        //p.setTransform(QTransform(1., 0., 0., -1., 0., resolution));
        p.drawText(x1, y1, QString::number(i));
    }

    auto penGrey = QPen(QColor(0, 0, 0, 225));
    penGrey.setWidthF(0.8);
    p.setPen(penGrey);


    if (dynamic_cast<TracerVisualiserWindow*>(this->parent())->tracerVisualiserWidget->clearFibers == true)
    {
        this->data->faceFibers.clear();
    }

    this->data->tetsWithFibers = vector<bool>(this->data->tetrahedra.size(), false);

    // Draw all triangles from the tets
    for(size_t tetId = 0 ; tetId < this->data->tetrahedra.size(); tetId++)
    {
        const auto tet = this->data->tetrahedra[tetId];

        // For every edge in the tet
        for(int i = 0 ; i < 4 ; i++)
        {
            for(int j = i + 1 ; j < 4 ; j++)
            {
                float x1 = (resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[tet[i]] - data->minF);
                float y1 = (resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[tet[i]] - data->minG);

                float x2 = (resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[tet[j]] - data->minF);
                float y2 = (resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[tet[j]] - data->minG);

                p.drawLine(x1, y1, x2, y2);
            }
        }

        if (mousePoints.size() == 0)
            continue;


        // For every triangle in every tet, get the fiber in it
        for(int i = 0 ; i < 4 ; i++)
        {
            for(int j = i + 1 ; j < 4 ; j++)
            {
                for(int k = j + 1 ; k < 4 ; k++)
                {
                    float x1 = this->data->vertexCoordinatesF[tet[i]];
                    float y1 = this->data->vertexCoordinatesG[tet[i]];

                    float x2 = this->data->vertexCoordinatesF[tet[j]];
                    float y2 = this->data->vertexCoordinatesG[tet[j]];

                    float x3 = this->data->vertexCoordinatesF[tet[k]];
                    float y3 = this->data->vertexCoordinatesG[tet[k]];

                    float x = this->data->minF + (polyPoints[0].x() / resolution) * (this->data->maxF - this->data->minF);
                    float y = this->data->minG + (polyPoints[0].y() / resolution) * (this->data->maxG - this->data->minG);

                    float det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

                    float alpha = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / det;
                    float betta = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / det;
                    float gamma = 1 - alpha - betta;

                    if (alpha >= 0 && betta >= 0 & gamma >= 0 && alpha <= 1 && betta <= 1 && gamma <= 1)
                    {
                        //printf("In triangle %d, %d, %d in tet %d.\n", tet[i], tet[j], tet[k], tetId);
                        //printf("In triangle (%f, %f) | (%f, %f) | (%f, %f) comparing with point (%f, %f) and alpha = %f, betta = %f, gamma = %f.\n", x1, y1, x2, y2, x3, y3, x, y, alpha, betta, gamma);
                        Data::FaceFiber fb;
                        fb.alpha = alpha;
                        fb.betta = betta;
                        fb.vertices = {tet[i], tet[j], tet[k]};

                        float pointBackX = alpha * x1 + betta * x2 + gamma * x3;
                        float pointBackY = alpha * y1 + betta * y2 + gamma * y3;

                        //printf ("Projecting back we get (%f, %f).\n", pointBackX, pointBackY);

                        const auto& visualiserWidget = dynamic_cast<TracerVisualiserWindow*>(this->parent())->tracerVisualiserWidget;
                        fb.colour = visualiserWidget->fiberColour;
                        this->data->faceFibers.push_back(fb);
                        this->data->tetsWithFibers[tetId] = true;

    //struct FaceFiber{
        //float alpha;
        //float betta;
        //std::vector<size_t> vertices;
    //};

                    }

                }
            }
        }
    }

    const float isovalue = 0.0;


    // Generate the fs meshes for all timesteps and the display list for the current timestep
    const auto& visualiserWidget = dynamic_cast<TracerVisualiserWindow*>(this->parent())->tracerVisualiserWidget;
    visualiserWidget->generateDisplayList();

    float iso = this->data->minF + (polyPoints[0].x() / resolution) * (this->data->maxF - this->data->minF);
    visualiserWidget->generateDisplayListTriangles(iso, this->data->vertexCoordinatesF, 1);

    float iso2 = this->data->minG + (polyPoints[0].y() / resolution) * (this->data->maxG - this->data->minG);
    visualiserWidget->generateDisplayListTriangles(iso2, this->data->vertexCoordinatesG, 2);
    visualiserWidget->update();

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
