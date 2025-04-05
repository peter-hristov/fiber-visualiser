#include <CGAL/number_utils.h>
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

#include "PlotWidget.h"
#include "src/Data.h"
#include "./Timer.h"


using namespace std;

const bool DRAW_GRIDLINES = true;


PlotWidget::PlotWidget(QWidget* parent, Data* _data, string _interpolationType)
    : QWidget(parent)
    , data(_data)
      , interpolationType(_interpolationType)
{
    setMouseTracking(true);
    this->setEnabled(true);








    //for (auto f = data->arr.faces_begin(); f != data->arr.faces_end(); ++f) 
    //{
        //if (f->is_unbounded()) 
        //{
            ////std::cout << "Unbounded face" << std::endl;
            //continue;
        //}

        //QVector<QPoint> points;

        //typename Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
        //typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
        //do {
            //typename Arrangement_2::Halfedge_const_handle e = curr;

            //// Get point from CGAL (and convert to double )
            //float u = CGAL::to_double(e->source()->point().x());
            //float v = CGAL::to_double(e->source()->point().y());

            //// Translate to the window frame
            //float u1 = (resolution / (data->maxF - data->minF)) * (u - data->minF);
            //float v1 = (resolution / (data->maxG - data->minG)) * (v - data->minG);


            //// Add to the polygon
            //points << QPoint(u1, v1);

            ////std::cout << "   (" << e->source()->point() << ")  -> " << "(" << e->target()->point() << ")" << std::endl;
        //} while (++curr != circ);



        //// We need to sort the colours by the sheetID, this way we get a nice colour blending
        //const int faceId = data->arrangementFacesIdices[f];
        //const vector<int> uniqueRoots = data->preimageGraphs[data->arrangementFacesIdices[f]].getUniqueRoots();

        //vector<int> sheetIds;

        //for (const int &fiberComponentID : uniqueRoots)
        //{
            //int vertexHIndex = data->vertexHtoIndex[{data->arrangementFacesIdices[f], fiberComponentID}];
            //int sheetID = data->reebSpace.findTriangle(vertexHIndex);

            //sheetIds.push_back(sheetID);
        //}

        //// Sort by the last element of each subvector
        //std::sort(sheetIds.begin(), sheetIds.end());

        //for (const int &sheetId : sheetIds)
        //{
            //this->arrangementPolygons << QPolygon(points);

            //int colourID = data->sheetToColour[sheetId];
            //vector<float> colorF = data->fiberColours[colourID % data->fiberColours.size()];
            //this->arrangementPolygonColours << QColor::fromRgbF(colorF[0], colorF[1], colorF[2], 0.392f);
        //}
    //}



    this->arrangementPolygons.resize(this->data->arrangementIndexToFace.size());
    this->arrangementPolygonColours.resize(this->data->arrangementIndexToFace.size());

    //
    // Used for lookup in the next loop
    //
    for (const auto &[vetexH, ufId] : data->reebSpace.data)
    {
        const auto &[faceID, fiberComponentID] = vetexH;

        Face_const_handle face = data->arrangementIndexToFace[faceID];

        if (face->is_unbounded()) 
        {
            continue;
        }

        int sheetID = data->reebSpace.findTriangle(vetexH);

        if (this->arrangementPolygons[faceID].size() == 0)
        {
            QVector<QPoint> points;

            typename Arrangement_2::Ccb_halfedge_const_circulator circ = face->outer_ccb();
            typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
            do {
                typename Arrangement_2::Halfedge_const_handle e = curr;

                // Get point from CGAL (and convert to double )
                float u = CGAL::to_double(e->source()->point().x());
                float v = CGAL::to_double(e->source()->point().y());

                // Translate to the window frame
                float u1 = (resolution / (data->maxF - data->minF)) * (u - data->minF);
                float v1 = (resolution / (data->maxG - data->minG)) * (v - data->minG);


                // Add to the polygon
                points << QPoint(u1, v1);

            } while (++curr != circ);

            this->arrangementPolygons[faceID] = QPolygon(points);
        }

        this->arrangementPolygonColours[faceID].push_back(sheetID);
    }

    // Sort so that the colours blend in properly
    for (int i = 0 ; i < this->arrangementPolygonColours.size() ; i++)
    {
        std::sort(this->arrangementPolygonColours[i].begin(), this->arrangementPolygonColours[i].end());
    }
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
    if (event->button() == Qt::NoButton) 
    {

        //
        // Old code for moving data points around
        //
        //if (this->data->vertexCoordinatesF.size() > 0) 
        //{

            //// Are we close to any of the vertices of the data?
            //for (int i = 0; i < this->data->vertexCoordinatesF.size(); i++) {
                //QPointF point;
                //point.setX((resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[i] - data->minF));
                //point.setY((resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[i] - data->minG));

                //QPointF bestPoint;
                //bestPoint.setX((resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[closestVertexPoint] - data->minF));
                //bestPoint.setY((resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[closestVertexPoint] - data->minG));

                //if (tv9k::geometry::getDistancePointPoint(mouseTransformedPoint, point) <
                        //tv9k::geometry::getDistancePointPoint(mouseTransformedPoint, bestPoint))
                //{
                    //// Set the closest vertex point
                    //closestVertexPoint = i;
                //}
            //}

            //QPointF bestPoint;
            //bestPoint.setX((resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[closestVertexPoint] - data->minF));
            //bestPoint.setY((resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[closestVertexPoint] - data->minG));

            //if (tv9k::geometry::getDistancePointPoint(mouseTransformedPoint, bestPoint) < sphereRadius * 4) {
                //this->closePointData = closestVertexPoint;
            //} else {
                //this->closePointData = -1;
            //}
        //}
        //else
        //{
            //this->closePointData = -1;
        //}


        // Have you clicked on a mouse point - the fiber point? Drag it.
        if (mousePoints.size() > 0 && this->closePointData == -1) 
        {
            // Find the closest mouse point
            for (int i = 0; i < mousePoints.size(); i++) 
            {
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

    if (event->buttons() == Qt::LeftButton) 
    {

        // Dragging a data point
        if (MouseDragMode::DataPoint == dragMode) 
        {
            //assert(movePoint >= 0 && movePoint < this->data->vertexCoordinatesF.size());

            //float fValue = this->data->minF + (mouseTransformedPoint.x() / resolution) * (this->data->maxF - this->data->minF);
            //float gValue = this->data->minG + (mouseTransformedPoint.y() / resolution) * (this->data->maxG - this->data->minG);

            //this->data->vertexCoordinatesF[movePoint] = fValue;
            //this->data->vertexCoordinatesG[movePoint] = gValue;
        }

        // Draggin a vertex
        if (MouseDragMode::Vertex == dragMode) {
            assert(movePoint >= 0 && movePoint < mousePoints.size());
            mousePoints[movePoint] = clickPoint;
            this->recomputeFiber = true;
        }

    }

    this->update();
}

    void
PlotWidget::mouseReleaseEvent(QMouseEvent* event)
{
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
            this->recomputeFiber = true;
            this->polygonLocked = false;
            this->update();
        }
    }

    this->update();
}

void PlotWidget::drawBackground(QPainter &p)
{
    // For each face
    for (int i = 0 ; i < this->arrangementPolygonColours.size() ; i++) 
    {
        // For each fiber component
        for (int j = 0 ; j < this->arrangementPolygonColours[i].size() ; j++) 
        {
            // Get the sheeID
            const int sheetId = this->arrangementPolygonColours[i][j];
            // Get the colour Id of the sheet ID
            int colourID = data->sheetToColour[sheetId];
            // The colour values
            vector<float> colorF = data->fiberColours[colourID % data->fiberColours.size()];

            // Draw the polygon with the appropriate colour
            p.setBrush(QColor::fromRgbF(colorF[0], colorF[1], colorF[2], 0.392f));
            p.setPen(Qt::NoPen);
            p.drawPolygon(this->arrangementPolygons[i]);
        }

    }


    //
    // Draw the polygons of the Reeb space
    //
    //for (int i = 0 ; i < this->arrangementPolygons.size() ; i++) 
    //{
        //// Set the random color for filling the polygon
        //p.setBrush(this->arrangementPolygonColours[i]);
        //p.setPen(Qt::NoPen);

        //// Draw the filled polygon with the random color
        //p.drawPolygon(this->arrangementPolygons[i]);

        ////qDebug() << "New polygon --- ";
        ////for (const QPoint& point : this->arrangementPolygons[i]) {
            ////qDebug() << "(" << point.x() << ", " << point.y() << ")";
        ////}
    //}

    //
    // Draw the Jacobi set
    //
    for (const auto &[edge, type] : data->jacobiType)
    {
        if (type != 1)
        {
            if (type == 0)
            {
                p.setPen(QPen(Qt::black, 0.1, Qt::DashLine));
            }
            else
            {
                p.setPen(QPen(Qt::black, 0.1));
            }

            float x1 = (resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[edge.first] - data->minF);
            float y1 = (resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[edge.first] - data->minG);

            float x2 = (resolution / (data->maxF - data->minF)) * (this->data->vertexCoordinatesF[edge.second] - data->minF);
            float y2 = (resolution / (data->maxG - data->minG)) * (this->data->vertexCoordinatesG[edge.second] - data->minG);

            p.drawLine(x1, y1, x2, y2);

        }

    }

    drawAxisLabels(p);
}

void PlotWidget::generateStaticCache()
{
    staticCache = std::make_unique<QPixmap>(resolution, resolution);
    staticCache->fill(Qt::white);

    QPainter p(staticCache.get());
    this->drawBackground(p);
}

void PlotWidget::resizeEvent(QResizeEvent* event)
{
    // If the widget size has changed, regenerate the static cache
    //if (!staticCache || staticCache->size() != size()) 
    //{
        //generateStaticCache();
    //}

    // You can call the base resizeEvent if needed, but it's not strictly required.
    // QWidget::resizeEvent(event); // Uncomment if needed
}

void PlotWidget::paintEvent(QPaintEvent*)
{


    QPainter p(this);
    p.setWindow(QRect(0, 0, resolution, resolution));
    // p.setViewport(QRect(0, 0, resolution, resolution));

    p.save();
    p.setTransform(QTransform(1., 0., 0., -1., 0., resolution));
    p.setPen(Qt::gray);

    this->painterCombinedTransform = p.combinedTransform();

    if (!staticCache || staticCache->size() != size()) 
    {
        generateStaticCache();
    }
    p.drawPixmap(0, 0, *(this->staticCache));

    //this->drawBackground(p);


    // Draw Fiber Surface Control Polygon
    drawAndRecomputeFS(p);

    p.restore();

}


    GLfloat
PlotWidget::rescaleScalar(const GLfloat min, const GLfloat max, const GLfloat value)
{
    return (this->resolution / (max - min)) * (value - min);
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

    // Leave a trail of fibers or not
    if (dynamic_cast<TracerVisualiserWindow*>(this->parent())->tracerVisualiserWidget->clearFibers == true)
    {
        this->data->faceFibers.clear();
    }

    // If we have selected a point in the scatterplot, compute the fiber surface
    if (this->recomputeFiber == true && polyPoints.size() == 1)
    {
        this->recomputeFiber = false;

        // Reside to the original range data dimensions
        float u = this->data->minF + (polyPoints[0].x() / resolution) * (this->data->maxF - this->data->minF);
        float v = this->data->minG + (polyPoints[0].y() / resolution) * (this->data->maxG - this->data->minG);

        const int currentFiberColour = 0;
        // Compute all tet exit points

        qDebug() << "Computing fiber ...";

        //Timer::start();
        //this->data->computeTetExitPoints(u, v, data->fiberColours[currentFiberColour]);
        //Timer::stop("Computed fiber                         :");

        //Timer::start();
        ////this->data->computeTetExitPoints(u, v, data->fiberColours[currentFiberColour]);
        //Timer::stop("Computed fiber old                     :");
        //qDebug() << "Next ...";

        //Timer::start();
        //this->data->computeTetExitPointsNew(u, v, data->fiberColours[currentFiberColour]);
        //Timer::stop("Computed fiber new                     :");
        //qDebug() << "Next ...";

        //Timer::start();
        this->data->computeTetExitPointsNewNew(u, v, data->fiberColours[currentFiberColour]);
        //Timer::stop("Computed fiber new new                 :");
        //qDebug() << "......";

        // Display fibers
        const auto& visualiserWidget = dynamic_cast<TracerVisualiserWindow*>(this->parent())->tracerVisualiserWidget;
        visualiserWidget->generateDisplayList();
        visualiserWidget->update();
    }
}
