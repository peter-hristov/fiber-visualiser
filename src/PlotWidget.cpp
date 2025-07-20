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


#include "./PlotWidget.h"
#include "./Data.h"
#include "./Timer.h"
#include "./CGALTypedefs.h"
#include "./utility/Geometry.h"
#include "./TracerVisualiserWindow.h"

using namespace std;

const bool DRAW_GRIDLINES = true;


PlotWidget::PlotWidget(QWidget* parent, Data &_data, string _interpolationType)
    : QWidget(parent)
    , data(_data)
      , interpolationType(_interpolationType)
{
    setMouseTracking(true);
    setEnabled(true);

    paddedMinF = data.tetMesh.minF - paddingScalingFactor * (data.tetMesh.maxF - data.tetMesh.minF);
    paddedMaxF = data.tetMesh.maxF + paddingScalingFactor * (data.tetMesh.maxF - data.tetMesh.minF);
    paddedMinG = data.tetMesh.minG - paddingScalingFactor * (data.tetMesh.maxG - data.tetMesh.minG);
    paddedMaxG = data.tetMesh.maxG + paddingScalingFactor * (data.tetMesh.maxG - data.tetMesh.minG);

    // Set up polygons for drawing
    //this->arrangementPolygons.resize(this->data.arrangementIndexToFace.size());
    //for (auto f = data.arr.faces_begin(); f != data.arr.faces_end(); ++f) 
    //{
        //if (f->is_unbounded()) 
        //{
            //continue;
        //}

        //QVector<QPoint> points;

        //typename Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
        //typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
        //do {
            //typename Arrangement_2::Halfedge_const_handle e = curr;

            //// Get point from CGAL (and convert to double )
            //const float u = CGAL::to_double(e->source()->point().x());
            //const float v = CGAL::to_double(e->source()->point().y());

            //// Translate to the window frame
            //const float u1 = (resolution / (data.maxF - data.minF)) * (u - data.minF);
            //const float v1 = (resolution / (data.maxG - data.minG)) * (v - data.minG);


            //// Add to the polygon
            //points << QPoint(u1, v1);

            ////std::cout << "   (" << e->source()->point() << ")  -> " << "(" << e->target()->point() << ")" << std::endl;
        //} while (++curr != circ);

        //const int faceId = data.arrangementFacesIdices[f];
        //this->arrangementPolygons[faceId] = QPolygon(points);
    //}
}

    void
PlotWidget::keyPressEvent(QKeyEvent* event)
{
    //if (event->key() == Qt::Key_I) {
    //this->data.mousePoints[0].setY(this->data.mousePoints[0].y() + 1);
    //this->update();
    //}
    //if (event->key() == Qt::Key_J) {
    //this->data.mousePoints[0].setX(this->data.mousePoints[0].x() - 1);
    //this->update();
    //}
    //if (event->key() == Qt::Key_K) {
    //this->data.mousePoints[0].setY(this->data.mousePoints[0].y() - 1);
    //this->update();
    //}
    //if (event->key() == Qt::Key_L) {
    //this->data.mousePoints[0].setX(this->data.mousePoints[0].x() + 1);
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

    auto &mousePoints = this->data.mousePoints;

    // Find out whether we're hovering on a vertex of the polygon and highlight it
    if (event->button() == Qt::NoButton) 
    {

        //
        // Old code for moving data points around
        //
        //if (this->data.vertexCoordinatesF.size() > 0) 
        //{

            //// Are we close to any of the vertices of the data?
            //for (int i = 0; i < this->data.vertexCoordinatesF.size(); i++) {
                //QPointF point;
                //point.setX((resolution / (data.maxF - data.minF)) * (this->data.vertexCoordinatesF[i] - data.minF));
                //point.setY((resolution / (data.maxG - data.minG)) * (this->data.vertexCoordinatesG[i] - data.minG));

                //QPointF bestPoint;
                //bestPoint.setX((resolution / (data.maxF - data.minF)) * (this->data.vertexCoordinatesF[closestVertexPoint] - data.minF));
                //bestPoint.setY((resolution / (data.maxG - data.minG)) * (this->data.vertexCoordinatesG[closestVertexPoint] - data.minG));

                //if (tv9k::geometry::getDistancePointPoint(mouseTransformedPoint, point) <
                        //tv9k::geometry::getDistancePointPoint(mouseTransformedPoint, bestPoint))
                //{
                    //// Set the closest vertex point
                    //closestVertexPoint = i;
                //}
            //}

            //QPointF bestPoint;
            //bestPoint.setX((resolution / (data.maxF - data.minF)) * (this->data.vertexCoordinatesF[closestVertexPoint] - data.minF));
            //bestPoint.setY((resolution / (data.maxG - data.minG)) * (this->data.vertexCoordinatesG[closestVertexPoint] - data.minG));

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

            if (tv9k::geometry::getDistancePointPoint(clickPoint, mousePoints[m]) <= sphereRadius) {
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
            //assert(movePoint >= 0 && movePoint < this->data.vertexCoordinatesF.size());

            //float fValue = this->data.minF + (mouseTransformedPoint.x() / resolution) * (this->data.maxF - this->data.minF);
            //float gValue = this->data.minG + (mouseTransformedPoint.y() / resolution) * (this->data.maxG - this->data.minG);

            //this->data.vertexCoordinatesF[movePoint] = fValue;
            //this->data.vertexCoordinatesG[movePoint] = gValue;
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

    auto &mousePoints = this->data.mousePoints;

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

    for (const auto &[sheetId, polygon] : data.reebSpace.sheetPolygon)
    {
        QVector<QPointF> points;

        // If the sheet is incomplete, the polygon will not be corret, just draww all the faces manually
        if (data.reebSpace.incompleteSheets.contains(sheetId))
        {

            // Loop through all faces to see which ones are in the sheet
            for (auto f = data.arrangement.arr.faces_begin(); f != data.arrangement.arr.faces_end(); ++f) 
            {
                const int currentFaceID = data.arrangement.arrangementFacesIdices[f];

                // For each fiber component in the face, see if one of those is in our sheet
                for (const auto &[triangleId, fiberComponentId] : this->data.reebSpace.fiberSeeds[currentFaceID])
                {
                    const int componentSheetId = data.reebSpace.reebSpace.findTriangle({currentFaceID, fiberComponentId});

                    // Now we can add the polygon
                    if (componentSheetId == sheetId)
                    {
                        typename Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
                        typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
                        do {
                            typename Arrangement_2::Halfedge_const_handle e = curr;

                            // Get point from CGAL (and convert to double )
                            const float u = CGAL::to_double(e->source()->point().x());
                            const float v = CGAL::to_double(e->source()->point().y());

                            // Add to the polygon
                            points << rescalePoint(u, v);

                            //std::cout << "   (" << e->source()->point() << ")  -> " << "(" << e->target()->point() << ")" << std::endl;
                        } while (++curr != circ);
                    }
                }
            }
        }

        // If the sheet is no incomplete, the polygon is valid, draw the directly
        else
        {
            for (const CartesianPoint &point : polygon) 
            {
                // Get point from CGAL (and convert to double )
                const float u = point.x();
                const float v = point.y();

                // Add to the polygon
                points << rescalePoint(u, v);
            }
        }

        QPolygonF qPolygon(points);

        const int colourID = data.reebSpace.sheetToColour[sheetId];
        const vector<float> colorF = data.fiberColours[colourID % data.fiberColours.size()];

        p.setBrush(QColor::fromRgbF(colorF[0], colorF[1], colorF[2], 0.392f));
        p.setPen(Qt::NoPen);
        p.drawPolygon(qPolygon);
    }

    //for (int i = 0 ; i < this->data.fiberSeeds.size() ; i++) 
    //{
        //const int currentFaceID = i;

        //// For each fiber component
        //for (int j = 0 ; j < this->data.fiberSeeds[i].size() ; j++) 
        //{
            ////const auto &[]
            //const auto &[triangleId, fiberComponentId] = this->data.fiberSeeds[i][j];
            //const int sheetId = this->data.reebSpace.findTriangle({currentFaceID, fiberComponentId});
            //const int colourID = data.sheetToColour[sheetId];
            //const vector<float> colorF = data.fiberColours[colourID % data.fiberColours.size()];

            //p.setBrush(QColor::fromRgbF(colorF[0], colorF[1], colorF[2], 0.392f));
            //p.setPen(Qt::NoPen);
            //p.drawPolygon(this->arrangementPolygons[i]);
        //}
    //}

    // We assume that the fiber seeds per face are sorted by their sheetId
    //for (int i = 0 ; i < this->data.fiberSeeds.size() ; i++) 
    //{
        //const int currentFaceID = i;

        //// For each fiber component
        //for (int j = 0 ; j < this->data.fiberSeeds[i].size() ; j++) 
        //{
            ////const auto &[]
            //const auto &[triangleId, fiberComponentId] = this->data.fiberSeeds[i][j];
            //const int sheetId = this->data.reebSpace.findTriangle({currentFaceID, fiberComponentId});
            //const int colourID = data.sheetToColour[sheetId];
            //const vector<float> colorF = data.fiberColours[colourID % data.fiberColours.size()];

            //p.setBrush(QColor::fromRgbF(colorF[0], colorF[1], colorF[2], 0.392f));
            //p.setPen(Qt::NoPen);
            //p.drawPolygon(this->arrangementPolygons[i]);
        //}
    //}


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
    //for (const auto &[edge, type] : data.reebSpace.jacobiType)
    //{
        //if (type != 1)
        //{
            //if (type == 0)
            //{
                ////p.setPen(QPen(Qt::black, 1, Qt::DashLine));
                //p.setPen(QPen(Qt::black, 0.2));
            //}
            //else
            //{
                //p.setPen(QPen(Qt::black, 0.2));
            //}

            //float x1 = (resolution / (data.tetMesh.maxF - data.tetMesh.minF)) * (this->data.tetMesh.vertexCoordinatesF[edge.first] - data.tetMesh.minF);
            //float y1 = (resolution / (data.tetMesh.maxG - data.tetMesh.minG)) * (this->data.tetMesh.vertexCoordinatesG[edge.first] - data.tetMesh.minG);

            //float x2 = (resolution / (data.tetMesh.maxF - data.tetMesh.minF)) * (this->data.tetMesh.vertexCoordinatesF[edge.second] - data.tetMesh.minF);
            //float y2 = (resolution / (data.tetMesh.maxG - data.tetMesh.minG)) * (this->data.tetMesh.vertexCoordinatesG[edge.second] - data.tetMesh.minG);

            //p.setRenderHint(QPainter::Antialiasing, true);
            //p.drawLine(x1, y1, x2, y2);
        //}
    //}

    // Draw all edges
    for (const auto &[edge, type] : data.tetMesh.edgeSingularTypes)
    {
        const float u1 = this->data.tetMesh.vertexCoordinatesF[edge[0]];
        const float v1 = this->data.tetMesh.vertexCoordinatesG[edge[0]];

        const float u2 = this->data.tetMesh.vertexCoordinatesF[edge[1]];
        const float v2 = this->data.tetMesh.vertexCoordinatesG[edge[1]];

        if (type == 0)
        {
            p.setPen(QPen(Qt::black, 4.2, Qt::DashLine));
        }
        else if (type == 1)
        {
            p.setPen(QPen(Qt::black, 3.2, Qt::SolidLine));
        }
        else
        {
            p.setPen(QPen(Qt::black, 10.2, Qt::SolidLine));
        }

        p.setRenderHint(QPainter::Antialiasing, true);
        p.drawLine(rescalePoint(u1, v1), rescalePoint(u2, v2));
    }

    // Draw all the vertex coordinates
    for(size_t i = 0 ; i <  this->data.tetMesh.vertexCoordinatesF.size() ; i++)
    {
        float u = this->data.tetMesh.vertexCoordinatesF[i];
        float v = this->data.tetMesh.vertexCoordinatesG[i];

        p.setPen(QPen(Qt::black, 6, Qt::SolidLine));
        p.setBrush(Qt::white);           // Fill color
        p.drawEllipse(rescalePoint(u, v), 20, 20);

        QFont font = p.font();
        font.setPointSize(70);
        p.setFont(font);
        p.drawText(rescalePoint(u, v), QString::number(i));
    }
}

void PlotWidget::generateStaticCache()
{
    if (!staticCache) 
    {
        staticCache = std::make_unique<QPixmap>(resolution, resolution);
        staticCache->fill(Qt::white);
        qDebug() << "Redrawing REEB SPACE ...";
        QPainter p(staticCache.get());
        this->drawBackground(p);
    }
}

void PlotWidget::resizeEvent(QResizeEvent* event)
{
    generateStaticCache();
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

    if (!staticCache) 
    {
        generateStaticCache();
    }
    p.drawPixmap(0, 0, *(this->staticCache));

    //this->drawBackground(p);


    // Draw Fiber Surface Control Polygon
    drawAndRecomputeFS(p);

    p.restore();
    drawAxisLabels(p);
}


    GLfloat
PlotWidget::rescaleScalar(const GLfloat min, const GLfloat max, const GLfloat value)
{
    return (this->resolution / (max - min)) * (value - min);
}

QPointF PlotWidget::rescalePoint(const float &u, const GLfloat &v)
{
    const float rescaledU = (resolution / (paddedMaxF - paddedMinF)) * (u - paddedMinF);
    const float rescaledV = (resolution / (paddedMaxG - paddedMinG)) * (v - paddedMinG);
    return QPointF(rescaledU, rescaledV);
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

    float xMin = paddedMinF;
    float yMin = paddedMinG;
    float xMax = paddedMaxF;
    float yMax = paddedMaxG;
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
        if (paddedMinF < number && number < paddedMaxF) {
            float ratio = (number - paddedMinF) / ((paddedMaxF - paddedMinF));
            float plotPosition = ratio * (resolution);
            p.drawLine(plotPosition, resolution - boxOffset - 20000, plotPosition, resolution - boxOffset);
        }
    }

    // Horizontal (constant Y)
    for (const auto number : this->horizontalLineNumbers) {
        if (paddedMinG < number && number < paddedMaxG) {
            // Invert the y axis
            float ratio = 1.0 - (number - paddedMinG) / ((paddedMaxG - paddedMinG));
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
            QString::fromStdString(data.tetMesh.longnameF) + " (" + QString::fromStdString(data.tetMesh.units) +
            ")");

    p.translate(boxOffset - 5, resolution / 2 + 20);
    p.rotate(-90);

    // y label
    p.drawText(
            0, 0, QString::fromStdString(data.tetMesh.longnameG) + " (" + QString::fromStdString(data.tetMesh.units) + ")");
}

    void
PlotWidget::drawAndRecomputeFS(QPainter& p)
{

    auto &mousePoints = this->data.mousePoints;

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
    //if (dynamic_cast<TracerVisualiserWindow*>(this->parent())->tracerVisualiserWidget->clearFibers == true)
    //{
        //this->data.faceFibers.clear();
    //}

    // If we have selected a point in the scatterplot, compute the fiber surface
    if (this->recomputeFiber == true && polyPoints.size() == 1)
    {
        this->recomputeFiber = false;

        // Reside to the original range data dimensions
        float u = this->paddedMinF + (polyPoints[0].x() / resolution) * (this->paddedMaxF - this->paddedMinF);
        float v = this->paddedMinG + (polyPoints[0].y() / resolution) * (this->paddedMaxG - this->paddedMinG);

        const int currentFiberColour = 0;
        // Compute all tet exit points

        qDebug() << "Computing fiber ...";

        //Timer::start();
        //this->data.computeTetExitPoints(u, v, data.fiberColours[currentFiberColour]);
        //Timer::stop("Computed fiber                         :");

        //Timer::start();
        ////this->data.computeTetExitPoints(u, v, data.fiberColours[currentFiberColour]);
        //Timer::stop("Computed fiber old                     :");
        //qDebug() << "Next ...";

        //Timer::start();
        //this->data.computeTetExitPointsNew(u, v, data.fiberColours[currentFiberColour]);
        //Timer::stop("Computed fiber new                     :");
        //qDebug() << "Next ...";

        //Timer::start();
        const bool clearFibers = dynamic_cast<TracerVisualiserWindow*>(this->parent())->tracerVisualiserWidget->clearFibers;
        this->data. computeTetExitPointsNewNew(u, v, clearFibers, -1, data.fiberColours[currentFiberColour]);
        //Timer::stop("Computed fiber new new                 :");
        //qDebug() << "......";

        // Display fibers
        const auto& visualiserWidget = dynamic_cast<TracerVisualiserWindow*>(this->parent())->tracerVisualiserWidget;
        visualiserWidget->generateDisplayList();
        visualiserWidget->update();
    }
}
