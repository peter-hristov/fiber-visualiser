#pragma once

#include <QApplication>
#include <QDesktopWidget>
#include <QLabel>
#include <QMainWindow>
#include <QMessageBox>
#include <QObject>
#include <QPainter>
#include <QVector>
#include <QtGui>
#include <QtMath>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <limits>
#include <map>
#include <qvector.h>
#include <utility>

#include "Data.h"
#include <cmath>

class PlotWidget : public QWidget
{
    Q_OBJECT

  signals:
    void redrawFiberSurface();

  public:
    // Index in the mousePoints array
    int closePoint = -1;
    int closePointData = -1;
    // Index in the mousePoints array
    int movePoint = -1;
    // Whether the fiber surface for the polygon has been rendered
    bool polygonLocked = false;
    // This is the initial point when tranclating the polygon
    QPointF initialMovePoint;
    
    //int 

    QTransform painterCombinedTransform;


    Data* data;
    float varianceScale = 0;
    float resolution = 300;
    // This is the radius of the sphere around a vertex in the fscp
    float sphereRadius = 3;
    std::string interpolationType = "";

    QPointF drawLinesPoint;

    GLfloat rescaleScalar(const GLfloat, const GLfloat, const GLfloat);

    float mult = 1.0;

    // Points on the FSCP polygon
    QVector<QPointF> polyPoints;

    enum MouseDragMode
    {
        Nothing,
        Vertex,
        Polygon,
        DataPoint
    } dragMode = MouseDragMode::Nothing;

    QVector<std::pair<QPointF, int>> points[2];
    QVector<std::pair<QVector<QPointF>, int>> triangles[2];

    // Deprecated
    // map<int, QImage> imageMap;
    // map<int, QPainter> painterMap;

    // map<int, QImage> imagePointMap;
    // map<int, QPainter> painterPointMap;

    // QVector<pair<QVector<QPointF>>> squares;

    QVector<QVector<float>> distanceField;

    //QVector<QPointF> mousePoints;
    QVector<QVector<int>> histogram;

    PlotWidget()
      : QWidget()
    {
        setFocusPolicy(Qt::StrongFocus);
    }

    PlotWidget(QWidget* parent, Data* _data, std::string _interpolationType);

    std::vector<float> verticalLineNumbers;
    std::vector<float> horizontalLineNumbers;

    void paintEvent(QPaintEvent* event);
    void resetDistanceField();
    void redoFS();

    void drawCube(float vertices[8][2]);
    void drawLine(float x1, float y1, float x2, float y2, int id, int);
    void addTriangle(float x1, float y1, float x2, float y2, float x3, float y3, int, int);
    void addPoint(float x, float y, const int c, int mode);

    void resetPoints(int);
    void resetTriangles(int);

    // void resetLines(int);

    void drawIsosurfaceTriangles(QPainter&);

    bool isPointInsidePolygon(std::vector<GLfloat>);

  protected:
    void mouseMoveEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mousePressEvent(QMouseEvent* event);
    void mouseDoubleClickEvent(QMouseEvent* event);
    void keyPressEvent(QKeyEvent* event);

    void computeFSCPDistanceField();
    float getFSCPDistance(float x, float y);
    void drawInteriorPointsImages(QPainter& p);
    void drawInteriorPoints(QPainter& p);
    void drawAxisLabels(QPainter& p);
    void drawAndRecomputeFS(QPainter& p);

    void compareImages();
    void generateTriangleImages();
};
