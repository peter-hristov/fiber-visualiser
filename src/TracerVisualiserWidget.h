#pragma once

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <QGLWidget>
#include <map>

#include "./ArcBall/Ball.h"
#include "./Data.h"
#include "./HistogramWidget.hpp"
#include "./PlotWidget.h"
#include "./utility/DisjointSet.hpp"
#include "./utility/MeshTriangle.h"
#include "./utility/SurfaceMesh.h"

class TracerVisualiserWidget : public QGLWidget
{
  public:
    TracerVisualiserWidget(QWidget*, QWidget*, QWidget*, Data*);
    GLfloat scale = 0;
    // GLfloat isovalue = -1;

    GLfloat isovalueMult = -1.0;

    int displayListIndex = 0;
    int displayListIndexF = 0;
    int displayListIndexC = 0;
    // bool drawLines = false;

    // int maxVisited = 0;
    int maxVisitedF = 0;

    bool shouldIntersect = true;
    bool shouldHideIsosurface = false;
    bool shouldHideFiberSurface = false;
    bool shouldHideCombinedSurface = false;
    bool shouldFocusSelectedObject = false;

    // std::vector<std::vector<std::vector<int>>> volumes;
    // std::vector<std::vector<std::vector<int>>> volumesF;

    void computeFiberSurface();

    double isosurfaceOpacity = 1.0;
    double fibersurfaceOpacity = 1.0;
    double cartesianSurfaceOpacity = 1.0;

    void generateDisplayList();

  protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

    // Events
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);
    void mouseDoubleClickEvent(QMouseEvent*);
    void keyPressEvent(QKeyEvent* event);

  private:
    Data* data;
    QWidget* sibling;
    HistogramWidget* histogramSibling;

    // Arcball stuff
    QPointF position;
    BallData theBall;

    // Render Triangles
    void drawSolidTriangle(GLfloat vertices[3][3]);

    // Render Various Functions
    void cube();
    void drawAxis(GLfloat, GLfloat);
    void drawWiredCube(const GLfloat vertices[8][3]);

    void drawScene();

    float translateX = 0.;
    float translateY = 0.;

    float initialX = -1.;
    float initialY = -1.;

    // Utility
    void setMaterial(GLfloat, GLfloat, GLfloat, GLfloat, GLfloat);

}; // class GLPolygonWidget
