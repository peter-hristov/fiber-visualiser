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

#include "./Data.h"
#include "./ArcBall/Ball.h"

class TracerVisualiserWidget : public QGLWidget
{
  public:
    TracerVisualiserWidget(QWidget*, QWidget*, Data*);
    GLfloat scale = 0;

    GLfloat isovalueMult = -1.0;

    int displayListIndex = 0;
    void generateDisplayList();

    int displayListIndexTriangles = 0;
    void generateDisplayListTriangles(const float);

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
