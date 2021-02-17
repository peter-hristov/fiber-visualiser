#include <qpoint.h>
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
#include <QKeyEvent>
#include <QtMath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include <qnamespace.h>
#include <queue>
#include <stdlib.h> /* srand, rand */
#include <string>
#include <time.h> /* time */
#include <typeinfo>
#include <vector>

#include "./TracerVisualiserWidget.h"
#include "./TracerVisualiserWindow.h"
#include "./utility/Geometry.h"
#include "./utility/TetrahedronDepth.h"
#include "./utility/Utility.h"
#include "src/utility/MarchingCubes.h"
#include "src/utility/MeshTriangle.h"
#include "src/utility/SurfaceMesh.h"

using namespace std;
using namespace tv9k::utility;

void
TracerVisualiserWidget::setMaterial(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha, GLfloat shininess)
{
    glMaterialfv(GL_FRONT, GL_AMBIENT, &(vector<GLfloat>({ red, green, blue, alpha })[0]));
    glMaterialfv(GL_FRONT, GL_DIFFUSE, &(vector<GLfloat>({ red, green, blue, alpha })[0]));
    glMaterialfv(GL_FRONT, GL_SPECULAR, &(vector<GLfloat>({ red, green, blue, alpha })[0]));
    glMaterialf(GL_FRONT, GL_SHININESS, shininess);
}

TracerVisualiserWidget::TracerVisualiserWidget(QWidget* parent,
                                               QWidget* _sibling,
                                               QWidget* _histogramSibling,
                                               Data* _data)
  : QGLWidget(parent)
{
    // Set points to other singletons
    this->data = _data;
    this->sibling = _sibling;
    this->histogramSibling = dynamic_cast<HistogramWidget*>(_histogramSibling);

    // Default values for paraters
    this->scale = data->zdim * 8;

    // Initialise Arcball
    Ball_Init(&theBall);
    Ball_Place(&theBall, qOne, 1);
}

void
TracerVisualiserWidget::computeFiberSurface()
{}

void
TracerVisualiserWidget::initializeGL()
{
    glClearColor(0.3, 0.3, 0.3, 0.0);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    // Enable alpha blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Do not enable this, doesn't work properly
    // glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90, 1, 0.01, 10000);
}

void
TracerVisualiserWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
}

void
TracerVisualiserWidget::generateDisplayList(const SurfaceMesh& mesh, const SurfaceType surfaceType)
{
    if (surfaceType == SurfaceType::isosurface) {
        glDeleteLists(displayListIndex, 1);
        displayListIndex = glGenLists(1);
        glNewList(displayListIndex, GL_COMPILE);

    } else if (surfaceType == SurfaceType::fibersurface) {
        glDeleteLists(displayListIndexF, 1);
        displayListIndexF = glGenLists(1);
        glNewList(displayListIndexF, GL_COMPILE);
    } else if (surfaceType == SurfaceType::combinedSurface) {
        glDeleteLists(displayListIndexC, 1);
        displayListIndexC = glGenLists(1);
        glNewList(displayListIndexC, GL_COMPILE);
    }
    else {
        // We shouldn't be here
        assert(false);
    }

    glBegin(GL_TRIANGLES);

    for (int i = 0; i < mesh.triangles.size(); i++) {

        double focusedOpacity;

        if (surfaceType == SurfaceType::isosurface) {
            focusedOpacity = this->isosurfaceOpacity;
        } else if (surfaceType == SurfaceType::fibersurface) {
            focusedOpacity = this->fibersurfaceOpacity;
        } else if (surfaceType == SurfaceType::combinedSurface) {
            focusedOpacity = this->cartesianSurfaceOpacity;
        }
        else {
            // We shouldn't be here
            assert(false);
        }

        const auto nonFocusedOpacity = 0.7;

        auto opacity = nonFocusedOpacity;

        // The statement isIsosurface == this->data->isIsoContourSelected holds when
        // 1. the selected object is from the isosurface and we are currently drawing the isosurface
        // 2. the selected object is from the fibersurface and we are currently drawing the fibersurface
        const bool areWeSelectedObject =
          (surfaceType == this->data->selectedSurfaceType[this->data->currentTimestep] &&
           mesh.triangles[i].triangleId == this->data->selectedID[this->data->currentTimestep]);

        // Full opacity if we're the selecte object
        if (
          // No selected object, so everything is selected
          (-1 == this->data->selectedID[this->data->currentTimestep]) ||
          // If we're drawing something not selected
          (true == areWeSelectedObject)) {
            opacity = focusedOpacity;
        }

        // Set colour (this is the ID of the object), +1 is to avoid (0, 0, 0) which
        // is ambiguous for both
        if (surfaceType == SurfaceType::isosurface) {
            glColor4ub(mesh.triangles[i].triangleId + 1, 0, 0, 255);
        } else if (surfaceType == SurfaceType::fibersurface) {
            glColor4ub(0, mesh.triangles[i].triangleId + 1, 0, 255);
        } else if (surfaceType == SurfaceType::combinedSurface) {
            glColor4ub(0, 0, mesh.triangles[i].triangleId + 1, 255);
        }

        // Set up the vertices and the normal of the triangle
        GLfloat vertices[3][3] = { mesh.triangles[i].vertices[0][0], mesh.triangles[i].vertices[0][1],
                                   mesh.triangles[i].vertices[0][2], mesh.triangles[i].vertices[1][0],
                                   mesh.triangles[i].vertices[1][1], mesh.triangles[i].vertices[1][2],
                                   mesh.triangles[i].vertices[2][0], mesh.triangles[i].vertices[2][1],
                                   mesh.triangles[i].vertices[2][2] };

        GLfloat normals[3][3] = { mesh.triangles[i].normals[0][0], mesh.triangles[i].normals[0][1],
                                  mesh.triangles[i].normals[0][2], mesh.triangles[i].normals[0][0],
                                  mesh.triangles[i].normals[1][1], mesh.triangles[i].normals[1][2],
                                  mesh.triangles[i].normals[0][0], mesh.triangles[i].normals[2][1],
                                  mesh.triangles[i].normals[2][2] };

        if (data->flatNormals) {
            drawSolidTriangle(vertices);
        } else {

            const auto plotWidget = dynamic_cast<PlotWidget*>(sibling);

            if (surfaceType == SurfaceType::isosurface && true == this->shouldIntersect &&
                plotWidget->isPointInsidePolygon(mesh.triangles[i].projectedVertices[0])) {
                setMaterial(1, 0.3, 0, opacity, 0.0);
            } else {
                setMaterial(
                  mesh.triangles[i].color[0], mesh.triangles[i].color[1], mesh.triangles[i].color[2], opacity, 30.0);
            }

            glNormal3fv(normals[0]);
            glVertex3fv(vertices[0]);

            if (surfaceType == SurfaceType::isosurface && true == this->shouldIntersect &&
                plotWidget->isPointInsidePolygon(mesh.triangles[i].projectedVertices[1])) {
                setMaterial(1, 0.3, 0, opacity, 0.0);
            } else {
                setMaterial(
                  mesh.triangles[i].color[0], mesh.triangles[i].color[1], mesh.triangles[i].color[2], opacity, 30.0);
            }

            glNormal3fv(normals[1]);
            glVertex3fv(vertices[1]);

            if (surfaceType == SurfaceType::isosurface && true == this->shouldIntersect &&
                plotWidget->isPointInsidePolygon(mesh.triangles[i].projectedVertices[2])) {
                setMaterial(1, 0.3, 0, opacity, 0.0);
            } else {
                setMaterial(
                  mesh.triangles[i].color[0], mesh.triangles[i].color[1], mesh.triangles[i].color[2], opacity, 30.0);
            }

            glNormal3fv(normals[2]);
            glVertex3fv(vertices[2]);
        }
    }
    glEnd();
    glEndList();
}

void
TracerVisualiserWidget::drawSolidTriangle(GLfloat vertices[3][3])
{
    GLfloat a[] = { vertices[1][0] - vertices[0][0], vertices[1][1] - vertices[0][1], vertices[1][2] - vertices[0][2] };

    GLfloat b[] = { vertices[2][0] - vertices[0][0], vertices[2][1] - vertices[0][1], vertices[2][2] - vertices[0][2] };

    // Invert Normal for sub/super level sets
    GLfloat normal[] = { isovalueMult * (a[2] * b[1] - b[1] * b[2]),
                         isovalueMult * (a[0] * b[2] - a[2] * b[0]),
                         isovalueMult * (a[1] * b[0] - a[0] * b[1]) };

    GLfloat normalFlipped[] = { -1 * normal[0], -1 * normal[1], -1 * normal[2] };

    glNormal3fv(normalFlipped);
    glVertex3fv(vertices[0]);
    glVertex3fv(vertices[1]);
    glVertex3fv(vertices[2]);
}

void
TracerVisualiserWidget::drawAxis(GLfloat length, GLfloat width)
{
    GLUquadric* line = gluNewQuadric();

    // Z Axis
    setMaterial(0., 0., 1., 1.0, 30.);
    gluCylinder(line, width, width, length, 100, 100);

    // Y Axis
    setMaterial(0., 1., 0., 1.0, 30.);
    glPushMatrix();
    glRotatef(-90, 1, 0, 0);
    gluCylinder(line, width, width, length, 100, 100);
    glPopMatrix();

    // X Axis
    setMaterial(1., 0., 0., 1.0, 30.);
    glPushMatrix();
    glRotatef(90, 0, 1, 0);
    gluCylinder(line, width, width, length, 100, 100);
    glPopMatrix();

    gluDeleteQuadric(line);
}

void
TracerVisualiserWidget::drawWiredCube(const GLfloat vertices[8][3])
{
    glBegin(GL_LINES);
    {
        // Front Side
        glVertex3fv(vertices[0]);
        glVertex3fv(vertices[2]);

        glVertex3fv(vertices[2]);
        glVertex3fv(vertices[3]);

        glVertex3fv(vertices[1]);
        glVertex3fv(vertices[3]);

        glVertex3fv(vertices[0]);
        glVertex3fv(vertices[1]);

        // Back Side
        glVertex3fv(vertices[4]);
        glVertex3fv(vertices[6]);

        glVertex3fv(vertices[6]);
        glVertex3fv(vertices[7]);

        glVertex3fv(vertices[5]);
        glVertex3fv(vertices[7]);

        glVertex3fv(vertices[4]);
        glVertex3fv(vertices[5]);

        // Sides
        glVertex3fv(vertices[0]);
        glVertex3fv(vertices[4]);

        glVertex3fv(vertices[2]);
        glVertex3fv(vertices[6]);

        glVertex3fv(vertices[3]);
        glVertex3fv(vertices[7]);

        glVertex3fv(vertices[1]);
        glVertex3fv(vertices[5]);
    }
    glEnd();
}

void
TracerVisualiserWidget::drawScene()
{
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Rescale normals when scaling
    glEnable(GL_RESCALE_NORMAL);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    GLfloat light_pos[] = { 0, 0, 10000, 1. };
    glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
    glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 180.);

    // Scrolling
    glTranslatef(0.0, 0.0, -1 * scale / 2);

    // Offset along x, y (with the right mouse button)
    glTranslatef(translateX, translateY, 0.0);

    // Arcball
    GLfloat mNow[16];
    Ball_Value(&theBall, mNow);
    glMultMatrixf(mNow);

    // GLUquadric* sphere = gluNewQuadric();
    // gluSphere(sphere, 2, 3, 3);

    glTranslatef((-1.0 * this->data->xdim) / 2.0, 0, 0);
    glTranslatef(0, 0, (1.0 * this->data->ydim) / 2.0);
    glTranslatef(0, (-1.0 * this->data->zdim) / 2.0, 0);

    glRotatef(-90., 1., 0., 0.);

    this->drawAxis(1000., 1.0 * this->data->xdim / 800.0);

    // Bounding Cube
    glPushMatrix();
    {
        glScalef(this->data->xdim - 1, this->data->ydim - 1, this->data->zdim - 1);
        setMaterial(255, 255, 255, 100, 30.0);
        this->drawWiredCube(tv9k::geometry::cubeVertices);
    }
    glPopMatrix();


    glFlush();
}

void
TracerVisualiserWidget::paintGL()
{
    this->drawScene();
}

void
TracerVisualiserWidget::wheelEvent(QWheelEvent* event)
{
    this->scale -= event->angleDelta().y() / 10;

    if (scale <= 0) {
        scale = 0;
    }

    this->update();
}

void
TracerVisualiserWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) {
        HVect vNow;
        vNow.x = (2.0 * event->x() - width()) / width();
        vNow.y = (height() - 2.0 * event->y()) / height();

        Ball_Mouse(&theBall, vNow);
        Ball_BeginDrag(&theBall);

        this->update();
    }
    if (event->button() == Qt::RightButton) {
        initialX = event->localPos().x();
        initialY = event->localPos().y();
    }
}

void
TracerVisualiserWidget::mouseMoveEvent(QMouseEvent* event)
{
    if (event->buttons() == Qt::LeftButton) {
        HVect vNow;
        vNow.x = (2.0 * event->localPos().x() - width()) / width();
        vNow.y = (height() - 2.0 * event->y()) / height();

        Ball_Mouse(&theBall, vNow);
        Ball_Update(&theBall);

        this->update();
    } else if (event->buttons() == Qt::RightButton) {
        float x = event->localPos().x();
        float y = event->localPos().y();

        translateX -= (initialX - x) / 10;
        translateY += (initialY - y) / 10;

        initialX = event->localPos().x();
        initialY = event->localPos().y();

        this->update();
    }
}

void
TracerVisualiserWidget::mouseReleaseEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) {
        Ball_EndDrag(&theBall);
        this->update();
    }
}

void
TracerVisualiserWidget::keyPressEvent(QKeyEvent* event)
{
    if (event->key() == Qt::Key_W) {
        translateX += 0.1 * this->data->xdim;
        this->update();
    }
    if (event->key() == Qt::Key_A) {
        translateY += 0.1 * this->data->xdim;
        this->update();
    }
    if (event->key() == Qt::Key_S) {
        translateX -= 0.1 * this->data->xdim;
        this->update();
    }
    if (event->key() == Qt::Key_D) {
        translateY -= 0.1 * this->data->xdim;
        this->update();
    }
}

void
TracerVisualiserWidget::mouseDoubleClickEvent(QMouseEvent* event)
{
}
