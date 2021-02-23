#include <cmath>
#include <iostream>

#include <QCheckBox>
#include <QComboBox>
#include <QPushButton>
#include <QTimer>
#include <qcheckbox.h>
#include <qgridlayout.h>
#include <qlabel.h>
#include <qnamespace.h>
#include <qprogressbar.h>
#include <qpushbutton.h>
#include <utility>

#include "./TracerVisualiserWindow.h"
#include "src/Data.h"

using namespace std;
void
TracerVisualiserWindow::keyPressEvent(QKeyEvent* event)
{
    if (this->data->mousePoints.empty())
        return;

    int speed = 2;

    if (event->key() == Qt::Key_I) {
        this->data->mousePoints[0].setY(this->data->mousePoints[0].y() - speed);
        this->update();
    }
    if (event->key() == Qt::Key_J) {
        this->data->mousePoints[0].setX(this->data->mousePoints[0].x() - speed);
        this->update();
    }
    if (event->key() == Qt::Key_K) {
        this->data->mousePoints[0].setY(this->data->mousePoints[0].y() + speed);
        this->update();
    }
    if (event->key() == Qt::Key_L) {
        this->data->mousePoints[0].setX(this->data->mousePoints[0].x() + speed);
        this->update();
    }
}

TracerVisualiserWindow::TracerVisualiserWindow(QWidget* parent,
                                               Data* _data,
                                               tv9k::InputInformation input
                                               )
  : QWidget(parent)
{
    this->data = _data;

    // Initialise Widgets
    plotWidget = new PlotWidget(this, data, "none", input);
    tracerVisualiserWidget = new TracerVisualiserWidget(this, plotWidget, data);

    // Set up layout
    windowLayout = new QGridLayout(this);
    windowLayout->addWidget(tracerVisualiserWidget, 1, 0);
    windowLayout->addWidget(plotWidget, 1, 1);
}

TracerVisualiserWindow::~TracerVisualiserWindow()
{
    delete plotWidget;
    delete windowLayout;
    delete tracerVisualiserWidget;
}
