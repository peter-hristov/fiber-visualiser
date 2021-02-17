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
#include "./utility/Utility.h"
#include "src/Data.h"
#include "src/HistogramWidget.hpp"

using namespace std;

void
TracerVisualiserWindow::keyPressEvent(QKeyEvent* event)
{
}

TracerVisualiserWindow::TracerVisualiserWindow(QWidget* parent,
                                               Data* _data,
                                               tv9k::InputInformation input,
                                               string interpolationType)
  : QWidget(parent)
{
    this->data = _data;

    //
    // Initialise Widgets
    //
    histogramWidget = new HistogramWidget(data, input.histogramResolution);
    plotWidget = new PlotWidget(this, data, interpolationType, input);
    tracerVisualiserWidget = new TracerVisualiserWidget(this, plotWidget, histogramWidget, data);


    windowLayout = new QGridLayout(this);
    windowLayout->addWidget(tracerVisualiserWidget, 1, 0);
    windowLayout->addWidget(plotWidget, 1, 1);

}

void
TracerVisualiserWindow::setTimestep(int timestep)
{
}

TracerVisualiserWindow::~TracerVisualiserWindow()
{
    delete plotWidget;
    delete windowLayout;
    delete tracerVisualiserWidget;
}
