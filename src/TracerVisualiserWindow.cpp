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

using namespace std;

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
