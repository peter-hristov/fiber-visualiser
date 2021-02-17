#pragma once

#include <QBoxLayout>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QGLWidget>
#include <QKeyEvent>
#include <QMenuBar>
#include <QObject>
#include <QProgressBar>
#include <QSlider>
#include <qcheckbox.h>
#include <qlabel.h>

#include "./Data.h"
#include "./GlobalConfig.h"
#include "./PlotWidget.h"
#include "./TracerVisualiserWidget.h"

class TracerVisualiserWindow : public QWidget
{
    Q_OBJECT
  public:
    TracerVisualiserWindow(QWidget*, Data*, tv9k::InputInformation);
    ~TracerVisualiserWindow();


    Data* data;

    PlotWidget* plotWidget;
    TracerVisualiserWidget* tracerVisualiserWidget;

    QGridLayout* windowLayout;
};
