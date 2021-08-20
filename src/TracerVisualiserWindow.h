#pragma once

#include <QKeyEvent>
#include <QWidget>
#include <QGridLayout>
#include <QSlider>


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

    void keyPressEvent(QKeyEvent* event);

    Data* data;

    PlotWidget* plotWidget;
    TracerVisualiserWidget* tracerVisualiserWidget;

    QGridLayout* windowLayout;

    QGridLayout *optionsLayout;
    QGridLayout * optionsLayout2;

    QCheckBox *checkboxShowVertices;
    QCheckBox *checkboxShowEdges;
    QCheckBox *checkboxShowFaces;

    QSlider *faceOpacitySlider;
    QSlider *edgeOpacitySlider;
    QSlider *vertexOpacitySlider;

    QSlider *fakeSlider;


    QCheckBox *checkbox2;
};
