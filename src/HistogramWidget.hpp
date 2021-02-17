#pragma once

#include "Data.h"
#include <qwidget.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <QPainter>
#include <QWidget>

#include <vector>

class HistogramWidget : public QWidget
{
  public:
    Data* data;
    int scaleX;
    int scaleY = 200;

    HistogramWidget()
      : QWidget()
    {
        setFocusPolicy(Qt::StrongFocus);
    }

    HistogramWidget(Data* _data, const size_t _resolution)
      : data(_data)
      , scaleX(_resolution)
    {}

    void paintEvent(QPaintEvent*)
    {
    }
};
