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
        QPainter p(this);
        p.save();
        p.setTransform(QTransform(1., 0., 0., -1., 0., (float)scaleY));

        p.setWindow(QRect(0, 0, scaleX, scaleY));

        p.setPen(Qt::black);

        // for (int i = 0 ; i < scaleX ; i++)
        //{
        // histogram[i] = 0;
        //}

        // for (int i = 0 ; i < scaleX/2 ; i++)
        //{
        // histogram[i] = i;
        // histogram[scaleX - i] = i;
        //}

        // for (int i = scaleX/2 ; i < scaleX ; i++)
        //{
        // histogram[i] = 100 - i;
        //}

        const auto& histogram = this->data->histograms[this->data->currentTimestep];

        unsigned int max = histogram[0];
        for (int i = 0; i < scaleX; i++) {
            max = std::max(max, histogram[i]);
        }

        // std::cout << "THE MAXIMUM IS " << max << endl;

        for (int i = 0; i < scaleX; i++) {
            int value = histogram[i];

            float scaledHeight = ((float)value / (float)max) * (float)scaleY;

            // cout << scaledHeight << endl;

            // cout << value << " " << max << " " << scaledHeight << endl;
            // std::cout << "H = " << scaledHeight << std::endl;

            p.fillRect(i, 0, 1, scaledHeight, Qt::SolidPattern);
        }

        p.restore();

        // auto pen = QPen(Qt::black);
        // pen.setWidth(5.5);

        // p.setPen(pen);
        // float boxOffset = 10;

        // p.drawLine(0, 100, scaleX, 100);
        // p.drawLine(boxOffset, scale - boxOffset , scale - boxOffset, scale -
        // boxOffset);
    }
};
