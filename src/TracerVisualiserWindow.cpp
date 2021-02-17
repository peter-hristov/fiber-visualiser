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
    if (event->key() == Qt::Key_D) {
        timeStepSlider->setValue(timeStepSlider->value() + 1);
    }

    else if (event->key() == Qt::Key_A) {
        timeStepSlider->setValue(timeStepSlider->value() - 1);
    }
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

    //
    // Initialise Controlls
    //
    saveISObjBtn = new QPushButton("Save Isosurface Objects");
    saveFSObjBtn = new QPushButton("Save Fiber Surface Objects");

    isovalueSpinBox = new QDoubleSpinBox();
    isovalueSpinBox->setMinimum(data->isoField->min);
    isovalueSpinBox->setMaximum(data->isoField->max);
    double step = (data->isoField->min - data->isoField->max) / 1000.;
    isovalueSpinBox->setSingleStep(step);
    isovalueSpinBox->setValue(data->isoField->min);
    isovalueSpinBox->setKeyboardTracking(false);
    // @TODO Hardcoded
    isovalueSpinBox->setDecimals(4);

    scatterPlotTypeComboBox = new QComboBox();
    scatterPlotTypeComboBox->insertItem(0, "Discrete");
    scatterPlotTypeComboBox->insertItem(1, "Continuous");

    isoFieldComboBox = new QComboBox();
    uFieldComboBox = new QComboBox();
    vFieldComboBox = new QComboBox();

    int counter = 0;
    for (const auto& field : this->data->scalarFields) {
        isoFieldComboBox->insertItem(counter++, field.first.c_str());
        uFieldComboBox->insertItem(counter, field.first.c_str());
        vFieldComboBox->insertItem(counter, field.first.c_str());

        if (field.first.c_str() == this->data->isoField->name) {
            isoFieldComboBox->setCurrentIndex(counter - 1);
        }
        if (field.first.c_str() == this->data->uField->name) {
            uFieldComboBox->setCurrentIndex(counter - 1);
        }
        if (field.first.c_str() == this->data->vField->name) {
            vFieldComboBox->setCurrentIndex(counter - 1);
        }
    }

    projectionTypeComboBox = new QComboBox();
    projectionTypeComboBox->insertItem(0, "None");
    projectionTypeComboBox->insertItem(1, "Interior Points");
    projectionTypeComboBox->insertItem(2, "Isosurface Triangles");
    projectionTypeComboBox->insertItem(3, "Tetrahedra");

    projectionTypeComboBox->setCurrentIndex(this->data->projectionType);

    mergeTreeTypeComboBox = new QComboBox();
    mergeTreeTypeComboBox->insertItem(0, "Join Tree");
    mergeTreeTypeComboBox->insertItem(1, "Split Tree");

    isosurfaceOpacitySlider = new QSlider(Qt::Horizontal);
    isosurfaceOpacitySlider->setValue(100);
    isosurfaceOpacitySlider->setTracking(false);

    toggleIsosurface = new QCheckBox("Show");
    toggleIsosurface->setChecked(true);
    isoSurfaceClear = new QPushButton("Clear");

    fibersurfaceOpacitySlider = new QSlider(Qt::Horizontal);
    fibersurfaceOpacitySlider->setValue(100);
    fibersurfaceOpacitySlider->setTracking(false);

    toggleFibersurface = new QCheckBox("Show");
    toggleFibersurface->setChecked(true);
    toggleIntersection = new QCheckBox("Intersection");
    toggleIntersection->setChecked(false);
    fiberSurfaceClear = new QPushButton("Clear");

    combinedSurfaceOpacitySlider = new QSlider(Qt::Horizontal);
    combinedSurfaceOpacitySlider->setValue(100);
    combinedSurfaceOpacitySlider->setTracking(false);
    toggleCombinedSurface = new QCheckBox("Show");
    toggleCombinedSurface->setChecked(true);
    combinedSurfaceComputeButton = new QPushButton("Compute");


    isosurfaceSimplificationSlider = new QSlider(Qt::Horizontal);
    isosurfaceSimplificationSlider->setTracking(false);
    isosurfaceSimplificationSlider->setMaximum(1000);
    isosurfaceSimplificationSlider->setValue(1000);

    fiberSurfaceSimplificationSlider = new QSlider(Qt::Horizontal);
    fiberSurfaceSimplificationSlider->setTracking(false);
    fiberSurfaceSimplificationSlider->setMaximum(1000);
    fiberSurfaceSimplificationSlider->setValue(1000);

    scaleSlider = new QSlider(Qt::Horizontal);
    scaleSlider->setTracking(false);

    isovalueSlider = new QSlider(Qt::Horizontal);
    isovalueSlider->setMinimum(0);
    isovalueSlider->setMaximum(500);
    isovalueSlider->setTickInterval(1);
    isovalueSlider->setSingleStep(1);
    isovalueSlider->setTracking(false);

    timeStepSlider = new QSlider(Qt::Horizontal);
    timeStepSlider->setMinimum(0);
    timeStepSlider->setMaximum(data->tdim - 1);
    timeStepSlider->setTracking(false);

    timeStepPreviousButton = new QPushButton("<");
    timeStepPreviousButton->setMaximumWidth(40);

    timeStepNextButton = new QPushButton(">");
    timeStepNextButton->setMaximumWidth(40);

    isovalueRangeLabel = new QLabel("Isovalue in [" + QString::number(data->isoField->min).mid(0, 6) + ", " +
                                    QString::number(data->isoField->max).mid(0, 6) + "]");

    // isovalueRangeLabel = new QLabel();

    ////QMovie* spinnerMovie = new QMovie("./external/loading.gif");
    // QMovie* spinnerMovie = new QMovie("~/Projects/tracer-visualiser/external/loading.gif");
    // isovalueRangeLabel->setMovie(spinnerMovie);
    // spinnerMovie->start();

    //
    // Set up Layout
    //

    // Left Layout

    auto fieldPickLayout = new QGridLayout();
    fieldPickLayout->addWidget(new QLabel("Iso Field"), 0, 0);
    fieldPickLayout->addWidget(isoFieldComboBox, 0, 1);
    fieldPickLayout->addWidget(new QLabel("Fiber U Field"), 0, 2);
    fieldPickLayout->addWidget(uFieldComboBox, 0, 3);
    fieldPickLayout->addWidget(new QLabel("Fiber V Field"), 0, 4);
    fieldPickLayout->addWidget(vFieldComboBox, 0, 5);

    comboBoxLayout = new QGridLayout();

    // comboBoxLayout->addWidget(saveISObjBtn, 0, 0);
    // comboBoxLayout->addWidget(saveFSObjBtn, 0, 1);

    comboBoxLayout->addWidget(new QLabel("Scatterplot Type"), 1, 0);
    comboBoxLayout->addWidget(scatterPlotTypeComboBox, 1, 1);
    comboBoxLayout->addWidget(new QLabel("Projection Type"), 2, 0);
    comboBoxLayout->addWidget(projectionTypeComboBox, 2, 1);
    comboBoxLayout->addWidget(new QLabel("Merge Tree Type"), 3, 0);
    comboBoxLayout->addWidget(mergeTreeTypeComboBox, 3, 1);

    isovalueLayout = new QGridLayout();
    isovalueLayout->addWidget(this->isovalueRangeLabel, 1, 0);
    isovalueLayout->addWidget(isovalueSpinBox, 1, 1);

    histogramLayout = new QGridLayout();

    // histogramLayout->addLayout(fieldPickLayout, 0, 0);
    histogramLayout->addWidget(histogramWidget, 0, 0);

    leftControlsLayout = new QGridLayout();
    leftControlsLayout->addLayout(comboBoxLayout, 0, 0);
    leftControlsLayout->addLayout(histogramLayout, 0, 1);
    leftControlsLayout->addLayout(isovalueLayout, 1, 0);
    leftControlsLayout->addWidget(isovalueSlider, 1, 1);

    auto leftControlsLayoutParent = new QGridLayout();
    leftControlsLayoutParent->addLayout(fieldPickLayout, 0, 0);
    leftControlsLayoutParent->addLayout(leftControlsLayout, 1, 0);

    // Right Layout
    timestepLayout = new QGridLayout();

    this->timeLabel = new QLabel(tr("Timestep: ") + QString::number(this->data->currentTimestep) + tr("/") +
                                 QString::number(this->data->tdim - 1) + " (" +
                                 QString::number(this->data->tVals[this->data->currentTimestep]) + ")");

    // isovalueRangeLabel->setText("Isovalue in [" + QString::number(data->isoField->min).mid(0, 6) + ", " +
    // QString::number(data->isoField->max).mid(0, 6) + "]");

    // timestepLayout->addWidget(timestepLabel, 0, 0);
    timestepLayout->addWidget(timeStepPreviousButton, 0, 0);
    timestepLayout->addWidget(timeStepNextButton, 0, 1);
    timestepLayout->addWidget(timeStepSlider, 0, 2);

    rightControlsLayout = new QGridLayout();
    rightControlsLayout->addWidget(timeLabel, 0, 0);
    rightControlsLayout->addLayout(timestepLayout, 0, 1);

    isosurfaceOpacityLayout = new QGridLayout();
    isosurfaceOpacityLayout->addWidget(isosurfaceOpacitySlider, 0, 0);
    isosurfaceOpacityLayout->addWidget(toggleIntersection, 0, 1);
    isosurfaceOpacityLayout->addWidget(isoSurfaceClear, 0, 2);
    isosurfaceOpacityLayout->addWidget(toggleIsosurface, 0, 3);

    fibersurfaceOpacityLayout = new QGridLayout();
    fibersurfaceOpacityLayout->addWidget(fibersurfaceOpacitySlider, 0, 0);
    fibersurfaceOpacityLayout->addWidget(fiberSurfaceClear, 0, 1);
    fibersurfaceOpacityLayout->addWidget(toggleFibersurface, 0, 2);

    combinedSurfaceOpacityLayout = new QGridLayout();
    combinedSurfaceOpacityLayout->addWidget(combinedSurfaceOpacitySlider, 0, 0);
    combinedSurfaceOpacityLayout->addWidget(combinedSurfaceComputeButton, 0, 1);
    combinedSurfaceOpacityLayout->addWidget(toggleCombinedSurface, 0, 2);

    this->loadingLabel = new QLabel("");

    rightControlsLayout->addWidget(new QLabel(tr("Isosurface Opacity:")), 1, 0);
    rightControlsLayout->addLayout(isosurfaceOpacityLayout, 1, 1);
    rightControlsLayout->addWidget(new QLabel(tr("Fiber Surface Opacity:")), 2, 0);
    rightControlsLayout->addLayout(fibersurfaceOpacityLayout, 2, 1);
    rightControlsLayout->addWidget(new QLabel(tr("Cartesian Surface Opacity:")), 3, 0);
    rightControlsLayout->addLayout(combinedSurfaceOpacityLayout, 3, 1);
    //rightControlsLayout->addWidget(isosurfaceSimplificationSlider, 3, 1);
    rightControlsLayout->addWidget(loadingLabel, 4, 0);

    this->progressBar = new QProgressBar();
    rightControlsLayout->addWidget(progressBar, 4, 1);

    windowLayout = new QGridLayout(this);
    windowLayout->addWidget(tracerVisualiserWidget, 1, 0);
    windowLayout->addWidget(plotWidget, 1, 1);
    windowLayout->addLayout(leftControlsLayoutParent, 2, 0);
    windowLayout->addLayout(rightControlsLayout, 2, 1);

    //
    // Signals and Slots
    //

    connect(scatterPlotTypeComboBox, QOverload<int>::of(&QComboBox::activated), [=](int index) {
        cout << "Opperation NOT SUPPORTED!" << index << endl;
    });

    connect(uFieldComboBox, QOverload<int>::of(&QComboBox::activated), [=](int index) {
        std::string selectedField = uFieldComboBox->currentText().toStdString();

        this->loadingLabel->setText("LOADING....");
        this->progressBar->setValue(50);
        // Change v field
        this->data->changeUField(selectedField);

        // Recompute all the scatterplots
        this->data->computeScatterplots(input.scatterplotResolution);

        // Delete FSCP
        //this->plotWidget->mousePoints.clear();
        this->plotWidget->polyPoints.clear();

        // Clear FS mesh
        this->data->fibersurfaceMeshes.clear();
        this->data->fibersurfaceMeshes.resize(this->data->tdim);

        // Clear combination mesh
        this->data->combinedMeshes.clear();
        this->data->combinedMeshes.resize(this->data->tdim);

        // Update 3D triangles
        // this->tracerVisualiserWidget->generateDisplayList(this->data->isosurfaceMeshes[this->data->currentTimestep],
        // true);
        // @TODO wasteful, just project triangles again
        this->data->computeMeshes(this->data->isoField->currentIsovalue);
        this->tracerVisualiserWidget->generateDisplayList(this->data->fibersurfaceMeshes[this->data->currentTimestep],
                                                          SurfaceType::fibersurface);
        this->tracerVisualiserWidget->generateDisplayList(this->data->combinedMeshes[this->data->currentTimestep],
                                                          SurfaceType::combinedSurface);

        this->loadingLabel->setText("");
        this->progressBar->setValue(0);
        this->tracerVisualiserWidget->update();

        // Update 2D triangles
        this->plotWidget->update();
    });

    connect(vFieldComboBox, QOverload<int>::of(&QComboBox::activated), [=](int index) {
        std::string selectedField = vFieldComboBox->currentText().toStdString();

        this->loadingLabel->setText("LOADING....");
        this->progressBar->setValue(50);
        // Change v field
        this->data->changeVField(selectedField);

        // Recompute all the scatterplots
        this->data->computeScatterplots(input.scatterplotResolution);

        // Delete FSCP
        //this->plotWidget->mousePoints.clear();
        this->plotWidget->polyPoints.clear();

        // Clear FS mesh
        this->data->fibersurfaceMeshes.clear();
        this->data->fibersurfaceMeshes.resize(this->data->tdim);

        // Clear combination mesh
        this->data->combinedMeshes.clear();
        this->data->combinedMeshes.resize(this->data->tdim);

        // Update 3D triangles
        // this->tracerVisualiserWidget->generateDisplayList(this->data->isosurfaceMeshes[this->data->currentTimestep],
        // true);
        // @TODO wasteful, just projecte triangles again
        this->data->computeMeshes(this->data->isoField->currentIsovalue);
        this->tracerVisualiserWidget->generateDisplayList(this->data->fibersurfaceMeshes[this->data->currentTimestep],
                                                          SurfaceType::fibersurface);
        this->tracerVisualiserWidget->generateDisplayList(this->data->combinedMeshes[this->data->currentTimestep],
                                                          SurfaceType::combinedSurface);
        this->tracerVisualiserWidget->update();
        this->loadingLabel->setText("");
        this->progressBar->setValue(0);

        // Update 2D triangles
        this->plotWidget->update();
    });

    connect(isoFieldComboBox, QOverload<int>::of(&QComboBox::activated), [=](int index) {
        std::string selectedField = isoFieldComboBox->currentText().toStdString();
        // cout << "You have selected - " <<
        // isoFieldComboBox->currentText().toStdString() << endl;

        this->loadingLabel->setText("LOADING....");
        this->progressBar->setValue(50);

        // Change currentIsofield
        this->data->changeIsoField(selectedField, mergeTreeTypeComboBox->currentIndex());

        // Redraw histogram widget
        this->data->computeHistograms(this->histogramWidget->scaleX);

        // Reset object selection
        for (int i = 0; i < data->tdim; i++) {
            this->data->selectedSurfaceType[i] = SurfaceType::none;
            this->data->selectedID[i] = -1;
            this->data->selectedObjectMinMax[i] = { { 0, 0, 0 }, { 0, 0, 0 } };
        }

        // Compute join trees in needed
        if (false == this->data->isoField->hasJoinTrees()) {
            this->data->isoField->computeJoinTrees();
        }

        // Compute isosurfaces
        this->data->computeMeshes(this->data->isoField->currentIsovalue);
        this->data->computeCombinedMeshes(this->data->currentSignedDistanceField[this->data->currentTimestep], this->data->isoField->currentIsovalue);

        this->tracerVisualiserWidget->generateDisplayList(this->data->isosurfaceMeshes[this->data->currentTimestep],
                SurfaceType::isosurface);
        this->tracerVisualiserWidget->generateDisplayList(this->data->combinedMeshes[this->data->currentTimestep], SurfaceType::combinedSurface);

        // Change isovalue range label
        isovalueRangeLabel->setText("Isovalue in [" + QString::number(data->isoField->min).mid(0, 6) + ", " +
                                    QString::number(data->isoField->max).mid(0, 6) + "]");

        // Change isofield selection stuff (without triggering them)
        this->isovalueSpinBox->blockSignals(true);
        isovalueSpinBox->setMinimum(data->isoField->min);
        isovalueSpinBox->setMaximum(data->isoField->max);
        double step = (data->isoField->min - data->isoField->max) / 1000.;
        isovalueSpinBox->setSingleStep(step);
        isovalueSpinBox->setValue(this->data->isoField->currentIsovalue);
        this->isovalueSpinBox->blockSignals(false);

        // Set isovalues
        this->isovalueSlider->blockSignals(true);
        float t =
          (this->data->isoField->currentIsovalue - data->isoField->min) / (data->isoField->max - data->isoField->min);
        this->isovalueSlider->setValue(static_cast<int>(t * this->isovalueSlider->maximum()));
        this->isovalueSlider->blockSignals(false);

        // Update stuff
        this->tracerVisualiserWidget->update();
        this->histogramWidget->update();
        this->plotWidget->update();

        this->loadingLabel->setText("");
        this->progressBar->setValue(0);
    });

    connect(projectionTypeComboBox, QOverload<int>::of(&QComboBox::activated), [=](int index) {
        // 0 = None
        // 1 = Interior Points
        // 2 = Isosurface Triangles
        // 3 = Interior Tetrahedra

         data->projectionType = index;
         //this->tracerVisualiserWidget->computeIsosurface();

         this->plotWidget->update();
    });

    connect(mergeTreeTypeComboBox, QOverload<int>::of(&QComboBox::activated), [=](int index) {
        cerr << "Operation currently not supported!";
        // @TODO Deprecate
        // if (index == 0)
        //{
        ////if (false == this->data->isoField->hasJoinTrees())
        ////{
        ////this->data->isoField->computeJoinTrees();
        ////}
        ////this->data->tree =
        ///&this->data->isoField->joinTrees[timeStepSlider->value()];
        // this->tracerVisualiserWidget->isovalueMult = -1;
        //}
        // if (index == 1)
        //{
        ////if (false == this->data->isoField->hasSplitTrees())
        ////{
        ////this->data->isoField->computeSplitTrees();
        ////}
        ////this->tracerVisualiserWidget->isovalueMult = 1;
        ////this->data->tree =
        ///&this->data->isoField->splitTrees[timeStepSlider->value()];

        ////cout << "The current tree is the split" << endl;
        //}

        // this->tracerVisualiserWidget->recomputeIsosurface = true;
        // this->tracerVisualiserWidget->update();
    });

    connect(combinedSurfaceComputeButton, &QPushButton::clicked, this, [=]() {
            this->data->computeCombinedMeshes(this->data->currentSignedDistanceField[this->data->currentTimestep], this->data->isoField->currentIsovalue);
            this->tracerVisualiserWidget->generateDisplayList(this->data->combinedMeshes[this->data->currentTimestep], SurfaceType::combinedSurface);
            });

    connect(isoSurfaceClear, &QPushButton::clicked, this, [=]() {
        // Set isovalue
        this->data->isoField->currentIsovalue = this->data->isoField->min;

        // Set spinbox isovalue to min
        this->isovalueSpinBox->blockSignals(true);
        {
            isovalueSpinBox->setValue(this->data->isoField->min);
        }
        this->isovalueSpinBox->blockSignals(false);

        // Set slider value to 0
        this->isovalueSlider->blockSignals(true);
        {
            this->isovalueSlider->setValue(0);
        }
        this->isovalueSlider->blockSignals(false);

        // Clear meshes
        this->data->isosurfaceMeshes.clear();
        this->data->isosurfaceMeshes.resize(this->data->tdim);
        this->tracerVisualiserWidget->generateDisplayList(this->data->isosurfaceMeshes[this->data->currentTimestep],
                                                          SurfaceType::isosurface);

        this->data->combinedMeshes.clear();
        this->data->combinedMeshes.resize(this->data->tdim);
        this->tracerVisualiserWidget->generateDisplayList(this->data->combinedMeshes[this->data->currentTimestep],
                                                          SurfaceType::combinedSurface);

        // Update visuals
        this->plotWidget->update();
        this->tracerVisualiserWidget->update();
    });

    connect(fiberSurfaceClear, &QPushButton::clicked, this, [=]() {
        // Clear FSCP
            this->data->mousePoints[this->data->uField->name][this->data->vField->name].clear();

        // Clear Fibermesh
        this->data->fibersurfaceMeshes.clear();
        this->data->fibersurfaceMeshes.resize(this->data->tdim);
        this->tracerVisualiserWidget->generateDisplayList(this->data->fibersurfaceMeshes[this->data->currentTimestep],
                                                          SurfaceType::fibersurface);

        this->data->combinedMeshes.clear();
        this->data->combinedMeshes.resize(this->data->tdim);
        this->tracerVisualiserWidget->generateDisplayList(this->data->combinedMeshes[this->data->currentTimestep],
                                                          SurfaceType::combinedSurface);

        this->plotWidget->update();
        this->tracerVisualiserWidget->update();
    });

    // @TODO Fix these
    connect(toggleIntersection, &QPushButton::clicked, this, [=]() {
        this->tracerVisualiserWidget->shouldIntersect = !this->tracerVisualiserWidget->shouldIntersect;
        this->tracerVisualiserWidget->generateDisplayList(this->data->isosurfaceMeshes[this->data->currentTimestep],
                                                          SurfaceType::isosurface);
        this->tracerVisualiserWidget->update();
    });

    connect(toggleIsosurface, &QPushButton::clicked, this, [=]() {
        this->tracerVisualiserWidget->shouldHideIsosurface = !this->tracerVisualiserWidget->shouldHideIsosurface;
        this->tracerVisualiserWidget->update();
    });

    connect(toggleFibersurface, &QPushButton::clicked, this, [=]() {
        this->tracerVisualiserWidget->shouldHideFiberSurface = !this->tracerVisualiserWidget->shouldHideFiberSurface;
        this->tracerVisualiserWidget->update();
    });

    connect(toggleCombinedSurface, &QPushButton::clicked, this, [=]() {
        this->tracerVisualiserWidget->shouldHideCombinedSurface = !this->tracerVisualiserWidget->shouldHideCombinedSurface;
        this->tracerVisualiserWidget->update();
    });

    connect(timeStepNextButton, &QPushButton::clicked, this, [=]() {
        timeStepSlider->setValue(timeStepSlider->value() + 1);
    });

    connect(timeStepPreviousButton, &QPushButton::clicked, this, [=]() {
        timeStepSlider->setValue(timeStepSlider->value() - 1);
    });

    connect(timeStepSlider, &QSlider::valueChanged, plotWidget, [=](int value) { this->setTimestep(value); });

    connect(isosurfaceOpacitySlider, &QSlider::valueChanged, plotWidget, [=]() {
        this->tracerVisualiserWidget->isosurfaceOpacity = static_cast<double>(isosurfaceOpacitySlider->value()) / 100.0;
        this->tracerVisualiserWidget->generateDisplayList(this->data->isosurfaceMeshes[this->data->currentTimestep],
                                                          SurfaceType::isosurface);
        this->tracerVisualiserWidget->update();
    });

    connect(fibersurfaceOpacitySlider, &QSlider::valueChanged, plotWidget, [=]() {
        this->tracerVisualiserWidget->fibersurfaceOpacity =
          static_cast<double>(fibersurfaceOpacitySlider->value()) / 100.0;
        // this->tracerVisualiserWidget->computeFiberSurface();
        this->tracerVisualiserWidget->generateDisplayList(this->data->fibersurfaceMeshes[this->data->currentTimestep],
                                                          SurfaceType::fibersurface);
        this->tracerVisualiserWidget->update();
    });

    connect(combinedSurfaceOpacitySlider, &QSlider::valueChanged, plotWidget, [=]() {
        this->tracerVisualiserWidget->cartesianSurfaceOpacity = static_cast<double>(combinedSurfaceOpacitySlider->value()) / 100.0;

        this->tracerVisualiserWidget->generateDisplayList(this->data->combinedMeshes[this->data->currentTimestep],
                                                          SurfaceType::combinedSurface);
        this->tracerVisualiserWidget->update();
    });

    connect(isosurfaceSimplificationSlider, &QSlider::valueChanged, plotWidget, [=](int value) {
        cerr << "Isosurface simplification currently not supported!";
        // this->tracerVisualiserWidget->manyShow = value;
        // this->tracerVisualiserWidget->recomputeIsosurface = true;
        // this->tracerVisualiserWidget->update();
    });

    connect(fiberSurfaceSimplificationSlider, &QSlider::valueChanged, plotWidget, [=](int value) {
        cerr << "Fiber Surface simplification currently not supported!";
        // this->tracerVisualiserWidget->manyShowF = value;
        // this->tracerVisualiserWidget->recomputeFiberSurface = true;
        // this->tracerVisualiserWidget->update();
    });

    connect(isovalueSlider, &QSlider::valueChanged, tracerVisualiserWidget, [=](float value) {
        // Update the new isovalue
        float t = value / this->isovalueSlider->maximum();
        double isovalue = t * data->isoField->max + (1 - t) * data->isoField->min;

        // Set isovalue
        this->data->isoField->currentIsovalue = isovalue;

        // Loading bar
        this->loadingLabel->setText("LOADING....");
        this->progressBar->setValue(50);
        {
            // Recompute isosurface for all timesteps
            this->data->computeMeshes(isovalue);

            this->tracerVisualiserWidget->generateDisplayList(this->data->isosurfaceMeshes[this->data->currentTimestep],
                                                              SurfaceType::isosurface);

            //this->data->computeCombinedMeshes(this->data->currentSignedDistanceField[this->data->currentTimestep], this->data->isoField->currentIsovalue);
            //this->tracerVisualiserWidget->generateDisplayList(this->data->combinedMeshes[this->data->currentTimestep], SurfaceType::combinedSurface);

            // Recompute the bounding box
            if (-1 != this->data->selectedID[this->data->currentTimestep] &&
                this->data->selectedSurfaceType[this->data->currentTimestep] == SurfaceType::isosurface) {

                this->data->selectedObjectMinMax[this->data->currentTimestep] =
                  this->data->isosurfaceMeshes[this->data->currentTimestep].getMinMaxPointsObject(
                    this->data->selectedID[this->data->currentTimestep]);
            }

            // Update Stuff
            this->tracerVisualiserWidget->update();
            this->plotWidget->update();
        }
        this->loadingLabel->setText("");
        this->progressBar->setValue(0);

        this->isovalueSpinBox->blockSignals(true);
        this->isovalueSpinBox->setValue(isovalue);
        this->isovalueSpinBox->blockSignals(false);
    });

    // This is used to update the value of the spinbox as we move the slider
    connect(isovalueSlider, &QSlider::sliderMoved, tracerVisualiserWidget, [=](float value) {
        float t = value / this->isovalueSlider->maximum();
        double isovalue = t * data->isoField->max + (1 - t) * data->isoField->min;

        this->isovalueSpinBox->blockSignals(true);
        this->isovalueSpinBox->setValue(isovalue);
        this->isovalueSpinBox->blockSignals(false);
    });

    connect(isovalueSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), [=](double isovalue) {
        plotWidget->update();

        // Set isovalue
        this->data->isoField->currentIsovalue = isovalue;

        this->loadingLabel->setText("LOADING....");
        this->progressBar->setValue(50);
        {
            // Recompute the isosurfaces for all timesteps
            this->data->computeMeshes(isovalue);
            this->tracerVisualiserWidget->generateDisplayList(this->data->isosurfaceMeshes[this->data->currentTimestep],
                                                              SurfaceType::isosurface);

            //this->data->computeCombinedMeshes(this->data->currentSignedDistanceField[this->data->currentTimestep], this->data->isoField->currentIsovalue);
            //this->tracerVisualiserWidget->generateDisplayList(this->data->combinedMeshes[this->data->currentTimestep], SurfaceType::combinedSurface);

            // Recompute the bounding box
            if (-1 != this->data->selectedID[this->data->currentTimestep] &&
                this->data->selectedSurfaceType[this->data->currentTimestep] == SurfaceType::isosurface) {
                this->data->selectedObjectMinMax[this->data->currentTimestep] =
                  this->data->isosurfaceMeshes[this->data->currentTimestep].getMinMaxPointsObject(
                    this->data->selectedID[this->data->currentTimestep]);
            }

            // Update Stuff
            this->tracerVisualiserWidget->update();
            this->plotWidget->update();
        }
        this->loadingLabel->setText("");
        this->progressBar->setValue(0);

        // Change the value of the linked slider
        this->isovalueSlider->blockSignals(true);
        {
            float t = (isovalue - data->isoField->min) / (data->isoField->max - data->isoField->min);
            this->isovalueSlider->setValue(static_cast<int>(t * this->isovalueSlider->maximum()));
        }
        this->isovalueSlider->blockSignals(false);
    });

    // connect(saveISObjBtn, &QPushButton::clicked, this, [=]() {
    // data->writeNcData("output-isosurface.nc", tracerVisualiserWidget->visited);
    // //data->writeNcDataVals("output-mesh.nc");
    //});
    // connect(saveFSObjBtn, &QPushButton::clicked, this, [=]() {
    // data->writeNcData("output-fiber-surface.nc",
    // tracerVisualiserWidget->visitedF); });
}

void
TracerVisualiserWindow::setTimestep(int timestep)
{
    assert(timestep >= 0);
    assert(timestep < this->data->tdim);

    this->loadingLabel->setText("LOADING....");
    this->progressBar->setValue(50);
    {
        this->data->changeTimeStep(timestep);

        this->tracerVisualiserWidget->generateDisplayList(this->data->isosurfaceMeshes[this->data->currentTimestep],
                                                          SurfaceType::isosurface);
        this->tracerVisualiserWidget->generateDisplayList(this->data->fibersurfaceMeshes[this->data->currentTimestep],
                                                          SurfaceType::fibersurface);
        this->tracerVisualiserWidget->generateDisplayList(this->data->combinedMeshes[this->data->currentTimestep],
                                                          SurfaceType::combinedSurface);

        this->timeLabel->setText(tr("Timestep: ") + QString::number(this->data->currentTimestep) + tr("/") +
                                 QString::number(this->data->tdim - 1) + " (" +
                                 QString::number(this->data->tVals[this->data->currentTimestep]) + ")");

        // Update visuals
        this->update();
        this->plotWidget->update();
        this->histogramWidget->update();
        this->tracerVisualiserWidget->update();
    }
    this->loadingLabel->setText("");
    this->progressBar->setValue(0);
}

TracerVisualiserWindow::~TracerVisualiserWindow()
{
    delete isovalueRangeLabel;
    delete plotWidget;
    delete tracerVisualiserWidget;
    delete histogramWidget;
    delete timestepLayout;
    delete comboBoxLayout;
    delete isovalueLayout;
    delete histogramLayout;
    delete rightControlsLayout;
    delete leftControlsLayout;
    delete windowLayout;
    delete scatterPlotTypeComboBox;
    delete isoFieldComboBox;
    delete projectionTypeComboBox;
    delete mergeTreeTypeComboBox;
    delete isosurfaceOpacitySlider;
    delete fibersurfaceOpacitySlider;
    delete isosurfaceSimplificationSlider;
    delete fiberSurfaceSimplificationSlider;
    delete timeStepSlider;
    delete isovalueSpinBox;
    delete saveISObjBtn;
    delete saveFSObjBtn;
    delete timeStepPreviousButton;
    delete timeStepNextButton;
    delete scaleSlider;
    delete isovalueSlider;
}
