#pragma once

//#include "Data.h"
#include "src/utility/ScalarField.h"
#include <QImage>

//#define VTK_DIR
//#define TTKBase_DIR

#ifdef VTK_DIR

#include <map>

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "src/utility/ScalarField.h"

#endif

#ifdef TTKBase_DIR

#include <ContinuousScatterPlot.h>
#include <Triangulation.h>

#endif

class ScatterPlot
{
  public:
    ScatterPlot() {}
    ScatterPlot(int);

    void initialize(int);

    //~ScatterPlot() {}

    QImage baseScatterPlotImage;
    QImage scatterPlotImage;

    // Resolution of the scatterPlotImage and grayScatterPlotImage images (both in
    // X and in Y, they are square)
    int resolution = 0;
    // void computeDensity(Data*, bool);
    void populateImage(int);

#ifdef TTKBase_DIR
    // void volumesContinuous(Data*, double);
    // void drawComponent(int);
    // void drawAllComponents();
    // std::map<int, QImage> imageMap;
#endif

    // Same but discrete
    void computeDensityDiscrete(tv9k::utility::ScalarField&, tv9k::utility::ScalarField&, const size_t);

  private:
    // The maximum value in the density array. Used to rescale to color space)
    double maxDensity = 0.0;

    // Array returned by TTK that represent the length of a fiber at points on a
    // grid (the future image)
    std::vector<std::vector<double>> density;

    // Array returned by TTK that represent the length of a fiber at points on a
    // grid (the future image)
    std::vector<std::vector<double>> volumeDensity;

    // Not really using this now, but it is needed for the continuous scatterplot
    std::vector<std::vector<char>> validPointMask;

    // Not really using this now, but it is needed for the continuous scatterplot
    std::vector<std::vector<char>> volumeValidPointMask;

    double computeMaxDensity();

#ifdef TTKBase_DIR
    // Compute the density and pointMask arrays using TTK
    // void computeDensityContinuous(Data*);

    //// Convert 3D cubic point cloud to a triangular mesh
    // vtkSmartPointer<vtkUnstructuredGrid> makeGrid(Data*);

    // std::vector<std::pair<int, vtkSmartPointer<vtkUnstructuredGrid>>>
    // extractSublevelSet(const vtkSmartPointer<vtkUnstructuredGrid>, Data*, const double, const bool);
    // std::vector<vtkSmartPointer<vtkUnstructuredGrid>> grid;
#endif
};
