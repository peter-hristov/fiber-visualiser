#include "./ScatterPlot.h"

#include "./utility/Utility.h"
#include <QPainter>
#include <cmath>
#include <iostream>

#ifdef VTK_DIR

#include <vtkActor.h>
#include <vtkAppendFilter.h>
#include <vtkCellArray.h>
#include <vtkContinuousScatterplot.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetTriangleFilter.h>
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
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#endif

#ifdef TTKBase_DIR

#include <ContinuousScatterPlot.h>
#include <Triangulation.h>

#endif

using namespace std;

ScatterPlot::ScatterPlot(int _resolution)
{
    this->initialize(_resolution);
}

void
ScatterPlot::initialize(int _resolution)
{
    this->resolution = _resolution;
    this->scatterPlotImage = QImage(resolution, resolution, QImage::Format_ARGB32);
    this->baseScatterPlotImage = QImage(resolution, resolution, QImage::Format_ARGB32);
    density = vector<vector<double>>(resolution, vector<double>(resolution, 0.0));
    validPointMask = vector<vector<char>>(resolution, vector<char>(resolution, 0.0));

    volumeDensity = vector<vector<double>>(resolution, vector<double>(resolution, 0.0));
    volumeValidPointMask = vector<vector<char>>(resolution, vector<char>(resolution, 0.0));
}

double
ScatterPlot::computeMaxDensity()
{
    double maximum = 0.0;

    for (int i = 0; i < resolution; i++) {
        for (int j = 0; j < resolution; j++) {
            maximum = max(maximum, density[i][j]);
        }
    }

    return maximum;
}

// void
// ScatterPlot::computeDensity(Data* data, bool continuous)
//{
// if (continuous) {
//#ifdef VTK_DIR
// this->computeDensityContinuous(data);
//#else
// cout << "Continuous ScatterPlot not supported." << endl;
//#endif
//} else {
// this->computeDensityDiscrete(data);
//}
//}

void
ScatterPlot::computeDensityDiscrete(tv9k::utility::ScalarField& uField,
                                    tv9k::utility::ScalarField& vField,
                                    const size_t timestep)
{
    // Make sure the two fields have the same dimensions
    assert(uField.xDim() == vField.xDim());
    assert(uField.yDim() == vField.yDim());
    assert(uField.zDim() == vField.zDim());
    assert(uField.tDim() == vField.tDim());

    density = vector<vector<double>>(resolution, vector<double>(resolution, 0.0));
    validPointMask = vector<vector<char>>(resolution, vector<char>(resolution, 0.0));

    for (int i = 0; i < uField.xDim(); i++) {
        for (int j = 0; j < uField.yDim(); j++) {
            for (int k = 0; k < uField.zDim(); k++) {

                float xFloat =
                  (resolution / (uField.max - uField.min)) * (uField.values[timestep][i][j][k] - uField.min);

                float yFloat =
                  (resolution / (vField.max - vField.min)) * (vField.values[timestep][i][j][k] - vField.min);

                int x = static_cast<int>(xFloat);
                int y = static_cast<int>(yFloat);

                // Clip do domain boundaries
                x = max(0, x);
                x = min(x, resolution - 1);

                y = max(0, y);
                y = min(y, resolution - 1);

                density[x][y] += 1.0;
            }
        }
    }

    maxDensity = computeMaxDensity();
}

void
ScatterPlot::populateImage(int factor)
{
    this->scatterPlotImage = QImage(resolution, resolution, QImage::Format_ARGB32);
    this->baseScatterPlotImage = QImage(resolution, resolution, QImage::Format_ARGB32);

    for (int i = 0; i < density.size(); i++) {
        for (int j = 0; j < density[i].size(); j++) {
            if (factor == 0) {
                if (density[i][j] > 0.0) {
                    baseScatterPlotImage.setPixelColor(QPoint(i, j), Qt::gray);
                } else {
                    baseScatterPlotImage.setPixelColor(QPoint(i, j), QColor(255, 255, 255));
                }
            } else {
                QColor color;

                // double a = min(factor * density[i][j]/maxDensity, 1.0);
                double a = min(log(factor * density[i][j]) / log(maxDensity), 1.0);
                // double a = min(sqrt(factor * density[i][j])/sqrt(maxDensity), 1.0);
                a = max(0.0, a);

                // if (a < 0.00001)
                //{
                // baseScatterPlotImage.setPixelColor(QPoint(i, j), QColor(255, 255,
                // 255));
                //}
                // else
                {
                    double c = 1 - a * (1 - 0.5);

                    color.setRedF(c);
                    color.setGreenF(c);
                    color.setBlueF(c);

                    baseScatterPlotImage.setPixelColor(QPoint(i, j), color);
                }

                // White -> Blue
                // if (a <= 0.25)
                //{
                // a = a / 0.25;

                //// White to Blue
                // color.setRedF(1.0 - a);
                // color.setGreenF(1.0 - a);
                // color.setBlueF(1.0);
                //}
                //// Blue -> Green
                // else if (0.25 < a && a <= 0.5)
                //{
                // a = (a - 0.25) / 0.25;

                //// Blue to Red
                // color.setRedF(0);
                // color.setGreenF(a);
                // color.setBlueF(1.0 - a);
                //}
                //// Green -> Yellow
                // else if (0.5 < a && a <= 0.75)
                //{
                // a = (a - 0.5) / 0.25;

                //// Blue to Red
                // color.setRedF(a);
                // color.setGreenF(1.0);
                // color.setBlueF(0.0);
                //}
                //// Yellow -> Red
                // else
                //{
                // a = (a - 0.75) / 0.25;

                //// Blue to Red
                // color.setRedF(1.0);
                // color.setGreenF(1.0 - a);
                // color.setBlueF(0.0);
                //}

                // baseScatterPlotImage.setPixelColor(QPoint(i, j), color);
            }
        }
    }
    scatterPlotImage = baseScatterPlotImage;
}

//#ifdef TTKBase_DIR

// void
// ScatterPlot::drawAllComponents()
//{
// scatterPlotImage = baseScatterPlotImage;

// QPainter p(&scatterPlotImage);

// for (auto i : imageMap) {
// p.drawImage(0, 0, i.second);
// i.second.save("./features/image_vol_" + QString::number(i.first) + ".png");
//}
//}

// void
// ScatterPlot::drawComponent(int component)
//{
// scatterPlotImage = baseScatterPlotImage;

// QPainter p(&scatterPlotImage);
// p.drawImage(0, 0, imageMap[component]);
//}

// void
// ScatterPlot::volumesContinuous(Data* data, double isovalue)
//{
// Utility::startTimer();
// printf("\n\n------------------------------------ Computing Sublevel set "
//"triangulation. \n");

// vector<pair<int, vtkSmartPointer<vtkUnstructuredGrid>>> grid =
// this->extractSublevelSet(this->makeGrid(data), data, isovalue, false);

// printf("------------------------------------ Done - %3ld.%06ld seconds.\n",
// Utility::endTimer().first,
// Utility::endTimer().second);

// cout << "THE SIZE IS " << grid.size();

// for (int t = 0; t < grid.size(); t++) {
// imageMap[grid[t].first] = QImage(resolution, resolution, QImage::Format_ARGB32);
//}

// for (int t = 0; t < grid.size(); t++) {
// vtkSmartPointer<vtkUnstructuredGrid> ugrid = grid[t].second;

// if (ugrid->GetNumberOfCells() < 1) {
// continue;
//}

// volumeDensity = vector<vector<double>>(resolution, vector<double>(resolution, 0.0));

// Utility::startTimer();
// printf("\n\n------------------------------------ Computing Projection %d. \n", t);

// ttk::Triangulation* triangulation = new ttk::Triangulation;
// triangulation->setInputPoints(ugrid->GetNumberOfPoints(),
//((vtkUnstructuredGrid*)ugrid)->GetPoints()->GetVoidPointer(0));
// triangulation->setInputCells(ugrid->GetNumberOfCells(), ugrid->GetCells()->GetPointer());

// ttk::ContinuousScatterPlot continuousScatterPlot;
// continuousScatterPlot.setVertexNumber(data->xdim * data->ydim * data->zdim);
// continuousScatterPlot.setDummyValue(true, 0);
// continuousScatterPlot.setTriangulation(triangulation);
// continuousScatterPlot.setResolutions(resolution, resolution);

// continuousScatterPlot.setInputScalarField1(ugrid->GetPointData()->GetArray("f1")->GetVoidPointer(0));
// continuousScatterPlot.setInputScalarField2(ugrid->GetPointData()->GetArray("f2")->GetVoidPointer(0));

// double scalarMin[2] = { static_cast<double>(data->uField->min), static_cast<double>(data->vField->min) };
// double scalarMax[2] = { static_cast<double>(data->uField->max), static_cast<double>(data->vField->max) };

// continuousScatterPlot.setScalarMin(scalarMin);
// continuousScatterPlot.setScalarMax(scalarMax);
// continuousScatterPlot.setOutputDensity(&volumeDensity);
// continuousScatterPlot.setOutputMask(&volumeValidPointMask);
// int returnCode = continuousScatterPlot.execute<double, double>();

// printf("------------------------------------ Done - %3ld.%06ld seconds.\n",
// Utility::endTimer().first,
// Utility::endTimer().second);

// for (int i = 0; i < volumeDensity.size(); i++) {
// for (int j = 0; j < volumeDensity[i].size(); j++) {
// if (volumeDensity[i][j] > 0.0) {
// auto color = Utility::getColorQt(grid[t].first);
//// color.setAlpha(100);
// imageMap[grid[t].first].setPixelColor(QPoint(i, j), color);

//// scatterPlotImage.setPixelColor(QPoint(i, j),
//// Utility::getColorQt(grid[t].first));
//// scatterPlotImage.setPixelColor(QPoint(i, j), Qt::blue);
//} else {
// imageMap[grid[t].first].setPixelColor(QPoint(i, j), QColor(255, 255, 255, 0));
//}
//}
//}
//}
// drawAllComponents();
//}

// void
// ScatterPlot::computeDensityContinuous(Data* data)
//{
// density = vector<vector<double>>(resolution, vector<double>(resolution, 0.0));
// validPointMask = vector<vector<char>>(resolution, vector<char>(resolution, 0.0));

// vtkSmartPointer<vtkUnstructuredGrid> ugrid = this->makeGrid(data);
//// vtkSmartPointer<vtkUnstructuredGrid> ugrid =
//// this->extractSublevelSet(this->makeGrid(data), data, -0.5, false)[0];

//// vtkSmartPointer<vtkUnstructuredGrid> ugrid =
//// this->extractSublevelSet(ugrid, 0.5, false);

//// vtkSmartPointer<vtkContinuousScatterplot> csp =
//// vtkSmartPointer<vtkContinuousScatterplot>::New(); csp->SetInputData(ugrid);
//// csp->SetField1("f1",resolution);
//// csp->SetField2("f2",resolution);
//// csp->Update();

//// for (int z = 0; z < 1; z++)
////{
//// for (int y = 0; y < resolution; y++)
////{
//// for (int x = 0; x < resolution; x++)
////{
//// double* pixel = static_cast<double*>(csp->GetOutput()->GetScalarPointer(x,
//// y, z));
//////cout<< (int)(pixel[0] + 0.5) << " , ";
//// density[x][y] = pixel[0] + 0.5;
////}
////}
////}

//// vtkSmartPointer<vtkXMLImageDataWriter> writer =
//// vtkSmartPointer<vtkXMLImageDataWriter>::New();
//// writer->SetFileName("image.vti");
//// writer->SetInputData(csp->GetOutput());
//// writer->Write();

//// return;

////
//// ret = getTriangulation(input);
////
// ttk::Triangulation* triangulation = new ttk::Triangulation;
//// Why does this not work?
//// triangulation->setInputPoints(ugrid->GetNumberOfPoints(),
//// ugrid->GetPoints()->GetVoidPointer(0), true);
//// @TODO WTF why do you even need the cast here?
// triangulation->setInputPoints(ugrid->GetNumberOfPoints(),
//((vtkUnstructuredGrid*)ugrid)->GetPoints()->GetVoidPointer(0));
// triangulation->setInputCells(ugrid->GetNumberOfCells(), ugrid->GetCells()->GetPointer());

////
//// Alternative way to set input fields
//// ret = getScalars(input);
////
//// vtkDataArray* inputScalars1_;
//// vtkDataArray* inputScalars2_;
//// vtkPointData* pointData = ugrid->GetPointData();
//// inputScalars1_ = pointData->GetArray("f1");
//// inputScalars2_ = pointData->GetArray("f2");
//// continuousScatterPlot.setInputScalarField1(inputScalars1_->GetVoidPointer(0));
//// continuousScatterPlot.setInputScalarField2(inputScalars2_->GetVoidPointer(0));

// ttk::ContinuousScatterPlot continuousScatterPlot;
//// continuousScatterPlot.setWrapper(this);
// continuousScatterPlot.setVertexNumber(data->xdim * data->ydim * data->zdim);
// continuousScatterPlot.setDummyValue(true, 0);
// continuousScatterPlot.setTriangulation(triangulation);
// continuousScatterPlot.setResolutions(resolution, resolution);

// continuousScatterPlot.setInputScalarField1(ugrid->GetPointData()->GetArray("f1")->GetVoidPointer(0));
// continuousScatterPlot.setInputScalarField2(ugrid->GetPointData()->GetArray("f2")->GetVoidPointer(0));

// double scalarMin[2] = { static_cast<double>(data->uField->min), static_cast<double>(data->vField->min) };
// double scalarMax[2] = { static_cast<double>(data->uField->max), static_cast<double>(data->vField->max) };

// continuousScatterPlot.setScalarMin(scalarMin);
// continuousScatterPlot.setScalarMax(scalarMax);
// continuousScatterPlot.setOutputDensity(&density);
// continuousScatterPlot.setOutputMask(&validPointMask);
// int returnCode = continuousScatterPlot.execute<double, double>();

// maxDensity = computeMaxDensity();
//}

// vtkSmartPointer<vtkUnstructuredGrid>
// ScatterPlot::makeGrid(Data* data)
//{
// vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
// vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

// vtkSmartPointer<vtkDoubleArray> isoWeights = vtkSmartPointer<vtkDoubleArray>::New();
// isoWeights->SetNumberOfValues(data->xdim * data->ydim * data->zdim * 2);
// isoWeights->SetName("i");

// vtkSmartPointer<vtkDoubleArray> weights = vtkSmartPointer<vtkDoubleArray>::New();
// weights->SetNumberOfValues(data->xdim * data->ydim * data->zdim * 2);
// weights->SetName("f1");

// vtkSmartPointer<vtkDoubleArray> weights2 = vtkSmartPointer<vtkDoubleArray>::New();
// weights2->SetNumberOfValues(data->xdim * data->ydim * data->zdim * 2);
// weights2->SetName("f2");

// for (int j = 0; j < data->ydim; j++) {
// for (int i = 0; i < data->xdim; i++) {
// for (int k = 0; k < data->zdim; k++) {
// const double point[3] = { i, j, k };
// points->InsertPoint(data->getIndex(i, j, k), point);

// isoWeights->SetValue(data->getIndex(i, j, k), static_cast<double>(data->vals[i][j][k]));
// weights->SetValue(data->getIndex(i, j, k), static_cast<double>(data->valsU()[i][j][k]));
// weights2->SetValue(data->getIndex(i, j, k), static_cast<double>(data->valsV()[i][j][k]));
//}
//}
//}

// for (int j = 0; j < data->ydim - 1; j++) {
// for (int i = 0; i < data->xdim - 1; i++) {
// for (int k = 0; k < data->zdim - 1; k++) {
//// vtkIdType cube[8] = {
//// data->getIndex(i, j, k),
//// data->getIndex(i, j, k + 1),
//// data->getIndex(i, j + 1, k),
//// data->getIndex(i, j + 1, k + 1),
//// data->getIndex(i + 1, j, k),
//// data->getIndex(i + 1, j, k + 1),
//// data->getIndex(i + 1, j + 1, k),
//// data->getIndex(i + 1, j + 1, k + 1),
////};

//// ugrid->InsertNextCell(VTK_HEXAHEDRON, 8, cube);

// vtkIdType tets[6][4] = {

//{ data->getIndex(i, j, k),
// data->getIndex(i + 1, j + 1, k + 1),
// data->getIndex(i, j, k + 1),
// data->getIndex(i + 1, j, k + 1) },
//{ data->getIndex(i, j, k),
// data->getIndex(i + 1, j + 1, k + 1),
// data->getIndex(i + 1, j, k),
// data->getIndex(i + 1, j, k + 1) },
//{ data->getIndex(i, j, k),
// data->getIndex(i + 1, j + 1, k + 1),
// data->getIndex(i + 1, j, k),
// data->getIndex(i + 1, j + 1, k) },
//{ data->getIndex(i, j, k),
// data->getIndex(i + 1, j + 1, k + 1),
// data->getIndex(i, j, k + 1),
// data->getIndex(i, j + 1, k + 1) },
//{ data->getIndex(i, j, k),
// data->getIndex(i + 1, j + 1, k + 1),
// data->getIndex(i, j + 1, k),
// data->getIndex(i, j + 1, k + 1) },
//{ data->getIndex(i, j, k),
// data->getIndex(i + 1, j + 1, k + 1),
// data->getIndex(i, j + 1, k),
// data->getIndex(i + 1, j + 1, k) }
//};

// for (int w = 0; w < 6; w++) {
// ugrid->InsertNextCell(VTK_TETRA, 4, tets[w]);
//}
//}
//}
//}

// ugrid->SetPoints(points);
// ugrid->GetPointData()->AddArray(isoWeights);
// ugrid->GetPointData()->AddArray(weights);
// ugrid->GetPointData()->AddArray(weights2);

//// auto tetra  = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
//// tetra->SetInputData(ugrid);
//// tetra->Update();
//// auto grid =  tetra->GetOutput();

//// return this->extractSublevelSet(ugrid, 0.5, false);

// return ugrid;
//}

// vector<pair<int, vtkSmartPointer<vtkUnstructuredGrid>>>
// ScatterPlot::extractSublevelSet(const vtkSmartPointer<vtkUnstructuredGrid> ugrid,
// Data* data,
// const double iso,
// const bool debug)
//{
// auto visitedArray = data->tree->computeVisitedForIsovalue(iso);

// vtkSmartPointer<vtkIntArray> visitedVtk = vtkSmartPointer<vtkIntArray>::New();
// visitedVtk->SetNumberOfValues(data->xdim * data->ydim * data->zdim * 2);
// visitedVtk->SetName("visited");

// for (int j = 0; j < data->ydim; j++) {
// for (int i = 0; i < data->xdim; i++) {
// for (int k = 0; k < data->zdim; k++) {
// visitedVtk->SetValue(data->getIndex(i, j, k), static_cast<int>(visitedArray[i][j][k]));
//}
//}
//}

// map<int, vtkSmartPointer<vtkAppendFilter>> appendMap;

//// For Tets
// for (vtkIdType tetraIndex = 0; tetraIndex < ugrid->GetNumberOfCells(); tetraIndex++) {
// int onePercent = 10 * ugrid->GetNumberOfCells() / 100;

// if (tetraIndex % onePercent == 0) {
// printf("Processing (%.2f%) tet %d from %d.\n",
// static_cast<double>(100.0 * tetraIndex / ugrid->GetNumberOfCells()),
// tetraIndex,
// ugrid->GetNumberOfCells());
//}

// int vis = -1;

//// cell is the current Tetrahedron
// vtkSmartPointer<vtkIdList> cell = vtkSmartPointer<vtkIdList>::New();
// ugrid->GetCellPoints(tetraIndex, cell);

// if (debug) {
// cout << endl << endl << endl << "Tet : " << endl;
//}

// vtkSmartPointer<vtkIdList> black = vtkSmartPointer<vtkIdList>::New();
// vtkSmartPointer<vtkIdList> white = vtkSmartPointer<vtkIdList>::New();

// vector<pair<double, double>> blackFiber;
//// vector<pair<double, double>> whiteFiber;
// vector<pair<double, double>> crossFiber;

//// For Vertices in a Tet
// for (vtkIdType cellIndex = 0; cellIndex < cell->GetNumberOfIds(); cellIndex++) {
// vtkIdType pointId = cell->GetId(cellIndex);

// double isovalue = ugrid->GetPointData()->GetArray("i")->GetComponent(pointId, 0);
//// Value at the vertices
// double value = ugrid->GetPointData()->GetArray("f1")->GetComponent(pointId, 0);
// double value2 = ugrid->GetPointData()->GetArray("f2")->GetComponent(pointId, 0);

// double point[3];
// ugrid->GetPoints()->GetPoint(pointId, point);
//// points->GetPoint(pointId, point);

// if (isovalue < iso) {
// vis = visitedVtk->GetComponent(pointId, 0);

// black->InsertNextId(pointId);
// blackFiber.push_back(make_pair(value, value2));
// if (debug) {
// printf("Black - ");
//}
//} else {
// white->InsertNextId(pointId);
//// whiteFiber.push_back(make_pair(value, value2));
// if (debug) {
// printf("White - ");
//}
//}

// if (debug) {
// printf("Id = %d, Coordinates = (%.1f, %.1f, %.1f), Values = (%.4f | "
//"%.4f, %.4f)\n",
// pointId,
// point[0],
// point[1],
// point[2],
// isovalue,
// value,
// value2);
//}

//// delete point;
//}

////
//// Find the crosed edges
////
// vtkSmartPointer<vtkPoints> crossPoints = vtkSmartPointer<vtkPoints>::New();
// for (vtkIdType b = 0; b < black->GetNumberOfIds(); b++) {
// for (vtkIdType w = 0; w < white->GetNumberOfIds(); w++) {
// long long int pointId = black->GetId(b);
// long long int pointId2 = white->GetId(w);

// double point[3];
// ugrid->GetPoints()->GetPoint(pointId, point);

// double point2[3];
// ugrid->GetPoints()->GetPoint(pointId2, point2);

// double isoval = ugrid->GetPointData()->GetArray("i")->GetComponent(pointId, 0);
// double isoval2 = ugrid->GetPointData()->GetArray("i")->GetComponent(pointId2, 0);

// double valueF = ugrid->GetPointData()->GetArray("f1")->GetComponent(pointId, 0);
// double valueF2 = ugrid->GetPointData()->GetArray("f2")->GetComponent(pointId, 0);

// double value2F = ugrid->GetPointData()->GetArray("f1")->GetComponent(pointId2, 0);
// double value2F2 = ugrid->GetPointData()->GetArray("f2")->GetComponent(pointId2, 0);

// double t = (iso - isoval) / (isoval2 - isoval);

//// Value at the cross points
// double f1 = value2F * t + (1 - t) * valueF;
// double f2 = value2F2 * t + (1 - t) * valueF2;

// crossFiber.push_back(make_pair(f1, f2));

// if (debug) {
// cout << "t is " << t << " ";
//}

// double crossPoint[3] = {
// point2[0] * t + (1 - t) * point[0],
// point2[1] * t + (1 - t) * point[1],
// point2[2] * t + (1 - t) * point[2],
//};

// crossPoints->InsertNextPoint(crossPoint);

// if (debug) {
// printf("Cross edge - %d %d (%.1f, %.1f, %.1f, (%f, %f, %f)) -> "
//"(%.1f, %.1f, %.1f, (%f, %f, %f) at (%.4f, %.4f, %.4f) fiber "
//"inter = (%f, %f)\n",
// pointId,
// pointId2,
// point[0],
// point[1],
// point[2],
// isoval,
// valueF,
// valueF2,
// point2[0],
// point2[1],
// point2[2],
// isoval2,
// value2F,
// value2F2,
// crossPoint[0],
// crossPoint[1],
// crossPoint[2],
// f1,
// f2);
//}
//}
//}

// if (debug) {
// cout << "Cross points are : " << endl;
//}
// for (vtkIdType p = 0; p < crossPoints->GetNumberOfPoints(); p++) {
// double point[3];
// crossPoints->GetPoint(p, point);

// if (debug) {
// printf("%d (%.4f, %.4f, %.4f)\n", p, point[0], point[1], point[2]);
//}
//}

// if (debug) {
// cout << endl << endl;
//}

// vtkSmartPointer<vtkPoints> blackPolyPoints = vtkSmartPointer<vtkPoints>::New();
//// vtkSmartPointer<vtkPoints> whitePolyPoints =
//// vtkSmartPointer<vtkPoints>::New();

// for (vtkIdType b = 0; b < black->GetNumberOfIds(); b++) {
// long long int pointId = black->GetId(b);
// double point[3];
// ugrid->GetPoints()->GetPoint(pointId, point);

// blackPolyPoints->InsertNextPoint(point);
//}

//// for (vtkIdType b = 0; b < white->GetNumberOfIds(); b++)
////{
//// long long int pointId = white->GetId(b);
//// double point[3];
//// ugrid->GetPoints()->GetPoint(pointId, point);

//// whitePolyPoints->InsertNextPoint(point);
////}

// for (vtkIdType p = 0; p < crossPoints->GetNumberOfPoints(); p++) {
// double point[3];
// crossPoints->GetPoint(p, point);

// blackPolyPoints->InsertNextPoint(point);
//// whitePolyPoints->InsertNextPoint(point);
//}

// if (debug) {
// cout << "Black Polyhedron points are : " << endl;
//}
// for (vtkIdType p = 0; p < blackPolyPoints->GetNumberOfPoints(); p++) {
// double point[3];
// blackPolyPoints->GetPoint(p, point);

// if (debug) {
// printf("%d (%.4f, %.4f, %.4f)\n", p, point[0], point[1], point[2]);
//}
//}

//// if (debug)
////{
//// cout << endl << "White Polyhedron points are : " << endl;
////}
//// for (vtkIdType p = 0; p < whitePolyPoints->GetNumberOfPoints(); p++)
////{
//// double point[3];
//// whitePolyPoints->GetPoint(p, point);

//// if (debug)
////{
//// printf("%d (%.4f, %.4f, %.4f)\n", p, point[0], point[1], point[2]);
////}
////}

//// vector<pair<double, double>> blackFiber;

//// It goes black, than cross

// vtkSmartPointer<vtkDoubleArray> f1 = vtkSmartPointer<vtkDoubleArray>::New();
// f1->SetNumberOfValues(blackPolyPoints->GetNumberOfPoints());
// f1->SetName("f1");

// vtkSmartPointer<vtkDoubleArray> f2 = vtkSmartPointer<vtkDoubleArray>::New();
// f2->SetNumberOfValues(blackPolyPoints->GetNumberOfPoints());
// f2->SetName("f2");

// for (vtkIdType i = 0; i < black->GetNumberOfIds(); i++) {
// f1->SetValue(i, blackFiber[i].first);
// f2->SetValue(i, blackFiber[i].second);

// if (debug) {
// cout << "Writing in " << i << " " << blackFiber[i].first << " " << blackFiber[i].second << endl;
//}
//}

// for (vtkIdType i = 0; i < crossPoints->GetNumberOfPoints(); i++) {
// f1->SetValue(i + black->GetNumberOfIds(), crossFiber[i].first);
// f2->SetValue(i + black->GetNumberOfIds(), crossFiber[i].second);

// if (debug) {
// cout << "Writing in " << i + black->GetNumberOfIds() << " " << crossFiber[i].first << " "
//<< crossFiber[i].second << endl;
//}
//}

// vtkSmartPointer<vtkUnstructuredGrid> blackGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
//// vtkSmartPointer<vtkUnstructuredGrid> whiteGrid =
//// vtkSmartPointer<vtkUnstructuredGrid>::New();

// blackGrid->SetPoints(blackPolyPoints);
//// whiteGrid->SetPoints(whitePolyPoints);

// blackGrid->GetPointData()->AddArray(f1);
// blackGrid->GetPointData()->AddArray(f2);

//// All White, disperse
// if (black->GetNumberOfIds() == 0) {
//}
// if (black->GetNumberOfIds() == 1) {
// vtkIdType tet[4] = { 0, 1, 2, 3 };

// blackGrid->InsertNextCell(VTK_TETRA, 4, tet);

//// vtkIdType tets[3][4] = {
////{0, 3, 4, 5},
////{0, 1, 2, 4},
////{0, 2, 4, 5}
////};

//// whiteGrid->InsertNextCell(VTK_TETRA, 4, tets[0]);
//// whiteGrid->InsertNextCell(VTK_TETRA, 4, tets[1]);
//// whiteGrid->InsertNextCell(VTK_TETRA, 4, tets[2]);
//} else if (black->GetNumberOfIds() == 3) {
//// vtkIdType tet[4] = {0, 1, 2, 3};

//// whiteGrid->InsertNextCell(VTK_TETRA, 4, tet);

// vtkIdType tets[3][4] = { { 0, 3, 4, 5 }, { 0, 1, 2, 4 }, { 0, 2, 4, 5 } };

// blackGrid->InsertNextCell(VTK_TETRA, 4, tets[0]);
// blackGrid->InsertNextCell(VTK_TETRA, 4, tets[1]);
// blackGrid->InsertNextCell(VTK_TETRA, 4, tets[2]);
//} else if (black->GetNumberOfIds() == 2) {
//// vtkIdType whiteTets[3][4] = {
////{0, 3, 4, 5},
////{0, 3, 4, 2},
////{0, 1, 3, 5}
////};

//// whiteGrid->InsertNextCell(VTK_TETRA, 4, whiteTets[0]);
//// whiteGrid->InsertNextCell(VTK_TETRA, 4, whiteTets[1]);
//// whiteGrid->InsertNextCell(VTK_TETRA, 4, whiteTets[2]);

// vtkIdType tets[3][4] = { { 0, 1, 4, 5 }, { 0, 2, 4, 5 }, { 0, 2, 3, 5 } };

// blackGrid->InsertNextCell(VTK_TETRA, 4, tets[0]);
// blackGrid->InsertNextCell(VTK_TETRA, 4, tets[1]);
// blackGrid->InsertNextCell(VTK_TETRA, 4, tets[2]);
//} else if (black->GetNumberOfIds() == 4) {
// blackGrid->SetPoints(blackPolyPoints);
// vtkIdType tet[4] = { 0, 1, 2, 3 };
// blackGrid->InsertNextCell(VTK_TETRA, 4, tet);
//}

// if (vis != -1) {
// if (appendMap.find(vis) == appendMap.end()) {
// cout << "New component " << vis << endl;
// appendMap.insert(make_pair(vis, vtkSmartPointer<vtkAppendFilter>::New()));
//}

// appendMap[vis]->AddInputData(blackGrid);
//}
//}

// vector<pair<int, vtkSmartPointer<vtkUnstructuredGrid>>> sublevelSets;

// for (auto i : appendMap) {
// i.second->Update();
// sublevelSets.push_back(make_pair(i.first, i.second->GetOutput()));
//}

//// for(int i = 0 ; i < active.size() ; i++)
////{
//// append[i]->Update();
//// sublevelSets.push_back(append[i]->GetOutput());
////}

// return sublevelSets;
//}

//#endif
