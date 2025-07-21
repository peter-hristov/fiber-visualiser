#include "./Data.h"
#include "./Timer.h"
#include "./CGALTypedefs.h"
#include "./FiberPoint.h"
#include "./Fiber.h"

#include <CGAL/enum.h>
#include <filesystem>
#include <cassert>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <queue>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <ranges>

#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>

#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

using namespace std;


void Data::saveFibers()
{
    std::cout << "Saving fibers in " << this->fibersFile << std::endl;
    //std::cout << "The fiber has size " << this->faceFibers.size() << std::endl;  

    // 1. Create the points
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto idArray = vtkSmartPointer<vtkIntArray>::New();
    auto colourArray = vtkSmartPointer<vtkFloatArray>::New();

    idArray->SetName("SheetId");
    idArray->SetNumberOfComponents(1);

    colourArray->SetName("Colour");
    colourArray->SetNumberOfComponents(3);

    for (const FiberPoint &p : this->faceFibers)
    {
        // Insert points and corresponding IDs
        points->InsertNextPoint(p.point[0], p.point[1], p.point[2]);
        idArray->InsertNextValue(p.sheetId);

        const array<float, 3> sheetColour = fiber::fiberColours[this->reebSpace.sheetConsequitiveIndices[p.sheetId] % fiber::fiberColours.size()];
        float color[3] = {sheetColour[0], sheetColour[1], sheetColour[2]};
        colourArray->InsertNextTuple(color);
    }

    // 3. Create the cells (wrap polyline in cell array)
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    for (int i = 1 ; i < this->faceFibers.size() ; i+=2)
    {
        //cout << i << " " << this->faceFibers[i-1].sheetId << " " << this->faceFibers[i].sheetId << endl;
        //cout << i << " " << this->faceFibers[i-1].triangleId << " " << this->faceFibers[i].triangleId << endl;
        //printf("(%f, %f, %f) -> (%f, %f, %f)\n", this->faceFibers[i-1].point[0], this->faceFibers[i-1].point[1], this->faceFibers[i-1].point[2], this->faceFibers[i].point[0], this->faceFibers[i].point[1], this->faceFibers[i].point[2]);

        if (this->faceFibers[i-1].sheetId == this->faceFibers[i].sheetId)
        {
            // One edge segment
            auto polyLine = vtkSmartPointer<vtkPolyLine>::New();
            polyLine->GetPointIds()->SetNumberOfIds(2);
            polyLine->GetPointIds()->SetId(0, i-1);
            polyLine->GetPointIds()->SetId(1, i);

            cells->InsertNextCell(polyLine);
        }
    }

    // 4. Create the polydata object
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(cells);

    // 5. Attach the VertexID array to the point data
    polyData->GetPointData()->AddArray(idArray);
    polyData->GetPointData()->AddArray(colourArray);
    polyData->GetPointData()->SetScalars(colourArray);  // optional: for coloring

    // 6. Write to .vtp file (XML format)
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(this->fibersFile.c_str());
    writer->SetInputData(polyData);
    writer->Write();
}

void Data::generatefFaceFibersForSheets(const int sheetOutputCount, const int numberOfFiberPoints, const std::string folderPath)
{
    namespace fs = std::filesystem;

    //string folderPath = "./sheetFibers";
    fs::path folderPathFs(folderPath);
    if (!fs::exists(folderPathFs)) 
    {
        fs::create_directory(folderPathFs);
    }

    for (const auto &[sheetId, colourId] : this->reebSpace.sheetConsequitiveIndices)
    {
        if (this->reebSpace.incompleteSheets.contains(sheetId))
        {
            printf("Skipping fiber %d, it's incomplete.",  sheetId);
        }

        if (colourId > sheetOutputCount || this->reebSpace.incompleteSheets.contains(sheetId))
        {
            continue;
        }

        std::cout << "-------------------------------------------------------------------------------------------- Generating fibers for sheet " << sheetId << "..." << std::endl;
        this->generatefFaceFibersForSheet(sheetId, numberOfFiberPoints);

        //std::cout << "Saving fibers..." << std::endl;
        this->fibersFile = folderPathFs.string() + "/fibers_" + std::to_string(sheetId) + ".vtp";
        this->saveFibers();
        this->faceFibers.clear();
    }

    this->fibersFile = "./fibers.vtp";

}

void Data::generatefFaceFibersForSheet(const int sheetId, const int numberOfFiberPoints)
{
    CartesianPolygon_2 &polygon = this->reebSpace.sheetPolygon[sheetId];

    if (polygon.size() == 0)
    {
        return;
    }


    // Compute the controid so that we can pull all verties towards it
    CartesianPoint centroid = CGAL::centroid(polygon.vertices_begin(), polygon.vertices_end());

    // To make sure we don't write a comma at the end of the array

    vector<vector<float>> fiberPoints;

    // If need only one, get it at the center
    if (numberOfFiberPoints == 1)
    {
        fiberPoints.push_back({(float)centroid.x(), (float)centroid.y()});
        this->computeTetExitPointsNewNew(fiberPoints[0][0], fiberPoints[0][1], false, sheetId);

        return;
    }


    // If we need more, sample along the boundary
    for (const CartesianPoint &point : polygon) 
    {
        // Get point from CGAL (and convert to double )
        float u = point.x();
        float v = point.y();

        // Interpolate closer to the centroid to make sure we are in the sheet ( if the sheet is "convex enough")
        float alpha = 0.2;
        u = (1 - alpha) * u + alpha * centroid.x();
        v = (1 - alpha) * v + alpha * centroid.y();

        fiberPoints.push_back({u, v});
        //fiberPoints.push_back({centroid.x(), centroid.y()});
    }


    // Calculate step size we only want some of the fiber points, not all
    double step = static_cast<double>(fiberPoints.size() - 1) / (numberOfFiberPoints - 1);

    for (int i = 0; i < numberOfFiberPoints; ++i) 
    {
        int index = static_cast<int>(i * step);
        this->computeTetExitPointsNewNew(fiberPoints[index][0], fiberPoints[index][1], false, sheetId);
    }
}


void Data::printSheetHistogram()
{
    Timer::start();
    if (this->currentFiberPoint.size() == 0)
    {
        return;
    }

    assert (this->currentFiberPoint.size() == 2);

    set<int> intersectedSheets;
    CartesianLine line(CartesianPoint(0, this->currentFiberPoint[0]), CartesianPoint(1, this->currentFiberPoint[1])); // Line through (0, 1) and (1, 0)
    for (const auto &[sheetId, polygon] : this->reebSpace.sheetPolygon)
    {
        QVector<QPoint> points;

        for (const CartesianSegment& segment : polygon.edges()) 
        {
            if (CGAL::do_intersect(segment, line)) 
            {
                intersectedSheets.insert(sheetId);
            }
        }
    }

    std::cout << "Intersected sheets: ";
    for (const auto &sheetId : intersectedSheets)
    {
        std::cout << sheetId << " (a = )" << this->reebSpace.sheetArea[sheetId] << std::endl;

    }
    std::cout << "\n";
    Timer::stop("Computed face intersected by the line  :");
}


void Data::computeTetExitPointsNewNew(const GLfloat u, const GLfloat v, const bool clearFibers, const int reebSheetIdOnly)
{
    if (true == clearFibers)
    {
        this->faceFibers.clear();
    }

    const auto newFiber = fiber::computeFiber(tetMesh, arrangement, reebSpace, {u, v}, reebSheetIdOnly); 

    // Append b to a
    this->faceFibers.insert(this->faceFibers.end(), newFiber.begin(), newFiber.end());


    //Timer::stop("Computed fiber in                      :");
}
