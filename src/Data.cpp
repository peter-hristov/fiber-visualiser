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


//void Data::printSheetHistogram()
//{
    //Timer::start();
    //if (this->currentFiberPoint.size() == 0)
    //{
        //return;
    //}

    //assert (this->currentFiberPoint.size() == 2);

    //set<int> intersectedSheets;
    //CartesianLine line(CartesianPoint(0, this->currentFiberPoint[0]), CartesianPoint(1, this->currentFiberPoint[1])); // Line through (0, 1) and (1, 0)
    //for (const auto &[sheetId, polygon] : this->reebSpace.sheetPolygon)
    //{
        //QVector<QPoint> points;

        //for (const CartesianSegment& segment : polygon.edges()) 
        //{
            //if (CGAL::do_intersect(segment, line)) 
            //{
                //intersectedSheets.insert(sheetId);
            //}
        //}
    //}

    //std::cout << "Intersected sheets: ";
    //for (const auto &sheetId : intersectedSheets)
    //{
        //std::cout << sheetId << " (a = )" << this->reebSpace.sheetArea[sheetId] << std::endl;

    //}
    //std::cout << "\n";
    //Timer::stop("Computed face intersected by the line  :");
//}
