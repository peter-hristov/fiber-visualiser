#include "./CGALTypedefs.h"

#include "./Data.h"
#include "./Timer.h"
#include "./FiberPoint.h"
#include "./Fiber.h"

#include <filesystem>
#include <cassert>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <queue>


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
