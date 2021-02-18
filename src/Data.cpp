#include "Data.h"
#include "./utility/Geometry.h"

#include <cassert>
#include <cstdio>
#include <omp.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_point.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

using namespace std;


void
Data::readNcData(tv9k::InputInformation input)
{
    // Read Dimensions
    this->originalXdim = 10;
    this->originalYdim = 10;
    this->originalZdim = 10;
    this->originalTdim = 10;

    this->xdim = 10;
    this->ydim = 10;
    this->zdim = 10;
    this->tdim = 10;

    this->min = 0;
    this->max = 10;
    this->longName = "Variable";

    // 1
    this->vertexDomainCoordinates.push_back({5,0,5});
    // 2
    this->vertexDomainCoordinates.push_back({10,5,5});
    // 3 
    this->vertexDomainCoordinates.push_back({5,10,5});
    // 4
    this->vertexDomainCoordinates.push_back({0,5,5});
    // 5
    this->vertexDomainCoordinates.push_back({5,5,10});
    // 6
    this->vertexDomainCoordinates.push_back({5,5,0});

    // 1
    this->vertexRangeCoordinates.push_back({6,3});
    // 2
    this->vertexRangeCoordinates.push_back({2,2});
    // 3
    this->vertexRangeCoordinates.push_back({5,5});
    // 4
    this->vertexRangeCoordinates.push_back({1,4});
    // 5
    this->vertexRangeCoordinates.push_back({3,6});
    // 6
    this->vertexRangeCoordinates.push_back({4,1});

    // {1,2,4,5}
    this->tetrahedra.push_back({0,1,3,4});
    // {2,3,4,5}
    this->tetrahedra.push_back({1,2,3,4});
    // {1,2,4,6}
    this->tetrahedra.push_back({0,1,3,5});
    // {2,3,4,6}
    this->tetrahedra.push_back({1,2,3,5});
}

