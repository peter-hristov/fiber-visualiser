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
}

