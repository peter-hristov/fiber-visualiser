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

#include <random>
#include <iomanip>

using std::cout;
using std::endl;
using std::setprecision;

using namespace std;

size_t Data::trippleToIndex(const size_t i, const size_t j, const size_t k)
{
    return this->xdim * this->ydim * k + this->xdim * j + i;
}

void Data::addTetsForCube(const size_t i, const size_t j, const size_t k)
{
    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i, j+1, k), this->trippleToIndex(i+1, j+1, k)});
    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i, j+1, k), this->trippleToIndex(i, j+1, k+1)});
    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i, j, k+1), this->trippleToIndex(i, j+1, k+1)});

    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i+1, j, k), this->trippleToIndex(i+1, j+1, k)});
    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i+1, j, k), this->trippleToIndex(i+1, j, k+1)});
    this->tetrahedra.push_back({this->trippleToIndex(i, j, k), this->trippleToIndex(i+1, j+1, k+1), this->trippleToIndex(i, j, k+1), this->trippleToIndex(i+1, j, k+1)});
}


void
Data::readNcData(tv9k::InputInformation input)
{
    // Read Dimensions
    this->originalXdim = 10;
    this->originalYdim = 10;
    this->originalZdim = 10;
    this->originalTdim = 10;

    this->xdim = 5;
    this->ydim = 5;
    this->zdim = 5;
    this->tdim = 5;

    this->longName = "Variable";

    // Add vertex range coordinates
    //for (int i = 0 ; i < this->xdim * this->ydim * this->zdim ; i++)
    //{
        ////std::random_device rd;
        ////std::default_random_engine eng(rd());
        ////std::uniform_real_distribution<float> distr(0, 10);
        ////this->vertexRangeCoordinates.push_back({distr(eng), distr(eng)});
    //}

    // Add vertex range coordinates
    for (int k = 0 ; k < this->zdim ; k++)
    {
        for (int j = 0 ; j < this->ydim ; j++)
        {
            for (int i = 0 ; i < this->xdim ; i++)
            {
                std::random_device rd;
                std::default_random_engine eng(rd());
                std::uniform_real_distribution<float> distr(0, 0.9);


                const float x = static_cast<float>(i);
                const float y = static_cast<float>(j);
                const float z = static_cast<float>(k);


                //this->vertexRangeCoordinates.push_back({static_cast<float>(i) + 2 + distr(eng), static_cast<float>(j) + 2 + distr(eng)});
                //this->vertexRangeCoordinates.push_back({static_cast<float>(i) + 2 , static_cast<float>(j) + 2});
                //this->vertexRangeCoordinates.push_back({x * x  + y*y + z * z, z});
                this->vertexRangeCoordinates.push_back({(sqrt(x * x  + y*y) - 3) * (sqrt(x * x  + y*y) - 3) + z * z - 2, z});
            }
        }
    }

    this->minF = this->vertexRangeCoordinates[0][0];
    this->maxF = this->vertexRangeCoordinates[0][0];

    this->minG = this->vertexRangeCoordinates[0][1];
    this->maxG = this->vertexRangeCoordinates[0][1];

    for (const auto &val : this->vertexRangeCoordinates)
    {
        this->minF = std::min(this->minF, val[0]);
        this->maxF = std::max(this->maxF, val[0]);

        this->minG = std::min(this->minG, val[1]);
        this->maxG = std::max(this->maxG, val[1]);
    }


    // Add vertex domain coordinates
    for (int k = 0 ; k < this->zdim ; k++)
    {
        for (int j = 0 ; j < this->ydim ; j++)
        {
            for (int i = 0 ; i < this->xdim ; i++)
            {
                this->vertexDomainCoordinates.push_back({static_cast<float>(i), static_cast<float>(j), static_cast<float>(k)});
            }
        }
    }

    // Add tets
    for (int i = 0 ; i < this->xdim - 1 ; i++)
    {
        for (int j = 0 ; j < this->ydim - 1 ; j++)
        {
            for (int k = 0 ; k < this->zdim - 1 ; k++)
            {
                this->addTetsForCube(i, j, k);
            }
        }
    }


    //this->tetrahedra.push_back({i,7,2,3});
    //this->tetrahedra.push_back({i,7,2,6});
    //this->tetrahedra.push_back({i,7,4,6});

    //this->tetrahedra.push_back({i,7,1,3});
    //this->tetrahedra.push_back({i,7,1,5});
    //this->tetrahedra.push_back({i,7,4,5});

    // 1
    //this->vertexDomainCoordinates.push_back({5,0,5});
    //// 2
    //this->vertexDomainCoordinates.push_back({10,5,5});
    //// 3 
    //this->vertexDomainCoordinates.push_back({5,10,5});
    //// 4
    //this->vertexDomainCoordinates.push_back({0,5,5});
    //// 5
    //this->vertexDomainCoordinates.push_back({5,5,10});
    //// 6
    //this->vertexDomainCoordinates.push_back({5,5,0});

    //this->vertexRangeCoordinates.push_back({6,3});
    //this->vertexRangeCoordinates.push_back({2,2});
    //this->vertexRangeCoordinates.push_back({5,5});
    //this->vertexRangeCoordinates.push_back({1,4});
    //this->vertexRangeCoordinates.push_back({3,6});
    //this->vertexRangeCoordinates.push_back({4,1});

    // {1,2,4,5}
    //this->tetrahedra.push_back({0,1,3,4});
    //// {2,3,4,5}
    //this->tetrahedra.push_back({1,2,3,4});
    //// {1,2,4,6}
    //this->tetrahedra.push_back({0,1,3,5});
    //// {2,3,4,6}
    //this->tetrahedra.push_back({1,2,3,5});
}

