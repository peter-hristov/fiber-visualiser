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

// Return a 4D vector
vector<float> inverseStereographicProjection(const vector<float> v)
{
    assert(v.size() == 3);

    float x = v[0];
    float y = v[1];
    float z = v[2];

    float sum = x*x + y*y + z*z;

    return {
        2*x / (1 + sum),
        2*y / (1 + sum),
        2*z / (1 + sum),
        (-1 + sum) / (1 + sum)
    };
}

vector<float> hopfMap(const vector<float> v)
{
    float a = v[0];
    float b = v[1];
    float c = v[2];
    float d = v[3];

    return {
        a*a + b*b - c*c - d*d,
        2 * (a * d + b * c),
        2 * (b * d - a * c)
    };
}


void
Data::readNcData(tv9k::InputInformation input)
{
    // Read Dimensions
    this->originalXdim = 10;
    this->originalYdim = 10;
    this->originalZdim = 10;
    this->originalTdim = 10;

    this->xdim = 2;
    this->ydim = 2;
    this->zdim = 2;
    this->tdim = 5;

    this->minX = 0;
    this->maxX = 2;
    this->minY = 0;
    this->maxY = 2;
    this->minZ = 0;
    this->maxZ = 2;

    this->longnameF = "f = tube size";
    this->longnameG = "g = height";

    // Add vertex domain coordinates
    //for (int k = 0 ; k < this->zdim ; k++)
    //{
        //for (int j = 0 ; j < this->ydim ; j++)
        //{
            //for (int i = 0 ; i < this->xdim ; i++)
            //{
                //this->vertexDomainCoordinates.push_back({static_cast<float>(i), static_cast<float>(j), static_cast<float>(k)});
            //}
        //}
    //}

    // Add vertex domain coordinates
    for (int k = 0 ; k < this->zdim ; k++)
    {
        float ratioZ = static_cast<float>(k)/static_cast<float>(this->zdim);
        for (int j = 0 ; j < this->ydim ; j++)
        {
            float ratioY = static_cast<float>(j)/static_cast<float>(this->ydim);
            for (int i = 0 ; i < this->xdim ; i++)
            {
                float ratioX = static_cast<float>(i)/static_cast<float>(this->xdim);

                // Interpolated the coordinates
                this->vertexDomainCoordinates.push_back({
                        (1 - ratioX) * this->minX + ratioX * this->maxX,
                        (1 - ratioY) * this->minY + ratioY * this->maxY,
                        (1 - ratioZ) * this->minZ + ratioZ * this->maxZ,
                        });
            }
        }
    }



    // Add vertex range coordinates
    for (int i = 0 ; i < this->xdim * this->ydim * this->zdim ; i++)
    {
        //std::random_device rd;
        //std::default_random_engine eng(rd());
        //std::uniform_real_distribution<float> distr(0, 10);
        //this->vertexRangeCoordinates.push_back({distr(eng), distr(eng)});
    }

    // Add vertex range coordinates
    for (int k = 0 ; k < this->zdim ; k++)
    {
        for (int j = 0 ; j < this->ydim ; j++)
        {
            for (int i = 0 ; i < this->xdim ; i++)
            {
                std::random_device rd;
                std::default_random_engine eng(rd());
                std::uniform_real_distribution<float> distr(0, 10);

                const auto coordinates = this->vertexDomainCoordinates[this->trippleToIndex(i, j, k)];

                float x = coordinates[0];
                float y = coordinates[1];
                float z = coordinates[2];

                float fValue = distr(eng); 
                float gValue = distr(eng); 

                this->vertexCoordinatesF.push_back(fValue);
                this->vertexCoordinatesG.push_back(gValue);

                cout << i << " " << j << " " << k << " | f = " << fValue << " g = " << gValue << endl;

                //const auto result = hopfMap(inverseStereographicProjection(coordinates));

                //this->vertexCoordinatesF.push_back(result[1]);
                //this->vertexCoordinatesG.push_back(result[2]);

                //if (x == 0 && y == 0 && z == 0)
                //{
                    //this->vertexCoordinatesF.push_back(0);
                    //this->vertexCoordinatesG.push_back(0);

                //}
                //else
                //{
                    //this->vertexCoordinatesF.push_back(x / sqrt(x*x + y*y + z*z));
                    //this->vertexCoordinatesG.push_back(y / sqrt(x*x + y*y + z*z));
                //}

                //this->vertexCoordinatesG.push_back(result[2]);

                //float x = static_cast<float>(i);
                //float y = static_cast<float>(j);
                //float z = static_cast<float>(k);

                //x -= 5;
                //y -= 5;
                //z -= 5;



                //this->vertexCoordinatesF.push_back(static_cast<float>(i) + 2 + distr(eng));
                //this->vertexCoordinatesG.push_back(static_cast<float>(j) + 2 + distr(eng));

                //this->vertexRangeCoordinates.push_back({static_cast<float>(i) + 2 + distr(eng), static_cast<float>(j) + 2 + distr(eng)});
                //this->vertexRangeCoordinates.push_back({static_cast<float>(i) + 2 , static_cast<float>(j) + 2});
                //this->vertexRangeCoordinates.push_back({x * x  + y*y + z * z, z});
                //this->vertexCoordinatesF.push_back(x + 5);
                //this->vertexCoordinatesF.push_back(x *x + y * y * z * x + z * z );
                //this->vertexCoordinatesG.push_back(x * x - y * y + pow(z, 4));

                //this->vertexCoordinatesF.push_back( x * y * z);
                //this->vertexCoordinatesG.push_back(y * y * x + z * y * z);

                // Torus
                //float torusBigRadius = 20;
                //float torusTubeRadius = 0;
                //this->vertexCoordinatesF.push_back(x + 5);
                //this->vertexCoordinatesG.push_back((sqrt(x*x  + y*y) - torusBigRadius) * (sqrt(x * x  + y*y) - torusBigRadius) + z*z - torusTubeRadius * torusTubeRadius);
            }
        }
    }

    this->minF = this->vertexCoordinatesF[0];
    this->maxF = this->vertexCoordinatesF[0];

    this->minG = this->vertexCoordinatesG[0];
    this->maxG = this->vertexCoordinatesG[0];

    for (int i = 0 ; i < this->vertexCoordinatesF.size() ; i++)
    {
        this->minF = std::min(this->minF, this->vertexCoordinatesF[i]);
        this->maxF = std::max(this->maxF, this->vertexCoordinatesF[i]);

        this->minG = std::min(this->minG, this->vertexCoordinatesG[i]);
        this->maxG = std::max(this->maxG, this->vertexCoordinatesG[i]);
    }

    cout << "The min height is and the max height is " << minG << " " << maxG << endl;

    this->minF -= 2;
    this->maxF += 2;

    this->minG -= 2;
    this->maxG += 2;

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

