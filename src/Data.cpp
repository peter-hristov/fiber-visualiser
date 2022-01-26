#include "Data.h"
#include "./utility/Geometry.h"
#include "./utility/utility.h"

#include <cassert>
#include <cstdio>
#include <omp.h>

//#include <boost/geometry.hpp>
//#include <boost/geometry/geometries/linestring.hpp>
//#include <boost/geometry/geometries/point_xy.hpp>
//#include <boost/geometry/geometries/polygon.hpp>
//#include <boost/geometry/geometries/multi_point.hpp>
//#include <boost/geometry/geometries/multi_polygon.hpp>

#include <random>
#include <iomanip>
#include <utility>

using std::cout;
using std::endl;

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

void Data::generateStandardGrid()
{
    // Read Dimensions
    this->originalXdim = 10;
    this->originalYdim = 10;
    this->originalZdim = 10;
    this->originalTdim = 10;

    this->xdim = 4;
    this->ydim = 4;
    this->zdim = 4;
    this->tdim = 5;

    this->minX = 0;
    this->maxX = this->xdim;
    this->minY = 0;
    this->maxY = this->ydim;
    this->minZ = 0;
    this->maxZ = this->zdim;

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

                //cout << i << " " << j << " " << k << " | f = " << fValue << " g = " << gValue << endl;

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

    //this->vertexCoordinatesF[0] = 9;
    //this->vertexCoordinatesG[0] = 8;
    //this->vertexCoordinatesF[1] = 1;
    //this->vertexCoordinatesG[1] = 4;
    //this->vertexCoordinatesF[2] = 5;
    //this->vertexCoordinatesG[2] = 3;
    //this->vertexCoordinatesF[3] = 3;
    //this->vertexCoordinatesG[3] = 3;
    //this->vertexCoordinatesF[4] = 9;
    //this->vertexCoordinatesG[4] = 3;
    //this->vertexCoordinatesF[5] = 0;
    //this->vertexCoordinatesG[5] = 2;
    //this->vertexCoordinatesF[6] = 6.5;
    //this->vertexCoordinatesG[6] = 3.5;
    //this->vertexCoordinatesF[7] = 0;
    //this->vertexCoordinatesG[7] = 8;

    this->vertexCoordinatesF[0] = 9;
    this->vertexCoordinatesG[0] = 8;
    this->vertexCoordinatesF[1] = 1;
    this->vertexCoordinatesG[1] = 4;
    this->vertexCoordinatesF[2] = 4;
    this->vertexCoordinatesG[2] = 9;
    this->vertexCoordinatesF[3] = 3;
    this->vertexCoordinatesG[3] = 3;

    this->vertexCoordinatesF[4] = 7;
    this->vertexCoordinatesG[4] = 9;

    this->vertexCoordinatesF[5] = 0;
    this->vertexCoordinatesG[5] = 2;
    this->vertexCoordinatesF[6] = 6.5;
    this->vertexCoordinatesG[6] = 3.5;
    this->vertexCoordinatesF[7] = 0;
    this->vertexCoordinatesG[7] = 8;


    // 8 1st exit
    //this->vertexDomainCoordinates.push_back({0,2,1});
    //this->vertexCoordinatesF.push_back(0);
    //this->vertexCoordinatesG.push_back(10);

    //this->tetrahedra.push_back({2,6,7,8});
    ////this->tetrahedra.push_back({2,3,7,8});


    ////this->vertexDomainCoordinates.push_back({0.5,1,2});
    ////this->vertexCoordinatesF.push_back(7);
    ////this->vertexCoordinatesG.push_back(10);

    ////this->tetrahedra.push_back({6,7,8,9});

    ////this->vertexDomainCoordinates.push_back({2,1,2});
    ////this->vertexCoordinatesF.push_back(7);
    ////this->vertexCoordinatesG.push_back(10);

    ////this->tetrahedra.push_back({7,8,9,10});


    //// 9 2nd exit
    //this->vertexDomainCoordinates.push_back({1,0,2});
    //this->vertexCoordinatesF.push_back(-1);
    //this->vertexCoordinatesG.push_back(11);

    //this->tetrahedra.push_back({4,5,7,9});


    //// 10 2nd of 1st exit
    //this->vertexDomainCoordinates.push_back({-1,1,1});
    //this->vertexCoordinatesF.push_back(-2);
    //this->vertexCoordinatesG.push_back(9);

    //this->tetrahedra.push_back({8,2,6,10});
    //this->tetrahedra.push_back({0,2,6,10});
    //this->tetrahedra.push_back({0,4,6,10});


    //// 11 2nd on 2nd exit
    //this->vertexDomainCoordinates.push_back({0,-1,1});
    //this->vertexCoordinatesF.push_back(6);
    //this->vertexCoordinatesG.push_back(2);

    //this->tetrahedra.push_back({9,4,5,11});
    ////this->tetrahedra.push_back({0,1,5,11});

    //// 12 3rd on 2nd exit
    //this->vertexDomainCoordinates.push_back({0,0,2});
    //this->vertexCoordinatesF.push_back(-1);
    //this->vertexCoordinatesG.push_back(10);

    //this->tetrahedra.push_back({9,4,12,11});
    ////this->tetrahedra.push_back({0,1,5,11});


    //// 13 4th on 2nd exit
    //this->vertexDomainCoordinates.push_back({-1,0,2});
    //this->vertexCoordinatesF.push_back(-2);
    //this->vertexCoordinatesG.push_back(5);

    //this->tetrahedra.push_back({13,4,12,11});


    //this->tetrahedra.push_back({0,4,11,13});
    //this->tetrahedra.push_back({0,4,10,13});

    //this->tetrahedra.push_back({6,4,10,13});

    //this->tetrahedra.push_back({13,4,12,});


    //this->tetrahedra.push_back({4,5,7,10});

    //this->tetrahedra.push_back({5,6,7,9});
    //this->tetrahedra.push_back({4,5,7,9});

    //this->tetrahedra.push_back({4,7,5,9});
    //this->tetrahedra.push_back({4,7,6,9});


    //this->vertexDomainCoordinates.push_back({0,2,2});
    //this->vertexCoordinatesF.push_back(4);
    //this->vertexCoordinatesG.push_back(4);

    //this->tetrahedra.push_back({4,7,5,9});
    //this->tetrahedra.push_back({4,7,6,9});



    //this->vertexDomainCoordinates.push_back({2,2,1});
    //this->vertexCoordinatesF.push_back(0);
    //this->vertexCoordinatesG.push_back(0);

    //this->tetrahedra.push_back({5,7,9,10});


    //this->tetrahedra.push_back({2,7,3,10});


    this->computeMinMaxFG();


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

void
Data::computeMinMaxFG()
{
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

    // cout << "The min height is and the max height is " << minG << " " << maxG << endl;

    //
    // Add some padding
    //
    this->minF -= .2;
    this->maxF += .2;
    this->minG -= .2;
    this->maxG += .2;
}

std::vector<GLfloat> getBarycenter(const std::vector<GLfloat> &a, const std::vector<GLfloat> &b, const std::vector<GLfloat> &c)
{
    return {
        {.3 * a[0] + .3 * b[0] + .3 * c[0]},
        {.3 * a[1] + .3 * b[1] + .3 * c[1]},
        {.3 * a[2] + .3 * b[2] + .3 * c[2]}
    };
}

std::string pointToString(const std::vector<GLfloat> point)
{
    string pointString = "[";

    for (int i = 0 ; i < point.size() ; i++)
    {
        pointString += std::to_string(point[i]) + ", ";
    }

    pointString.pop_back();
    pointString.pop_back();
    pointString += "]";

    return pointString;
}

std::vector<GLfloat> getBarycenter(const std::vector<std::vector<GLfloat>> vertices)
{
    const float ratio = 1.0 / static_cast<float>(vertices.size());

    std::vector<float> barycenter = {0,0,0};

    for (int i = 0 ; i < vertices.size() ; i++)
    {
        barycenter[0] += ratio * vertices[i][0];
        barycenter[1] += ratio * vertices[i][1];
        barycenter[2] += ratio * vertices[i][2];
    }

    return barycenter;
}


void 
Data::generateSphereMesh()
{
    // Read Dimensions
    this->originalXdim = 10;
    this->originalYdim = 10;
    this->originalZdim = 10;
    this->originalTdim = 10;

    this->xdim = 2;
    this->ydim = 2;
    this->zdim = 2;
    this->tdim = 2;

    this->minX = 0;
    this->maxX = this->xdim;
    this->minY = 0;
    this->maxY = this->ydim;
    this->minZ = 0;
    this->maxZ = this->zdim;

    this->longnameF = "f";
    this->longnameG = "g";

    this->vertexDomainCoordinates.push_back({0, 0, 0});
    this->vertexDomainCoordinates.push_back({1, 0, 0});
    //this->vertexDomainCoordinates.push_back({0, 1, 0});
    this->vertexDomainCoordinates.push_back({.5, 0.8660, 0});
    //this->vertexDomainCoordinates.push_back({.25, .25, 1});

    // Forth vertex of the tetrahedron
    auto fourth = getBarycenter({this->vertexDomainCoordinates[0], this->vertexDomainCoordinates[1], this->vertexDomainCoordinates[2]});
    fourth[2] += 1;

    this->vertexDomainCoordinates.push_back(fourth);

    // Subdividing the sides
    this->vertexDomainCoordinates.push_back(getBarycenter({this->vertexDomainCoordinates[0], this->vertexDomainCoordinates[1], this->vertexDomainCoordinates[2]}));
    this->vertexDomainCoordinates.push_back(getBarycenter({this->vertexDomainCoordinates[0], this->vertexDomainCoordinates[1], this->vertexDomainCoordinates[3]}));
    this->vertexDomainCoordinates.push_back(getBarycenter({this->vertexDomainCoordinates[0], this->vertexDomainCoordinates[2], this->vertexDomainCoordinates[3]}));

    // Center
    this->vertexDomainCoordinates.push_back(getBarycenter({this->vertexDomainCoordinates[0], this->vertexDomainCoordinates[1], this->vertexDomainCoordinates[2], this->vertexDomainCoordinates[3]}));

    //this->vertexDomainCoordinates.push_back(getBarycenter({this->vertexDomainCoordinates[0], this->vertexDomainCoordinates[1], this->vertexDomainCoordinates[5]}));

    // Rebecca (3cc)
    //vector<string> tophat = {"v0", "V4", "V5", "v1", "V2", "v3"};

    vector<string> tophat = {
//"v0", "V4", "v2", "V3", "v1", "v5",
//"v0", "V4", "v3", "v2", "V5", "v1", "V6"
//"v0", "V4", "v3", "v2", "V6", "v1", "V5"
"v0", "v3", "v2", "v1", "V4", "V6", "V5"
};


    //vector<string> tophat = {"v0", "v4", "v5", "v1", "v2", "v3"};

    std::tie(this->vertexCoordinatesF, vertexCoordinatesG) = utility::generateCoordinatesFromTophat(tophat);

    this->computeMinMaxFG();

    //this->tetrahedra.push_back({4, 0, 1, 2});

    //this->tetrahedra.push_back({4, 0, 5, 7});
    //this->tetrahedra.push_back({4, 0, 1, 7});
    //this->tetrahedra.push_back({4, 1, 5, 7});

    // This is for one subdivision
    // {
    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 1, 4});
    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 2, 4});
    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 1, 2, 4});

    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 2, 3});
    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 1, 3});
    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 1, 2, 3});
    //// }


    // This is for two subdivisions
    // {
    this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 1, 4});
    this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 2, 4});
    this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 1, 2, 4});

    this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 1, 5});
    this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 3, 5});
    this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 1, 3, 5});

    this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 2, 6});
    this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 6, 3});
    this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 6, 2, 3});

    this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 1, 2, 3});
    // }

    this->tetsWithFibers = vector<bool>(this->tetrahedra.size(), false);

    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 3, 5});
    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 1, 3, 5});

    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 2, 6});
    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 0, 3, 6});
    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 2, 3, 6});

    //this->tetrahedra.push_back({this->vertexDomainCoordinates.size() - 1, 1, 2, 3});
}

void
Data::readNcData(tv9k::InputInformation input)
{
    //this->generateStandardGrid();
    this->generateSphereMesh();
}
