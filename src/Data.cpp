#include "Data.h"

#include <fstream>
#include <utility>

using namespace std;

void
Data::computeMinMaxFG()
{
    this->minX = this->vertexDomainCoordinates[0][0];
    this->maxX = this->vertexDomainCoordinates[0][0];

    this->minY = this->vertexDomainCoordinates[0][1];
    this->maxY = this->vertexDomainCoordinates[0][1];

    this->minZ = this->vertexDomainCoordinates[0][2];
    this->maxZ = this->vertexDomainCoordinates[0][2];

    for (int i = 0 ; i < this->vertexDomainCoordinates.size() ; i++)
    {
        this->minX = std::min(this->minX, this->vertexDomainCoordinates[i][0]);
        this->maxX = std::max(this->maxX, this->vertexDomainCoordinates[i][0]);

        this->minY = std::min(this->minY, this->vertexDomainCoordinates[i][1]);
        this->maxY = std::max(this->maxY, this->vertexDomainCoordinates[i][1]);

        this->minZ = std::min(this->minZ, this->vertexDomainCoordinates[i][2]);
        this->maxZ = std::max(this->maxZ, this->vertexDomainCoordinates[i][2]);
    }

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


void
Data::readNcData(tv9k::InputInformation input)
{
    //this->generateStandardGrid();
    //this->generateSphereMesh();

    // Read Dimensions
    this->longnameF = "f";
    this->longnameG = "g";

    std::ifstream dataFile (input.filename);

    if (false == dataFile.is_open()) { throw "Could not open data file."; }

    // Read in data in a string and skip the comments
    string rawStringData;
    string myline;
    while (dataFile) {
        std::getline (dataFile, myline);
        if (myline[0] == '#')
        {
            //std::cout << myline << '\n';
        }
        else
        {
            rawStringData += " " + myline;
        }
    }

    // Set up the inputstream
    std::istringstream dataStream(rawStringData);

    // Read in the number of vertices and tets
    int numVertices, numTets;
    dataStream >> numVertices >> numTets;

    cout << "This is the number of vertices and tets " << numVertices  << " " << numTets << endl;

    // Initialize all the arrays I need
    this->vertexCoordinatesF = std::vector<GLfloat>(numVertices, 0);
    this->vertexCoordinatesG = std::vector<GLfloat>(numVertices, 0);
    this->tetrahedra = std::vector<std::vector<size_t>>(numTets, {0, 0, 0, 0});
    this->vertexDomainCoordinates = std::vector<std::vector<GLfloat>>(numVertices, {0, 0, 0});

    // Read in the domain coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> this->vertexDomainCoordinates[i][0];
        dataStream >> this->vertexDomainCoordinates[i][1];
        dataStream >> this->vertexDomainCoordinates[i][2];
    }

    // Read in the range coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> this->vertexCoordinatesF[i];
        dataStream >> this->vertexCoordinatesG[i];
    }
    
    // Read in the tetrahedron configuration
    for  (int i = 0 ; i < numTets ; i++)
    {
        dataStream >> this->tetrahedra[i][0];
        dataStream >> this->tetrahedra[i][1];
        dataStream >> this->tetrahedra[i][2];
        dataStream >> this->tetrahedra[i][3];
    }

    this->computeMinMaxFG();
}
