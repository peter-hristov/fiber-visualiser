#include "Data.h"

#include <fstream>
#include <sstream>
#include <utility>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>


using namespace std;

void Data::computeTetExitPoints(const GLfloat u, const GLfloat v, const std::vector<float> color)
{
    this->tetsWithFibers = vector<bool>(this->tetrahedra.size(), false);


    // Find the point in the arrangement

    // Create the point-location query structure
    Point_location pl(arr);

    // The query point (u, v)
    Point_2 query_point(u, v);

    // Locate the point in the arrangement
    CGAL::Object result = pl.locate(query_point);

    // Try to assign to a face, edge or a vertex
    Arrangement_2::Face_const_handle face;
    Arrangement_2::Halfedge_const_handle edge;
    Arrangement_2::Vertex_const_handle vertex;

    int currentFaceID = 0;


    if (CGAL::assign(face, result)) 
    {
        currentFaceID = this->arrangementFacesIdices[face];
    } 
    // If we are on an edge, just grad an adjacent face
    else if (CGAL::assign(edge, result)) 
    {
        face = edge->face();
        currentFaceID = this->arrangementFacesIdices[face];
    } 
    // If we are on a vertex grab an indicent edge and get its face
    else if (CGAL::assign(vertex, result)) 
    {
        edge = vertex->incident_halfedges();
        face = edge->face();
        currentFaceID = this->arrangementFacesIdices[face];
    } else {
        assert(false);
    }



    //cout << endl << endl;

    // For every tet, compute the two exit points
    for(size_t tetId = 0 ; tetId < this->tetrahedra.size(); tetId++)
    {
        const auto tet = this->tetrahedra[tetId];

        // For every triangle in every tet, get the fiber in it
        for(int i = 0 ; i < 4 ; i++)
        {
            for(int j = i + 1 ; j < 4 ; j++)
            {
                for(int k = j + 1 ; k < 4 ; k++)
                {
                    float x1 = this->vertexCoordinatesF[tet[i]];
                    float y1 = this->vertexCoordinatesG[tet[i]];

                    float x2 = this->vertexCoordinatesF[tet[j]];
                    float y2 = this->vertexCoordinatesG[tet[j]];

                    float x3 = this->vertexCoordinatesF[tet[k]];
                    float y3 = this->vertexCoordinatesG[tet[k]];

                    float det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

                    float alpha = ((y2 - y3) * (u - x3) + (x3 - x2) * (v - y3)) / det;
                    float beta = ((y3 - y1) * (u - x3) + (x1 - x3) * (v - y3)) / det;
                    float gamma = 1 - alpha - beta;

                    // Are we inside the triangle. We exclude the 0 and 1 because weird things happen there
                    if (alpha > 0 && beta > 0 && gamma > 0 && alpha < 1 && beta < 1 && gamma < 1)
                    {
                        //printf("In triangle %ld, %ld, %ld in tet %ld.\n", tet[i], tet[j], tet[k], tetId);
                        //printf("In triangle (%f, %f) | (%f, %f) | (%f, %f) comparing with point (%f, %f) and alpha = %f, betta = %f, gamma = %f.\n", x1, y1, x2, y2, x3, y3, x, y, alpha, betta, gamma);

                        //const std::set<int> triangle({tet[i], tet[j], tet[k]});

                        const int triangleVertexA = tet[i];
                        const int triangleVertexB = tet[j];
                        const int triangleVertexC = tet[k];

                        //printf("Current triangle (%d, %d, %d)\n", triangleVertexA, triangleVertexB, triangleVertexC);


                        //for (const auto [key, value] : this->faceDisjointSets[currentFaceID].data)
                        //{
                            //cout << "Triangle ";
                            //for (const auto v : key)
                            //{
                                //cout << v << " ";
                            //}

                            //cout << " with root " << this->faceDisjointSets[currentFaceID].find(value) << " ( " << this->faceDisjointSets[currentFaceID].findTriangle(key) << ") " << endl;

                        //}


                        //for (const auto &[key, value] : this->reebSpace.data)
                        //{
                            //printf("Face ID = %d, fiber component root = %d, SheetID = %d\n", key.first, key.second, this->reebSpace.findTriangle(key));

                        //}


                        // Which sheets does this fiber belong to?
                        // 1. Triangle -> Face ComponentID
                        const int componentID = this->faceDisjointSets[currentFaceID].findTriangle(std::set<int>({triangleVertexA, triangleVertexB, triangleVertexC}));
                        //printf("The face ID is %d and the component ID is = %d\n", currentFaceID, componentID);

                        // 2. Fac ComponentID -> Reeb Space Sheet
                        const int sheetID = this->reebSpace.findTriangle({currentFaceID, componentID});

                        //printf("The Sheet ID is = %d\n", sheetID);

                        const int sheetColourID = this->sheetToColour[sheetID];

                        // 3. Get the colou of the sheet
                        const vector<float> sheetColour = this->fiberColours[sheetColourID];


                        FaceFiberPoint fb(alpha, beta, {
                                this->vertexDomainCoordinates[tet[i]],
                                this->vertexDomainCoordinates[tet[j]],
                                this->vertexDomainCoordinates[tet[k]],
                                },
                                sheetColour);

                        this->faceFibers.push_back(fb);
                        this->tetsWithFibers[tetId] = true;
                    }
                }
            }
        }
    }
}

void
Data::computeMinMaxRangeDomainCoordinates()
{
    // Compute the min/max domain coordinates
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


    // Compute the min/max range coordinates
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

    // Add some padding to the range coordinates for better visibility
    //this->minF -= .2;
    //this->maxF += .2;
    //this->minG -= .2;
    //this->maxG += .2;

    this->minF -= 0.2 * (this->maxF - this->minF);
    this->maxF += 0.2 * (this->maxF - this->minF);
    this->minG -= 0.2 * (this->maxG - this->minG);
    this->maxG += 0.2 * (this->maxG - this->minG);
}

void Data::readDataVTK(string filename)
{
    // Read the VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();




    // Set deault names for the range axis
    this->longnameF = "f";
    this->longnameG = "g";

    int numVertices = mesh->GetPoints()->GetNumberOfPoints(); 
    int numTets = mesh->GetNumberOfCells();

    // Initialize all the data arrays
    this->vertexCoordinatesF = std::vector<GLfloat>(numVertices, 0);
    this->vertexCoordinatesG = std::vector<GLfloat>(numVertices, 0);
    this->tetrahedra = std::vector<std::vector<size_t>>(numTets, {0, 0, 0, 0});
    this->vertexDomainCoordinates = std::vector<std::vector<GLfloat>>(numVertices, {0, 0, 0});



    // Print vertices
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
    std::cout << "Vertices:\n";
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
        double p[3];
        points->GetPoint(i, p);
        std::cout << "Vertex " << i << ": (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";

        this->vertexDomainCoordinates[i][0] = p[0];
        this->vertexDomainCoordinates[i][1] = p[1];
        this->vertexDomainCoordinates[i][2] = p[2];
    }

    // Print tetrahedra
    std::cout << "\nTetrahedra:\n";
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++) {
        vtkCell* cell = mesh->GetCell(i);
        if (cell->GetNumberOfPoints() == 4) { // Tetrahedron check
            std::cout << "Tetrahedron " << i << ": ";
            for (vtkIdType j = 0; j < 4; j++) {
                std::cout << cell->GetPointId(j) << " ";
                this->tetrahedra[i][j] = cell->GetPointId(j);
            }
            std::cout << "\n";
        }
    }

    // Print vertex data arrays
    std::cout << "\nVertex Data Arrays:\n";
    vtkPointData* pointData = mesh->GetPointData();

    assert(pointData->GetNumberOfArrays() >= 2);

    vtkDataArray* fDataArray = pointData->GetArray(0);
    vtkDataArray* gDataArray = pointData->GetArray(1);

    assert(fDataArray->GetNumberOfTuples() == numVertices);
    assert(gDataArray->GetNumberOfTuples() == numVertices);

    for (vtkIdType i = 0; i < fDataArray->GetNumberOfTuples(); i++) 
    {

        this->vertexCoordinatesF[i] = fDataArray->GetTuple1(i);
    }

    for (vtkIdType i = 0; i < gDataArray->GetNumberOfTuples(); i++) 
    {
        this->vertexCoordinatesG[i] = gDataArray->GetTuple1(i);
    }

}

void
Data::readData(string filename)
{
    // Set deault names for the range axis
    this->longnameF = "f";
    this->longnameG = "g";

    // Open data file
    std::ifstream dataFile (filename);
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

    // Set up the inputstream from the string
    std::istringstream dataStream(rawStringData);

    // Read in the number of vertices and tets
    int numVertices, numTets;
    dataStream >> numVertices >> numTets;

    // Initialize all the data arrays
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
}
