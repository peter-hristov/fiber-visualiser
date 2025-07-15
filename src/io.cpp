#include <filesystem>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <random>
#include <ranges>

#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>

#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#include "./io.h"


Data* io::readData(const std::string &filename)
{
    std::filesystem::path filePath(filename);
    
    if (!std::filesystem::exists(filePath)) 
    {
        std::cerr << "File does not exist: " << filename;
        return nullptr;
    }

    std::string extension = filePath.extension().string();
    if (extension == ".vtu") 
    {
        return io::readDataVtu(filename);
    } 
    else if (extension == ".txt") 
    {
        return io::readDataTxt(filename);
    } 

    std::cerr << "Unsupported file type: " <<  extension;
    return nullptr;
}

Data* io::readDataVtu(const std::string &filename)
{
    Data *data = new Data();

    // Read the VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());

    try
    {
        reader->Update();
    }
    catch (const std::exception& e)
    {
        std::cerr << "VTK failed to read the file: " << e.what() << std::endl;
        return nullptr;
    }

    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();
    if (!mesh)
    {
        std::cerr << "Failed to get mesh output from the file: " << filename << std::endl;
        return nullptr;
    }

    if (mesh->GetNumberOfPoints() == 0)
    {
        std::cerr << "Mesh contains no points: " << filename << std::endl;
        return nullptr;
    }

    if (mesh->GetNumberOfCells() == 0)
    {
        std::cerr << "Mesh contains no cells: " << filename << std::endl;
        return nullptr;
    }

    // Set deault names for the range axis
    data->longnameF = "f";
    data->longnameG = "g";

    int numVertices = mesh->GetPoints()->GetNumberOfPoints(); 
    int numTets = mesh->GetNumberOfCells();

    // Initialize all the data arrays
    data->vertexCoordinatesF = std::vector<GLfloat>(numVertices, 0);
    data->vertexCoordinatesG = std::vector<GLfloat>(numVertices, 0);
    data->tetrahedra = std::vector<std::vector<size_t>>(numTets, {0, 0, 0, 0});
    data->vertexDomainCoordinates = std::vector<std::vector<GLfloat>>(numVertices, {0, 0, 0});

    // Print vertices
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
    //std::cout << "Vertices:\n";
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
        double p[3];
        points->GetPoint(i, p);
        //std::cout << "Vertex " << i << ": (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";

        data->vertexDomainCoordinates[i][0] = p[0];
        data->vertexDomainCoordinates[i][1] = p[1];
        data->vertexDomainCoordinates[i][2] = p[2];
    }

    // Print tetrahedra
    //std::cout << "\nTetrahedra:\n";
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++) {
        vtkCell* cell = mesh->GetCell(i);
        if (cell->GetNumberOfPoints() == 4) { // Tetrahedron check
            //std::cout << "Tetrahedron " << i << ": ";
            for (vtkIdType j = 0; j < 4; j++) {
                //std::cout << cell->GetPointId(j) << " ";
                data->tetrahedra[i][j] = cell->GetPointId(j);
            }
            //std::cout << "\n";
        }
    }

    // Print vertex data arrays
    //std::cout << "\nVertex Data Arrays:\n";
    vtkPointData* pointData = mesh->GetPointData();

    assert(pointData->GetNumberOfArrays() >= 2);

    vtkDataArray* fDataArray = pointData->GetArray(1);
    vtkDataArray* gDataArray = pointData->GetArray(0);

    assert(fDataArray->GetNumberOfTuples() == numVertices);
    assert(gDataArray->GetNumberOfTuples() == numVertices);

    for (vtkIdType i = 0; i < fDataArray->GetNumberOfTuples(); i++) 
    {
        data->vertexCoordinatesF[i] = fDataArray->GetTuple1(i);
    }

    for (vtkIdType i = 0; i < gDataArray->GetNumberOfTuples(); i++) 
    {
        data->vertexCoordinatesG[i] = gDataArray->GetTuple1(i);
    }

    return data;
}


Data* io::readDataTxt(const std::string &filename)
{
    Data *data = new Data();
    
    // Set deault names for the range axis
    data->longnameF = "f";
    data->longnameG = "g";

    // Open data file
    std::ifstream dataFile (filename);
    if (false == dataFile.is_open()) 
    { 
        std::cerr << "Could not open file: " << filename;
        return nullptr;
    }


    // Read in data in a string and skip the comments
    std::string rawStringData;
    std::string myline;
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
    data->vertexCoordinatesF = std::vector<GLfloat>(numVertices, 0);
    data->vertexCoordinatesG = std::vector<GLfloat>(numVertices, 0);
    data->tetrahedra = std::vector<std::vector<size_t>>(numTets, {0, 0, 0, 0});
    data->vertexDomainCoordinates = std::vector<std::vector<GLfloat>>(numVertices, {0, 0, 0});

    // Read in the domain coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> data->vertexDomainCoordinates[i][0];
        dataStream >> data->vertexDomainCoordinates[i][1];
        dataStream >> data->vertexDomainCoordinates[i][2];
    }

    // Read in the range coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> data->vertexCoordinatesF[i];
        dataStream >> data->vertexCoordinatesG[i];
    }
    
    // Read in the tetrahedron configuration
    for  (int i = 0 ; i < numTets ; i++)
    {
        dataStream >> data->tetrahedra[i][0];
        dataStream >> data->tetrahedra[i][1];
        dataStream >> data->tetrahedra[i][2];
        dataStream >> data->tetrahedra[i][3];
    }

    return data;
}
