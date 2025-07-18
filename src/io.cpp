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
#include "src/TetMesh.h"


TetMesh io::readData(const std::string &filename)
{
    std::filesystem::path filePath(filename);
    
    if (!std::filesystem::exists(filePath)) 
    {
        throw std::runtime_error("File does not exist: " + filename);
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

    throw std::runtime_error("Unsupported file type: " + extension);
}

TetMesh io::readDataVtu(const std::string &filename)
{
    TetMesh tetMesh;

    // Read the VTU file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());

    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();
    if (!mesh)
    {
        throw std::runtime_error("Failed to get mesh output from the file: " + filename);
    }

    if (mesh->GetNumberOfPoints() == 0)
    {
        throw std::runtime_error("Mesh contains no points: " + filename);
    }

    if (mesh->GetNumberOfCells() == 0)
    {
        throw std::runtime_error("Mesh contains no cells: " + filename);
    }

    // Set deault names for the range axis
    tetMesh.longnameF = "f";
    tetMesh.longnameG = "g";

    int numVertices = mesh->GetPoints()->GetNumberOfPoints(); 
    int numTets = mesh->GetNumberOfCells();

    // Initialize all the data arrays
    tetMesh.vertexCoordinatesF = std::vector<float>(numVertices, 0);
    tetMesh.vertexCoordinatesG = std::vector<float>(numVertices, 0);
    tetMesh.tetrahedra = std::vector<std::array<size_t, 4>>(numTets, {0, 0, 0, 0});
    tetMesh.vertexDomainCoordinates = std::vector<std::vector<float>>(numVertices, {0, 0, 0});

    // Print vertices
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
    //std::cout << "Vertices:\n";
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
        double p[3];
        points->GetPoint(i, p);
        //std::cout << "Vertex " << i << ": (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";

        tetMesh.vertexDomainCoordinates[i][0] = p[0];
        tetMesh.vertexDomainCoordinates[i][1] = p[1];
        tetMesh.vertexDomainCoordinates[i][2] = p[2];
    }

    // Print tetrahedra
    //std::cout << "\nTetrahedra:\n";
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++) {
        vtkCell* cell = mesh->GetCell(i);
        if (cell->GetNumberOfPoints() == 4) { // Tetrahedron check
            //std::cout << "Tetrahedron " << i << ": ";
            for (vtkIdType j = 0; j < 4; j++) {
                //std::cout << cell->GetPointId(j) << " ";
                tetMesh.tetrahedra[i][j] = cell->GetPointId(j);
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
        tetMesh.vertexCoordinatesF[i] = fDataArray->GetTuple1(i);
    }

    for (vtkIdType i = 0; i < gDataArray->GetNumberOfTuples(); i++) 
    {
        tetMesh.vertexCoordinatesG[i] = gDataArray->GetTuple1(i);
    }

    return tetMesh;
}


TetMesh io::readDataTxt(const std::string &filename)
{
    TetMesh tetMesh;
    
    // Set deault names for the range axis
    tetMesh.longnameF = "f";
    tetMesh.longnameG = "g";

    // Open data file
    std::ifstream dataFile (filename);
    if (false == dataFile.is_open()) 
    { 
        throw std::runtime_error("Could not open file: " + filename);
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
    tetMesh.vertexCoordinatesF = std::vector<float>(numVertices, 0);
    tetMesh.vertexCoordinatesG = std::vector<float>(numVertices, 0);
    tetMesh.tetrahedra = std::vector<std::array<size_t, 4>>(numTets, {0, 0, 0, 0});
    tetMesh.vertexDomainCoordinates = std::vector<std::vector<float>>(numVertices, {0, 0, 0});

    // Read in the domain coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> tetMesh.vertexDomainCoordinates[i][0];
        dataStream >> tetMesh.vertexDomainCoordinates[i][1];
        dataStream >> tetMesh.vertexDomainCoordinates[i][2];
    }

    // Read in the range coordinates
    for  (int i = 0 ; i < numVertices ; i++)
    {
        dataStream >> tetMesh.vertexCoordinatesF[i];
        dataStream >> tetMesh.vertexCoordinatesG[i];
    }
    
    // Read in the tetrahedron configuration
    for  (int i = 0 ; i < numTets ; i++)
    {
        dataStream >> tetMesh.tetrahedra[i][0];
        dataStream >> tetMesh.tetrahedra[i][1];
        dataStream >> tetMesh.tetrahedra[i][2];
        dataStream >> tetMesh.tetrahedra[i][3];
    }

    return tetMesh;
}



void io::saveSheets(const TetMesh &tetMesh, const Arrangement &arrangement, const ReebSpace &reebSpace, const std::string &outputSheetPolygonsFilename)
{
    // Save all the polygons
    std::filesystem::path filePathOutput(outputSheetPolygonsFilename);

    // Write to the file
    std::ofstream outFile(filePathOutput);
    if (!outFile) 
    {
        throw std::runtime_error("Error: Could not open file for writing: " + filePathOutput.string());
    }

    outFile << reebSpace.sheetPolygon.size() << std::endl;

    for (const auto &[sheetId, polygon] : reebSpace.sheetPolygon)
    {
        outFile << "SheetId = " << sheetId << std::endl;

        // Compute the controid so that we can pull all verties towards it
        CartesianPoint centroid = CGAL::centroid(polygon.vertices_begin(), polygon.vertices_end());

        outFile << "[";
        for (int i = 0 ; i < polygon.size() ; i++)
        {
            const CartesianPoint &point = polygon[i];
            double u = point.x();
            double v = point.y();

            // Interpolate closer to the centroid
            const double alpha = 0.5;
            u = (1 - alpha) * u + alpha * centroid.x();
            v = (1 - alpha) * v + alpha * centroid.y();

            outFile << u << ", " << v << ", " << 0;
            if (i < polygon.size() - 1)
            {
                outFile << ", ";

            }
        }
        outFile << "]" << std::endl;
    }

    outFile.close();
}
