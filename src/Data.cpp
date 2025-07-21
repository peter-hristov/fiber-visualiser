#include "Data.h"
#include "./Timer.h"
#include "./DisjointSet.h"
#include "src/CGALTypedefs.h"
#include "src/FiberPoint.h"

#include <CGAL/enum.h>
#include <filesystem>
#include <cassert>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <queue>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <ranges>

#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>

#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

using namespace std;


void Data::saveFibers()
{
    std::cout << "Saving fibers in " << this->fibersFile << std::endl;
    //std::cout << "The fiber has size " << this->faceFibers.size() << std::endl;  

    // 1. Create the points
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto idArray = vtkSmartPointer<vtkIntArray>::New();
    auto colourArray = vtkSmartPointer<vtkFloatArray>::New();

    idArray->SetName("SheetId");
    idArray->SetNumberOfComponents(1);

    colourArray->SetName("Colour");
    colourArray->SetNumberOfComponents(3);

    for (const FiberPoint &p : this->faceFibers)
    {
        // Insert points and corresponding IDs
        points->InsertNextPoint(p.point[0], p.point[1], p.point[2]);
        idArray->InsertNextValue(p.sheetId);

        const array<float, 3> sheetColour = this->fiberColours[this->reebSpace.sheetConsequitiveIndices[p.sheetId] % this->fiberColours.size()];
        float color[3] = {sheetColour[0], sheetColour[1], sheetColour[2]};
        colourArray->InsertNextTuple(color);
    }

    // 3. Create the cells (wrap polyline in cell array)
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    for (int i = 1 ; i < this->faceFibers.size() ; i+=2)
    {
        //cout << i << " " << this->faceFibers[i-1].sheetId << " " << this->faceFibers[i].sheetId << endl;
        //cout << i << " " << this->faceFibers[i-1].triangleId << " " << this->faceFibers[i].triangleId << endl;
        //printf("(%f, %f, %f) -> (%f, %f, %f)\n", this->faceFibers[i-1].point[0], this->faceFibers[i-1].point[1], this->faceFibers[i-1].point[2], this->faceFibers[i].point[0], this->faceFibers[i].point[1], this->faceFibers[i].point[2]);

        if (this->faceFibers[i-1].sheetId == this->faceFibers[i].sheetId)
        {
            // One edge segment
            auto polyLine = vtkSmartPointer<vtkPolyLine>::New();
            polyLine->GetPointIds()->SetNumberOfIds(2);
            polyLine->GetPointIds()->SetId(0, i-1);
            polyLine->GetPointIds()->SetId(1, i);

            cells->InsertNextCell(polyLine);
        }
    }

    // 4. Create the polydata object
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetLines(cells);

    // 5. Attach the VertexID array to the point data
    polyData->GetPointData()->AddArray(idArray);
    polyData->GetPointData()->AddArray(colourArray);
    polyData->GetPointData()->SetScalars(colourArray);  // optional: for coloring

    // 6. Write to .vtp file (XML format)
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(this->fibersFile.c_str());
    writer->SetInputData(polyData);
    writer->Write();
}

void Data::generatefFaceFibersForSheets(const int sheetOutputCount, const int numberOfFiberPoints, const std::string folderPath)
{
    namespace fs = std::filesystem;

    //string folderPath = "./sheetFibers";
    fs::path folderPathFs(folderPath);
    if (!fs::exists(folderPathFs)) 
    {
        fs::create_directory(folderPathFs);
    }

    for (const auto &[sheetId, colourId] : this->reebSpace.sheetConsequitiveIndices)
    {
        if (this->reebSpace.incompleteSheets.contains(sheetId))
        {
            printf("Skipping fiber %d, it's incomplete.",  sheetId);
        }

        if (colourId > sheetOutputCount || this->reebSpace.incompleteSheets.contains(sheetId))
        {
            continue;
        }

        std::cout << "-------------------------------------------------------------------------------------------- Generating fibers for sheet " << sheetId << "..." << std::endl;
        this->generatefFaceFibersForSheet(sheetId, numberOfFiberPoints);

        //std::cout << "Saving fibers..." << std::endl;
        this->fibersFile = folderPathFs.string() + "/fibers_" + std::to_string(sheetId) + ".vtp";
        this->saveFibers();
        this->faceFibers.clear();
    }

    this->fibersFile = "./fibers.vtp";

}

void Data::generatefFaceFibersForSheet(const int sheetId, const int numberOfFiberPoints)
{
    CartesianPolygon_2 &polygon = this->reebSpace.sheetPolygon[sheetId];

    if (polygon.size() == 0)
    {
        return;
    }


    // Compute the controid so that we can pull all verties towards it
    CartesianPoint centroid = CGAL::centroid(polygon.vertices_begin(), polygon.vertices_end());

    // To make sure we don't write a comma at the end of the array

    vector<vector<float>> fiberPoints;

    // If need only one, get it at the center
    if (numberOfFiberPoints == 1)
    {
        fiberPoints.push_back({(float)centroid.x(), (float)centroid.y()});
        this->computeTetExitPointsNewNew(fiberPoints[0][0], fiberPoints[0][1], false, sheetId);

        return;
    }


    // If we need more, sample along the boundary
    for (const CartesianPoint &point : polygon) 
    {
        // Get point from CGAL (and convert to double )
        float u = point.x();
        float v = point.y();

        // Interpolate closer to the centroid to make sure we are in the sheet ( if the sheet is "convex enough")
        float alpha = 0.2;
        u = (1 - alpha) * u + alpha * centroid.x();
        v = (1 - alpha) * v + alpha * centroid.y();

        fiberPoints.push_back({u, v});
        //fiberPoints.push_back({centroid.x(), centroid.y()});
    }


    // Calculate step size we only want some of the fiber points, not all
    double step = static_cast<double>(fiberPoints.size() - 1) / (numberOfFiberPoints - 1);

    for (int i = 0; i < numberOfFiberPoints; ++i) 
    {
        int index = static_cast<int>(i * step);
        this->computeTetExitPointsNewNew(fiberPoints[index][0], fiberPoints[index][1], false, sheetId);
    }
}


void Data::printSheetHistogram()
{
    Timer::start();
    if (this->currentFiberPoint.size() == 0)
    {
        return;
    }

    assert (this->currentFiberPoint.size() == 2);

    set<int> intersectedSheets;
    CartesianLine line(CartesianPoint(0, this->currentFiberPoint[0]), CartesianPoint(1, this->currentFiberPoint[1])); // Line through (0, 1) and (1, 0)
    for (const auto &[sheetId, polygon] : this->reebSpace.sheetPolygon)
    {
        QVector<QPoint> points;

        for (const CartesianSegment& segment : polygon.edges()) 
        {
            if (CGAL::do_intersect(segment, line)) 
            {
                intersectedSheets.insert(sheetId);
            }
        }
    }

    std::cout << "Intersected sheets: ";
    for (const auto &sheetId : intersectedSheets)
    {
        std::cout << sheetId << " (a = )" << this->reebSpace.sheetArea[sheetId] << std::endl;

    }
    std::cout << "\n";
    Timer::stop("Computed face intersected by the line  :");
}


void Data::computeTetExitPointsNewNew(const GLfloat u, const GLfloat v, const bool clearFibers, const int reebSheetIdOnly)
{
    if (true == clearFibers)
    {
        this->faceFibers.clear();
    }


    //
    // Get the ID of the face we are intersecting
    //

    //Timer::start();

    // Store the current fiber point
    this->currentFiberPoint = {u, v};

    // The query point (u, v)
    Point_2 query_point(u, v);


    // Locate the point in the arrangement
    CGAL::Object result = this->arrangement.pl->locate(query_point);

    // Try to assign to a face, edge or a vertex
    Arrangement_2::Face_const_handle face;
    Arrangement_2::Halfedge_const_handle edge;
    Arrangement_2::Vertex_const_handle vertex;

    int currentFaceID = 0;

    if (CGAL::assign(face, result)) 
    {
        currentFaceID = this->arrangement.arrangementFacesIdices[face];
        printf("Found query point.");
    } 
    // If we are on an edge, just grad an adjacent face
    else if (CGAL::assign(edge, result)) 
    {
        face = edge->face();
        currentFaceID = this->arrangement.arrangementFacesIdices[face];
        printf("Found query point.");
    } 
    // If we are on a vertex grab an indicent edge and get its face
    else if (CGAL::assign(vertex, result)) 
    {
        edge = vertex->incident_halfedges();
        face = edge->face();
        currentFaceID = this->arrangement.arrangementFacesIdices[face];
        printf("Found query point.");
    } else 
    {
        printf("NOT Found query point.");
        assert(false);
    }
    //Timer::stop("Computed active arrangement face       :");

    printf("The current faces ID is %d\n", currentFaceID);


    //Timer::start();

    // The sizes of these data structures are linear in the size of the fiber, not an issue
    std::queue<int> bfsQueue;
    // This also acts as the visited array
    std::unordered_map<int, int> triangleSheetId;
    // Needed to close loops closed fibers, hack because we are not visiting tets, but triangles
    std::unordered_set<pair<int, int>, MyHash<pair<int, int>>> activeAdjacentTrianglesConnected;
    // Cache barycentric coordintes, they are expensive to compute
    std::unordered_map<int, std::array<double, 3>> triangleBarycentricCoordinates;

    if (reebSheetIdOnly == -1)
    {
        std::cout << "There are " << this->reebSpace.fiberSeeds[currentFaceID].size() << " fiber components with (sheet IDs, sorted IDs): ";
    }

    //vector<int> sheetIds;
    for (const auto &[fiberComponentId, triangleId] : this->reebSpace.fiberSeeds[currentFaceID])
    {
        const int sheetId = this->reebSpace.correspondenceGraph.findElement({currentFaceID, fiberComponentId});

        if (reebSheetIdOnly == -1 || sheetId == reebSheetIdOnly)
        {
            bfsQueue.push(triangleId);
            cout << "Adding seed \n";
        }

        //const int sheetId = this->reebSpace.findTriangle({currentFaceID, fiberComponentId});
        //triangleColour[triangleId] = this->sheetConsequitiveIndices[sheetId] % this->fiberColours.size();
        triangleSheetId[triangleId] = sheetId;

        //sheetIds.push_back(sheetId);

        if (reebSheetIdOnly == -1)
        {
            printf("(%d, %d) ", sheetId, this->reebSpace.sheetConsequitiveIndices[sheetId]);
        }
    }
    //std::cout << std::endl;


    // Define query point
    CartesianPoint P(u, v);

    while (false == bfsQueue.empty())
    {
        const int currentTriangleId = bfsQueue.front();
        const int currentSheeId = triangleSheetId[currentTriangleId];
        bfsQueue.pop();

        const std::array<float, 3> sheetColour = this->fiberColours[this->reebSpace.sheetConsequitiveIndices[currentSheeId] % this->fiberColours.size()];

        const set<int> triangleUnpacked = this->tetMesh.triangles[currentTriangleId];
        const vector<int> triangleIndices = std::vector<int>(triangleUnpacked.begin(), triangleUnpacked.end());

        std::array<double, 3> barycentricCoordinatesCurrent;

        if (triangleBarycentricCoordinates.contains(currentTriangleId))
        {
            barycentricCoordinatesCurrent = triangleBarycentricCoordinates[currentTriangleId];
        }
        else
        {
            const CartesianPoint A(this->tetMesh.vertexCoordinatesF[triangleIndices[0]], this->tetMesh.vertexCoordinatesG[triangleIndices[0]]);
            const CartesianPoint B(this->tetMesh.vertexCoordinatesF[triangleIndices[1]], this->tetMesh.vertexCoordinatesG[triangleIndices[1]]);
            const CartesianPoint C(this->tetMesh.vertexCoordinatesF[triangleIndices[2]], this->tetMesh.vertexCoordinatesG[triangleIndices[2]]);
            CGAL::Barycentric_coordinates::triangle_coordinates_2(A, B, C, P, barycentricCoordinatesCurrent.begin());
            triangleBarycentricCoordinates[currentTriangleId] = barycentricCoordinatesCurrent;
        }

        // Sanity check
        assert(barycentricCoordinatesCurrent[0] > 0 && barycentricCoordinatesCurrent[1] > 0 && barycentricCoordinatesCurrent[2] > 0);

        // Look at the neighbours
        for (const int &neighbourTriagleId : this->tetMesh.tetIncidentTriangles[currentTriangleId])
        {
            if (neighbourTriagleId == currentTriangleId)
            {
                continue;
            }

            // We can't skip visited neighbours, because there may be a fiber between us (completing a circle)
            //if (triangleColour.contains(neighbourTriagle)) { continue; }

            const set<int> triangle2Unpacked = this->tetMesh.triangles[neighbourTriagleId];
            const vector<int> triangle2Indices = std::vector<int>(triangle2Unpacked.begin(), triangle2Unpacked.end());



            // The neighbour is active if we've already seen it
            bool isActive = triangleSheetId.contains(neighbourTriagleId);

            // Or if the image of the triangle contains the fiber points
            // We use a fast test to avoid having to use barycentric coordinates all the time
            if (false == isActive)
            {
                CartesianPoint A(this->tetMesh.vertexCoordinatesF[triangle2Indices[0]], this->tetMesh.vertexCoordinatesG[triangle2Indices[0]]);
                CartesianPoint B(this->tetMesh.vertexCoordinatesF[triangle2Indices[1]], this->tetMesh.vertexCoordinatesG[triangle2Indices[1]]);
                CartesianPoint C(this->tetMesh.vertexCoordinatesF[triangle2Indices[2]], this->tetMesh.vertexCoordinatesG[triangle2Indices[2]]);

                std::vector<CartesianPoint> triangle = {A, B, C};
                const auto result = CGAL::bounded_side_2(triangle.begin(), triangle.end(), P);

                isActive = (result == CGAL::ON_BOUNDED_SIDE);
            }

            // Determine if the triangle is active
            if (isActive)
            {
                // Only add the neighbour if we have not already visited it
                if (false == triangleSheetId.contains(neighbourTriagleId))
                {
                    // BFS things
                    bfsQueue.push(neighbourTriagleId);
                    triangleSheetId[neighbourTriagleId] = currentSheeId;
                }

                // Even if we have aleady added a neighbour, maybe there still isn't a fiber between us (for finishing loops)

                // At this point, we know that both us and we neighbour are active, is there already a fiber between us? Then skip
                if (activeAdjacentTrianglesConnected.contains({currentTriangleId, neighbourTriagleId}))
                {
                    continue;
                }
                else
                {
                    activeAdjacentTrianglesConnected.insert({currentTriangleId, neighbourTriagleId});
                    activeAdjacentTrianglesConnected.insert({neighbourTriagleId, currentTriangleId});
                }




                // Compute barycentric coordinates for drawing
                std::array<double, 3> barycentricCoordinatesNeighbour;
                if (triangleBarycentricCoordinates.contains(neighbourTriagleId))
                {
                    barycentricCoordinatesNeighbour = triangleBarycentricCoordinates[neighbourTriagleId];
                }
                else
                {
                    CartesianPoint A(this->tetMesh.vertexCoordinatesF[triangle2Indices[0]], this->tetMesh.vertexCoordinatesG[triangle2Indices[0]]);
                    CartesianPoint B(this->tetMesh.vertexCoordinatesF[triangle2Indices[1]], this->tetMesh.vertexCoordinatesG[triangle2Indices[1]]);
                    CartesianPoint C(this->tetMesh.vertexCoordinatesF[triangle2Indices[2]], this->tetMesh.vertexCoordinatesG[triangle2Indices[2]]);

                    CGAL::Barycentric_coordinates::triangle_coordinates_2(A, B, C, P, barycentricCoordinatesNeighbour.begin());
                    triangleBarycentricCoordinates[neighbourTriagleId] = barycentricCoordinatesNeighbour;
                }

                


                //
                // Add a fiber segment
                //
                FiberPoint fb(
                        barycentricCoordinatesCurrent[0], 
                        barycentricCoordinatesCurrent[1], 
                        {
                            this->tetMesh.vertexDomainCoordinates[triangleIndices[0]],
                            this->tetMesh.vertexDomainCoordinates[triangleIndices[1]],
                            this->tetMesh.vertexDomainCoordinates[triangleIndices[2]],
                        },
                        sheetColour);
                fb.sheetId = currentSheeId;
                fb.triangleId = currentTriangleId;
                this->faceFibers.push_back(fb);

                FiberPoint fb2(barycentricCoordinatesNeighbour[0], barycentricCoordinatesNeighbour[1], {
                        this->tetMesh.vertexDomainCoordinates[triangle2Indices[0]],
                        this->tetMesh.vertexDomainCoordinates[triangle2Indices[1]],
                        this->tetMesh.vertexDomainCoordinates[triangle2Indices[2]],
                        },
                        sheetColour);
                fb2.sheetId = currentSheeId;
                fb2.triangleId = neighbourTriagleId;
                this->faceFibers.push_back(fb2);

                //printf("Adding fiber between %d -> %d\n", currentTriangleId, neighbourTriagleId);
            }
        }
    }

    //Timer::stop("Computed fiber in                      :");
}

//void Data::computeTetExitPointsNew(const GLfloat u, const GLfloat v, const std::vector<float> color)
//{
    ////this->faceFibers.clear();

    ////
    //// Get the ID of the face we are intersecting
    ////

    ////Timer::start();

    //// The query point (u, v)
    //Point_2 query_point(u, v);


    //// Locate the point in the arrangement
    //CGAL::Object result = this->arrangement.pl->locate(query_point);

    //// Try to assign to a face, edge or a vertex
    //Arrangement_2::Face_const_handle face;
    //Arrangement_2::Halfedge_const_handle edge;
    //Arrangement_2::Vertex_const_handle vertex;

    //int currentFaceID = 0;

    //if (CGAL::assign(face, result)) 
    //{
        //currentFaceID = this->arrangement.arrangementFacesIdices[face];
    //} 
    //// If we are on an edge, just grad an adjacent face
    //else if (CGAL::assign(edge, result)) 
    //{
        //face = edge->face();
        //currentFaceID = this->arrangement.arrangementFacesIdices[face];
    //} 
    //// If we are on a vertex grab an indicent edge and get its face
    //else if (CGAL::assign(vertex, result)) 
    //{
        //edge = vertex->incident_halfedges();
        //face = edge->face();
        //currentFaceID = this->arrangement.arrangementFacesIdices[face];
    //} else 
    //{
        //assert(false);
    //}
    ////Timer::stop("Computed active arrangement face       :");



    //int i = -1;
    //for (const auto &[triangle, triangleId] : this->reebSpace.preimageGraphs[currentFaceID].data)
    //{
        //i++;
        //int j = -1;
        //for (const auto&[triangle2, triangleId2] : this->reebSpace.preimageGraphs[currentFaceID].data)
        //{
            //j++;
            //if (j <= i) { continue; }

            //const set<int> triangleUnpacked = this->tetMesh.triangles[triangle];
            //const set<int> triangle2Unpacked = this->tetMesh.triangles[triangle2];

            //if (this->tetMesh.connectedTriangles.contains({triangleUnpacked, triangle2Unpacked}))
            //{
                //const int componentID = this->reebSpace.preimageGraphs[currentFaceID].find(triangleId);
                ////const int pairToHIndex = this->vertexHtoIndex[{currentFaceID, componentID}];
                //const int sheetID = this->reebSpace.reebSpace.findTriangle({currentFaceID, componentID});
                //const int sheetColourID = this->reebSpace.sheetConsequitiveIndices[sheetID] % this->fiberColours.size();
                //const vector<float> sheetColour = this->fiberColours[sheetColourID];

                ////
                //// Get the IDs and barycentri coordinates for the first point
                ////
                //vector<int> vertexIds;

                //for(const int &vertexId : triangleUnpacked)
                //{
                    //vertexIds.push_back(vertexId);
                //}

                //CartesianPoint A(this->tetMesh.vertexCoordinatesF[vertexIds[0]], this->tetMesh.vertexCoordinatesG[vertexIds[0]]);
                //CartesianPoint B(this->tetMesh.vertexCoordinatesF[vertexIds[1]], this->tetMesh.vertexCoordinatesG[vertexIds[1]]);
                //CartesianPoint C(this->tetMesh.vertexCoordinatesF[vertexIds[2]], this->tetMesh.vertexCoordinatesG[vertexIds[2]]);

                //// Define query point
                //CartesianPoint P(u, v);

                //std::array<double, 3> coordinates;
                //CGAL::Barycentric_coordinates::triangle_coordinates_2(A, B, C, P, coordinates.begin());

                //assert(coordinates[0] >= 0 && coordinates[1] >= 0 && coordinates[2] >= 0);

                //FaceFiberPoint fb(coordinates[0], coordinates[1], {
                        //this->tetMesh.vertexDomainCoordinates[vertexIds[0]],
                        //this->tetMesh.vertexDomainCoordinates[vertexIds[1]],
                        //this->tetMesh.vertexDomainCoordinates[vertexIds[2]],
                        //},
                        //sheetColour);
                //this->faceFibers.push_back(fb);




                //// Get the IDs and barycentri coordinates for the second point
                //vector<int> vertexIds2;

                //for(const int &vertexId : triangle2Unpacked)
                //{
                    //vertexIds2.push_back(vertexId);
                //}


                //// Define triangle vertices
                //CartesianPoint A2(this->tetMesh.vertexCoordinatesF[vertexIds2[0]], this->tetMesh.vertexCoordinatesG[vertexIds2[0]]);
                //CartesianPoint B2(this->tetMesh.vertexCoordinatesF[vertexIds2[1]], this->tetMesh.vertexCoordinatesG[vertexIds2[1]]);
                //CartesianPoint C2(this->tetMesh.vertexCoordinatesF[vertexIds2[2]], this->tetMesh.vertexCoordinatesG[vertexIds2[2]]);

                //std::array<double, 3> coordinates2;
                //CGAL::Barycentric_coordinates::triangle_coordinates_2(A2, B2, C2, P, coordinates2.begin());
                //assert(coordinates2[0] >= 0 && coordinates2[1] >= 0 && coordinates2[2] >= 0);

                //FaceFiberPoint fb2(coordinates2[0], coordinates2[1], {
                        //this->tetMesh.vertexDomainCoordinates[vertexIds2[0]],
                        //this->tetMesh.vertexDomainCoordinates[vertexIds2[1]],
                        //this->tetMesh.vertexDomainCoordinates[vertexIds2[2]],
                        //},
                        //sheetColour);

                //this->faceFibers.push_back(fb2);
            //}
        //}
    //}
//}







void Data::computeTetExitPoints(const GLfloat u, const GLfloat v, const std::vector<float> color)
{
    this->faceFibers.clear();
    this->tetsWithFibers = vector<bool>(this->tetMesh.tetrahedra.size(), false);

    //
    // Get the ID of the face we are intersecting
    //

    // The query point (u, v)
    Point_2 query_point(u, v);

    // Locate the point in the arrangement
    CGAL::Object result = this->arrangement.pl->locate(query_point);

    // Try to assign to a face, edge or a vertex
    Arrangement_2::Face_const_handle face;
    Arrangement_2::Halfedge_const_handle edge;
    Arrangement_2::Vertex_const_handle vertex;

    int currentFaceID = 0;

    if (CGAL::assign(face, result)) 
    {
        currentFaceID = this->arrangement.arrangementFacesIdices[face];
    } 
    // If we are on an edge, just grad an adjacent face
    else if (CGAL::assign(edge, result)) 
    {
        face = edge->face();
        currentFaceID = this->arrangement.arrangementFacesIdices[face];
    } 
    // If we are on a vertex grab an indicent edge and get its face
    else if (CGAL::assign(vertex, result)) 
    {
        edge = vertex->incident_halfedges();
        face = edge->face();
        currentFaceID = this->arrangement.arrangementFacesIdices[face];
    } else 
    {
        assert(false);
    }

    // For every tet, compute the two exit points
    for(size_t tetId = 0 ; tetId < this->tetMesh.tetrahedra.size(); tetId++)
    {
        const auto tet = this->tetMesh.tetrahedra[tetId];

        // For every triangle in every tet, get the fiber point in it
        for(int i = 0 ; i < 4 ; i++)
        {
            for(int j = i + 1 ; j < 4 ; j++)
            {
                for(int k = j + 1 ; k < 4 ; k++)
                {
                    float x1 = this->tetMesh.vertexCoordinatesF[tet[i]];
                    float y1 = this->tetMesh.vertexCoordinatesG[tet[i]];

                    float x2 = this->tetMesh.vertexCoordinatesF[tet[j]];
                    float y2 = this->tetMesh.vertexCoordinatesG[tet[j]];

                    float x3 = this->tetMesh.vertexCoordinatesF[tet[k]];
                    float y3 = this->tetMesh.vertexCoordinatesG[tet[k]];


                    const float xmin = std::min({x1, x2, x3});
                    const float xmax = std::max({x1, x2, x3});
                    const float ymin = std::min({y1, y2, y3});
                    const float ymax = std::max({y1, y2, y3});

                    // This triangle is not relevant, point is outside the bounding box
                    if (u < xmin || u > xmax || v < ymin || v > ymax) 
                    {
                        continue;
                    }


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
                        const int triangleID = this->tetMesh.triangleIndices[std::set<int>({triangleVertexA, triangleVertexB, triangleVertexC})];
                        const int componentID = this->reebSpace.preimageGraphs[currentFaceID].findElement(triangleID);
                        //printf("The face ID is %d and the component ID is = %d\n", currentFaceID, componentID);

                        // 2. Fac ComponentID -> Reeb Space Sheet
                        //const int pairToHIndex = this->vertexHtoIndex[{currentFaceID, componentID}];
                        //const int sheetID = this->reebSpace.findTriangle(pairToHIndex);
                        const int sheetID = this->reebSpace.correspondenceGraph.findElement({currentFaceID, componentID});

                        //printf("The Sheet ID is = %d\n", sheetID);

                        const int sheetColourID = this->reebSpace.sheetConsequitiveIndices[sheetID];

                        // 3. Get the colou of the sheet
                        const array<float, 3> sheetColour = this->fiberColours[sheetColourID];


                        FiberPoint fb(alpha, beta, {
                                this->tetMesh.vertexDomainCoordinates[tet[i]],
                                this->tetMesh.vertexDomainCoordinates[tet[j]],
                                this->tetMesh.vertexDomainCoordinates[tet[k]],
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



