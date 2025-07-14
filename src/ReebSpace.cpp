#include "./ReebSpace.h"
#include "./DisjointSet.h"
#include "./CGALTypedefs.h"
#include "./Timer.h"
#include "./utility/indicators.hpp"

#include <cassert>
#include <cstddef>
#include <map>
#include <queue>
#include <set>
#include <iterator>
#include <unordered_map>

std::pair<std::set<int>, std::set<int>> ReebSpace::getMinusPlusTrianglesIndex(Arrangement_2::Halfedge_const_handle currentHalfEdge, Data *data)
{
    // Step 1. Initialize lists
    std::set<int> plusTriangles;
    std::set<int> minusTriangles;

    // Step 2. Find the edge in the mesh corresponding to the segment corresponding to the half edge
    const Segment_2 &segment = *data->arr.originating_curves_begin(currentHalfEdge);
    //std::cout << "Half-edge   from: " << currentHalfEdge->source()->point() << " to " << currentHalfEdge->target()->point() << std::endl;
    //std::cout << "Source-edge from: " << segment.source() << " to " << segment.target() << std::endl;

    // These will always be sorted, it's how we created the segments
    const int aIndex = data->arrangementPointsIdices[segment.source()];
    const int bIndex = data->arrangementPointsIdices[segment.target()];

    // Sanity check
    assert(aIndex < bIndex);

    //printf("The original indices are %d and %d", data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]);
    //printf("\n");


    // Check to see if the segment and half edge have the same orientation
    const bool isSegmentLeftToRight = segment.source() < segment.target(); 
    const bool isCurrentHalfEdgeLeftToRight = (currentHalfEdge->direction() == CGAL::ARR_LEFT_TO_RIGHT);

    // The half edge has the same direction as the original edge
    if (isSegmentLeftToRight == isCurrentHalfEdgeLeftToRight)
    {
        //std::cout << "Upper link becomes lower link." << std::endl;

        //printf("The upper link of (%d, %d) : ", aIndex, bIndex);
        for (const int v: data->upperLink[std::pair<int, int>({aIndex, bIndex})]) 
        {
            assert(data->triangleToIndex.contains({aIndex, bIndex, v}));
            //minusTriangles.push_back({aIndex, bIndex, v});
            minusTriangles.insert(data->triangleToIndex[{aIndex, bIndex, v}]);

            //preimageGraphs[twinFaceID].erase({aIndex, bIndex, v});
            //printf("%d ", v);
        }
        //printf("\n");

        //printf("The lower link of (%d, %d) : ", aIndex, bIndex);
        for (const int v: data->lowerLink[std::pair<int, int>({aIndex, bIndex})]) 
        {
            assert(data->triangleToIndex.contains({aIndex, bIndex, v}));
            //plusTriangles.push_back({aIndex, bIndex, v});
            plusTriangles.insert(data->triangleToIndex[{aIndex, bIndex, v}]);
            //preimageGraphs[twinFaceID].insert({aIndex, bIndex, v});
            //printf("%d ", v);
        }
        //printf("\n");
    }
    else
    {
        //std::cout << "Lower link becomes upper link." << std::endl;

        //printf("The upper link of (%d, %d) : ", aIndex, bIndex);
        for (const int v: data->upperLink[std::pair<int, int>({aIndex, bIndex})]) 
        {
            assert(data->triangleToIndex.contains({aIndex, bIndex, v}));
            //plusTriangles.push_back({aIndex, bIndex, v});
            plusTriangles.insert(data->triangleToIndex[{aIndex, bIndex, v}]);
            //preimageGraphs[twinFaceID].insert({aIndex, bIndex, v});
            //printf("%d ", v);
        }
        //printf("\n");

        //printf("The lower link of (%d, %d) : ", aIndex, bIndex);
        for (const int v: data->lowerLink[std::pair<int, int>({aIndex, bIndex})]) 
        {
            assert(data->triangleToIndex.contains({aIndex, bIndex, v}));
            //minusTriangles.push_back({aIndex, bIndex, v});
            minusTriangles.insert(data->triangleToIndex[{aIndex, bIndex, v}]);
            //preimageGraphs[twinFaceID].erase({aIndex, bIndex, v});
            //printf("%d ", v);
        }
        //printf("\n");
    }

    return {minusTriangles, plusTriangles};
}

bool ReebSpace::isUpperLinkEdgeVertex(int aIndex, int bIndex, int vIndex, Data *data)
{
    // Make sure the vertices of the edge are in sorted order to have consistent orientation
    if (aIndex > bIndex)
    {
        std::swap(aIndex, bIndex);
    }

    // Define the two points that form the line
    const Point_2 a(data->vertexCoordinatesF[aIndex], data->vertexCoordinatesG[aIndex]);
    const Point_2 b(data->vertexCoordinatesF[bIndex], data->vertexCoordinatesG[bIndex]);

    // Define the test point
    const Point_2 v(data->vertexCoordinatesF[vIndex], data->vertexCoordinatesG[vIndex]);  // Change this to test different locations

    // Determine which half-plane r is in
    const CGAL::Orientation result = CGAL::orientation(a, b, v);

    //printf("Checking line (%d, %d) against vertex %d\n", aIndex, bIndex, vIndex);

    //std::cout << "a coords = " << a << std::endl;
    //std::cout << "b coords = " << b << std::endl;
    //std::cout << "v coords = " << v << std::endl;


    // Upper link = left
    if (result == CGAL::LEFT_TURN) {
        return true;
        //std::cout << "Point r is in the LEFT half-plane.\n";
    // Lower link = right
    } else if (result == CGAL::RIGHT_TURN) {
        return false;
        //std::cout << "Point r is in the RIGHT half-plane.\n";
    // This should not happen for generic maps
    } else {
        //std::cout << "Point r is on the line.\n";
        throw std::runtime_error("Input data is degeneate, a triangle is mapped to a line.");       
        //assert(false);
    }

    // Paranoid assert
    assert(false);
}


void ReebSpace::computeUpperLowerLink(Data *data)
{
    // All the edges from the mesh
    std::set<std::pair<int, int>> edges;

    // otherwise the lower and upper link flip around
    for (const std::vector<size_t> tet : data->tetrahedra)
    {
        // All pairs give you all six edges
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                // Get the indices of the vertices for the edge
                int aIndex = tet[a];
                int bIndex = tet[b];

                // Make sure the vertices of the edge are in sorted order to have consistent orientation
                if (aIndex > bIndex)
                {
                    std::swap(aIndex, bIndex);
                }

                // Add to the list of edge
                edges.insert({aIndex, bIndex});

                // Search though the other two unused vertices
                for (int v = 0 ; v < 4 ; v++)
                {
                    const int vIndex = tet[v];

                    // If the currect vertex is not one of two used to define the edge
                    if (vIndex != aIndex && vIndex != bIndex) 
                    {
                        const bool isUpperLink = isUpperLinkEdgeVertex(aIndex, bIndex, vIndex, data);

                        if (true == isUpperLink) {
                            data->upperLink[std::pair<int, int>({aIndex, bIndex})].insert(vIndex);
                        } else {
                            data->lowerLink[std::pair<int, int>({aIndex, bIndex})].insert(vIndex);
                        }
                    }
                }
            }
        }
    }

    // Initialize the disjoint set for the upper and lower link 
    std::map<std::pair<int, int>, DisjointSet<int>> upperLinkComponentsDS;
    for (const auto &[edge, vertices] : data->upperLink)
    {
        upperLinkComponentsDS[edge].initialize(vertices);
    }

    std::map<std::pair<int, int>, DisjointSet<int>> lowerLinkComponentsDS;
    for (const auto &[edge, vertices] : data->lowerLink)
    {
        lowerLinkComponentsDS[edge].initialize(vertices);
    }

    // Add all edges to compute the connected components

    for (const auto &[edge, vertices] : data->upperLink)
    {
        for (const auto &v1: vertices)
        {
            for (const auto &v2: vertices)
            {
                if ((v1 < v2 && edges.contains({v1, v2})) || (v1 > v2 && edges.contains({v2, v1})))
                {
                    upperLinkComponentsDS[edge].union_setsTriangle(v1, v2);
                }
            }
        }
    }

    for (const auto &[edge, vertices] : data->lowerLink)
    {
        for (const auto v1: vertices)
        {
            for (const auto v2: vertices)
            {
                if ((v1 < v2 && edges.contains({v1, v2})) || (v1 > v2 && edges.contains({v2, v1})))
                {
                    lowerLinkComponentsDS[edge].union_setsTriangle(v1, v2);
                }
            }
        }
    }

    for (const auto &edge : edges)
    {
        printf("Currently at the edge [%d, %d].\n", edge.first, edge.second);

        printf("The upper link had %d components and these vertices: ", upperLinkComponentsDS[edge].countConnectedComponents());
        for (const auto &v : data->upperLink[edge])
        {
            printf("%d ", v);
        }
        printf("\n");

        printf("The lower link had %d components and these vertices: ", lowerLinkComponentsDS[edge].countConnectedComponents());
        for (const auto &v : data->lowerLink[edge])
        {
            printf("%d ", v);
        }
        printf("\n");
    }
}


void ReebSpace::computeTriangleAdjacency(Data *data)
{

    std::set<std::set<int>> allTriangles;

    for (const std::vector<size_t> tet : data->tetrahedra)
    {
        // The triangles of the tet
        std::set<std::set<int>> triangles;

        // All pairs give you all six edges
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                for (int c = b + 1 ; c < 4 ; c++)
                {
                    int aIndex = tet[a];
                    int bIndex = tet[b];
                    int cIndex = tet[c];
                    allTriangles.insert({aIndex, bIndex, cIndex});
                }
            }
        }
    }

    //
    // Set up the indices for all triangles
    //

    for (const std::set<int> triangle : allTriangles)
    {
        data->indexToTriangle.push_back(triangle);
        data->triangleToIndex[triangle] = data->indexToTriangle.size() - 1;
    }

    data->adjacentTrianglesIndex.resize(allTriangles.size());



    // Compute the adjacency of triangles in the mesh, two triangles are adjacent when they are the faces of the same tet
    for (const std::vector<size_t> tet : data->tetrahedra)
    {
        // The triangles of the tet
        std::set<std::set<int>> triangles;

        // All pairs give you all six edges
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                for (int c = b + 1 ; c < 4 ; c++)
                {
                    int aIndex = tet[a];
                    int bIndex = tet[b];
                    int cIndex = tet[c];
                    triangles.insert({aIndex, bIndex, cIndex});
                }
            }
        }

        // Connect all the triangles together
        for(const std::set<int> t1 : triangles)
        {
            for(const std::set<int> t2 : triangles)
            {
                // Create a pair of sets
                std::pair<std::set<int>, std::set<int>> pairOfTriangles = {t1, t2};

                // Insert the pair into the set
                data->connectedTriangles.insert(pairOfTriangles);

                //data->adjacentTriangles[t1].push_back(t2);
                //data->adjacentTriangles[t2].push_back(t1);
                
                int t1Index = data->triangleToIndex[t1];
                int t2Index = data->triangleToIndex[t2];

                data->adjacentTrianglesIndex[t1Index].push_back(t2Index);
                data->adjacentTrianglesIndex[t2Index].push_back(t1Index);
            }
        }
    }
}

void ReebSpace::computeArrangement(Data *data) 
{
    Timer::start();

    // Add in the vertices of the mesh 
    data->arrangementPoints.resize(data->vertexCoordinatesF.size());
    for (int i = 0 ; i < data->vertexCoordinatesF.size() ; i++)
    {
        const float u = data->vertexCoordinatesF[i];
        const float v = data->vertexCoordinatesG[i];
        const Point_2 point(u, v);

        data->arrangementPoints[i] = point;
        data->arrangementPointsIdices[point] = i;
    };

    Timer::stop("Converted vertices to points           :");

    Timer::start();
    // Make sure you don't add duplicate edge to the arrangement
    std::map<std::set<size_t>, bool> uniqueEdges;
    for (int i = 0 ; i < data->tetrahedra.size() ; i++)
    {
        // Bottom face
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][0], data->tetrahedra[i][1]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][1], data->tetrahedra[i][2]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][2], data->tetrahedra[i][0]})] = true;

        // Connect bottom face to top vertex
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][0], data->tetrahedra[i][3]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][1], data->tetrahedra[i][3]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][2], data->tetrahedra[i][3]})] = true;
    }
    Timer::stop("Computed unique edges                  :");

    Timer::start();
    // Add the unique edges as setments to the arrangement
    std::vector<Segment_2> segments;
    for (const auto& edge : uniqueEdges) 
    {
        // Put in a vector for easy access
        std::vector<int> edgeVector(edge.first.begin(), edge.first.end());
        assert(edgeVector.size() == 2);

        segments.push_back(Segment_2(data->arrangementPoints[edgeVector[0]], data->arrangementPoints[edgeVector[1]]));
        //std::cout << "Adding edge " << edgeVector[0] << " - " << edgeVector[1] << std::endl;
    }
    Timer::stop("Converted edges to segments            :");



    //Timer::start();
    //int intersections = 0;
    //for (int i = 0 ; i < segments.size(); i++)
    //{
        //for (int j = i+1 ; j < segments.size(); j++)
        //{
            //// Check for intersection
            //if (CGAL::do_intersect(segments[i], segments[j])) 
            //{
                //intersections++;
            //}
        //}
    //}

    //std::cout << "There are " << intersections << " out of " << segments.size() << " segments." << std::endl;
    //Timer::stop("Checking intersections                 :");


    Timer::start();
    // Insert all the segments into the arrangement using the range insert method
    CGAL::insert(data->arr, segments.begin(), segments.end());
    Timer::stop("Computed arrangement                   :");

    //Timer::start();
    //Arrangement_2 arr;
    //for (const auto& segment : segments) {
        //CGAL::insert(arr, Curve_2(segment));
    //}
    //Timer::stop("Computed arrangement sequantially      :");


    std::cout << "The arrangement size:"
        << "   |V| = " << data->arr.number_of_vertices()
        << ",  |E| = " << data->arr.number_of_edges()
        << ",  |F| = " << data->arr.number_of_faces() << std::endl << std::endl;

    //std::cout << std::endl << std::endl << "The sequantial arrangement size:\n"
        //<< "   |V| = " << arr.number_of_vertices()
        //<< ",  |E| = " << arr.number_of_edges()
        //<< ",  |F| = " << arr.number_of_faces() << std::endl;



    // Print out all the faces in the arrangement
    //std::cout << "Faces in the arrangement:" << std::endl;

    int counter = 0;

    data->arrangementIndexToFace.resize(data->arr.number_of_faces());

    for (auto f = data->arr.faces_begin(); f != data->arr.faces_end(); ++f) 
    {

        Arrangement_2::Face_const_handle a = f;
        data->arrangementFacesIdices[a] = counter;
        data->arrangementIndexToFace[counter] = f;
        counter++;

        //std::cout << data->arrangementFacesIdices[a] << std::endl;

        //if (f->is_unbounded()) {
            //std::cout << "Unbounded face" << std::endl;
            //continue;
        //}

        //std::cout << "Bounded face with " << f->number_of_holes() << " holes " << std::endl;


        //std::cout << "inner     = " << f->number_of_inner_ccbs() << std::endl;
        //std::cout << "outer     = " << f->number_of_outer_ccbs() << std::endl;
        //std::cout << "holes     = " << f->number_of_holes() << std::endl;
        //std::cout << "isolated  = " << f->number_of_isolated_vertices() << std::endl;

        //print_ccb<Arrangement_2>(f->outer_ccb());


        //f->number_of_inner_ccbs

        //for (auto e = f->outer_ccbs_begin(); e != f->outer_ccbs_end(); ++e) {
        //std::cout << "E" << std::endl;



        //}
    }

    


    // Use to find reverse edges
    //for(auto e = arr.edges_begin(); e != arr.edges_end(); ++e) {
        //std::cout << "Edge: " << e->curve() << "\n";
        //std::cout << "  Originating curve(s):\n";
        //for(auto oc = arr.originating_curves_begin(e);
                //oc != arr.originating_curves_end(e); ++oc)
        //{
            //std::cout << "    " << *oc << "\n";
        //}
        //std::cout << std::endl;
    //}
}

void ReebSpace::testTraverseArrangement(Data *data)
{
    // Find the unbounded face (hold the boundary of the arrangement)
    Face_const_handle outerFace;

    // Iterate over all faces and find the unbounded one
    for (Face_const_iterator fit = data->arr.faces_begin(); fit != data->arr.faces_end(); ++fit) 
    {
        if (fit->is_unbounded()) 
        {
            // If this has already been set, we have two outerFace, this should never happen.
            assert(outerFace == Face_const_handle());
            outerFace = fit;
            break;
        }
    }

    // Sanity check, make sure the outer face is simple
    assert(outerFace->number_of_inner_ccbs() == 1);
    assert(outerFace->number_of_outer_ccbs() == 0);
    assert(outerFace->number_of_holes() == 1);
    assert(outerFace->number_of_isolated_vertices() == 0);

    // Get the the first half edge of the outerFace (could be any edge, this is a matter of convention)
    Halfedge_const_handle outerHalfEdge = *outerFace->holes_begin();

    // Make sure there is only one originating curve, something has gone wrong otherwise (edge overlap)
    assert(std::distance(data->arr.originating_curves_begin(outerHalfEdge), data->arr.originating_curves_end(outerHalfEdge)) == 1);

    // Starting from an outer half edge
    std::queue<Arrangement_2::Halfedge_const_handle> traversalQueue;
    traversalQueue.push(outerHalfEdge);

    // Which faces are visited
    std::set<Arrangement_2::Face_const_handle> visited;
    visited.insert(outerHalfEdge->face());

    while (false == traversalQueue.empty())
    {
        // Pop an half edge out
        Arrangement_2::Halfedge_const_handle currentHalfEdge = traversalQueue.front();
        traversalQueue.pop();
        Arrangement_2::Face_const_handle currentFace = currentHalfEdge->face();

        // Get the twin
        Arrangement_2::Halfedge_const_handle twin = currentHalfEdge->twin();
        Arrangement_2::Face_const_handle twinFace = twin->face();

        // Get ids of the current face and the twin face
        int currentFaceID = data->arrangementFacesIdices[currentFace];
        int twinFaceID = data->arrangementFacesIdices[twinFace];

        // If we have never visited this face, then we have never visited any of the half edges.
        if (visited.find(twinFace) == visited.end())
        {
            //printf("NEW FACE ------------------------------------------ %d -> %d \n", data->arrangementFacesIdices[currentFace], data->arrangementFacesIdices[twinFace]);
            visited.insert(twinFace);

            // Sanity check we should only have onbounded face, the outside face.
            assert(false == twinFace->is_unbounded());


            Arrangement_2::Ccb_halfedge_const_circulator start = twinFace->outer_ccb();
            Arrangement_2::Ccb_halfedge_const_circulator curr = start;

            do {
                traversalQueue.push(curr);

                // Make sure there is only one originating curve (sanity check)
                //const Segment_2 &segment = *data->arr.originating_curves_begin(curr);
                //std::cout << "Half-edge   from: " << curr->source()->point() << " to " << curr->target()->point() << std::endl;
                //std::cout << "Source-edge from: " << segment.source() << " to " << segment.target() << std::endl;
                //printf("The original indices are %d and %d", data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]);
                //printf("\n\n");

                ++curr;
            } while (curr != start);
        }
    }
}





void ReebSpace::computeTwinFacePreimageGraph(Data *data, Arrangement_2::Halfedge_const_handle &currentHalfEdge)
{
    // Get the face
    Arrangement_2::Face_const_handle currentFace = currentHalfEdge->face();

    // Get the twin
    Arrangement_2::Halfedge_const_handle twin = currentHalfEdge->twin();
    Arrangement_2::Face_const_handle twinFace = twin->face();

    // Get ids of the current face and the twin face
    int currentFaceID = data->arrangementFacesIdices[currentFace];
    int twinFaceID = data->arrangementFacesIdices[twinFace];

    auto [minusTriangles, plusTriangles] = ReebSpace::getMinusPlusTrianglesIndex(currentHalfEdge, data);

    //std::cout << "Computed minus/plus triangles..." << std::endl;
    //printf("The current preimage graph has %d triangles...\n", data->preimageGraphs[currentFaceID].data.size());

    // Set the current preimage graph to be the preimage graph of the parent
    std::set<int> preimageGraph;
    for (const auto &[t, id] : data->preimageGraphs[currentFaceID].data)
    {
        preimageGraph.insert(t);
    }

    //printf("Minus triangles...\n");

    // Add and remove triangles of the upper/lower link of the crossed edge
    for (const auto &triangle: minusTriangles)
    {
        preimageGraph.erase(triangle);
        //std::cout << triangle << std::endl;
    }

    //printf("Plus triangles...\n");
    for (const auto &triangle: plusTriangles)
    {
        preimageGraph.insert(triangle);
        //std::cout << triangle << std::endl;
    }

    //std::cout << "Computed preimage graph soup..." << std::endl;
    //for (const auto &triangle: preimageGraph)
    //{
        //std::cout << triangle << std::endl;
    //}

    //
    // Step 3. Compute disjointSets[twinFaceID]
    //

    data->preimageGraphs[twinFaceID].initialize(preimageGraph);

    for (const auto &[t1, id1] : data->preimageGraphs[twinFaceID].data)
    {
        for (const auto &t2 : data->adjacentTrianglesIndex[t1])
        {
            if (data->preimageGraphs[twinFaceID].data.contains(t2))
            {
                data->preimageGraphs[twinFaceID].union_setsTriangle(t1, t2);
            }
        }
    }

    //std::cout << "Unioned sets..." << std::endl;

    // Finaly make sure everyon points to their root
    data->preimageGraphs[twinFaceID].update();

    // Used when drawing the arrangement
    //data->arrangementFiberComponents[twinFaceID] = data->preimageGraphs[twinFaceID].countConnectedComponents();
}

void ReebSpace::computePreimageGraphs(Data *data, const bool discardFiberSeedsSets)
{

    using namespace indicators;
    ProgressBar bar{
        option::BarWidth{50},
            option::Start{"["},
            option::Fill{"■"},
            option::Lead{"■"},
            option::Remainder{"-"},
            option::End{" ]"},
            option::PostfixText{"Computing Reeb space."},
            option::ShowPercentage{true},  // Show percentage on the bar
            option::ShowElapsedTime{true},
            //option::ForegroundColor{Color::cyan},
            //option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
    };









    // Find the unbounded face (hold the boundary of the arrangement)
    Face_const_handle outerFace;

    // Iterate over all faces and find the unbounded one
    for (Face_const_iterator fit = data->arr.faces_begin(); fit != data->arr.faces_end(); ++fit) 
    {
        if (fit->is_unbounded()) 
        {
            // If this has already been set, we have two outerFace, this should never happen.
            assert(outerFace == Face_const_handle());
            outerFace = fit;
            break;
        }
    }

    // Sanity check, make sure the outer face is simple
    assert(outerFace->number_of_inner_ccbs() == 1);
    assert(outerFace->number_of_outer_ccbs() == 0);
    assert(outerFace->number_of_holes() == 1);
    assert(outerFace->number_of_isolated_vertices() == 0);

    // Loop the edges of the boundary
    //Arrangement_2::Ccb_halfedge_const_circulator inner_ccb = *outerFace->holes_begin();
    //Arrangement_2::Ccb_halfedge_const_circulator he = inner_ccb;

    // Iterate over the halfedges forming the inner CCB
    //do {
        //std::cout << he->source()->point() << " -> ";
        //std::cout << he->target()->point() << std::endl;

        // Get the originating curve of the half edge
        // Maybe need a function for this
        //const Segment_2& segment = *data->arr.originating_curves_begin(he);
        //std::cout << data->arrangementPointsIdices[segment.source()] << ", ";
        //std::cout << data->arrangementPointsIdices[segment.target()] << std::endl;
    //} while (++he != inner_ccb);

    //printf("\n\n");

    // Get the the first half edge of the outerFace (could be any edge, this is a matter of convention)
    //Halfedge_const_handle outerHalfEdge = *outerFace->holes_begin();

    // Make sure there is only one originating curve, something has gone wrong otherwise (edge overlap)
    //assert(std::distance(data->arr.originating_curves_begin(outerHalfEdge), data->arr.originating_curves_end(outerHalfEdge)) == 1);

    // Extract the original curve
    //const Segment_2& segment = *data->arr.originating_curves_begin(outerHalfEdge);

    //printf("The half edge came from this edge:\n");
    //std::cout << segment.source() << " -> " << segment.target() << std::endl;

    //printf("The original indices are %d and %d", data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]);
    //printf("\n\n");







    // Starting from an outer half edge
    std::queue<Face_const_handle> traversalQueue;
    // This is the order in which faces are added, also servers as a visited flag
    std::vector<int> order(data->arrangementFacesIdices.size(), -1);


    traversalQueue.push(outerFace);

    // This is the order in which a face has been processed, note this is different than level
    // This is used as an index for faces, so that we don't do double work when checking edges for correspondence
    int orderIndex = 0;
    order[data->arrangementFacesIdices[outerFace]] = orderIndex;

    // The disjoint set to track the connected components of the preimage graph
    data->preimageGraphs.resize(data->arrangementFacesIdices.size());

    // The number of connected components for each preimage graph (computed from the disjoint set)
    //data->arrangementFiberComponents.resize(data->arrangementFacesIdices.size(), -1);

    // If we want the fiber seeds, initialize them
    if (false == discardFiberSeedsSets)
    {
        data->fiberSeeds.resize(data->arrangementFacesIdices.size());
    }

    int graphsInMemory = 0;
    float averageAraphsInMemory = 0;
    int barTickThreshold = data->arrangementFacesIdices.size() / 50;

    while (false == traversalQueue.empty())
    {
        // Pop an half edge out
        Face_const_handle currentFace = traversalQueue.front();
        traversalQueue.pop();

        // Get ids of the current face and the twin face
        int currentFaceID = data->arrangementFacesIdices[currentFace];

        //std::cout << "Current face " << currentFaceID << std::endl;

        // Sanity check, this should always be true
        if (false == currentFace->is_unbounded()) 
        {
            //assert(false == data->preimageGraphs[currentFaceID].isEmpty());
        }

        //
        // For each neighbouring face
        //

        // Unbounded face (the starting one) has a different way of addressing its neighbours
        Halfedge_const_handle start;
        if (currentFace->is_unbounded()) 
        {  
            start = *currentFace->holes_begin();
        }
        else
        {
            start = *currentFace->outer_ccbs_begin();
        }

        Arrangement_2::Ccb_halfedge_const_circulator curr = start;
        do {

            Face_const_handle twinFace = curr->twin()->face();
            const int twinFaceID = data->arrangementFacesIdices[twinFace];

            // If the neighbour has not been visited, we enqueue it and also compute its preimage graph
            if (-1 == order[twinFaceID])
            {
                traversalQueue.push(twinFace);
                order[twinFaceID] = ++orderIndex;

                //std::cout << "Computing for neighbour " << twinFaceID << std::endl;

                // Compute the preimage graph of this unvisited face
                ReebSpace::computeTwinFacePreimageGraph(data, curr);
                graphsInMemory++;

                // Get all unique roots and a representative
                const std::vector<std::pair<int, int>> representativesAndRoots = data->preimageGraphs[twinFaceID].getUniqueRepresentativesAndRoots();

                // Initialize the vertices of H with the connected components of the twin graph
                for (const auto &[representative, root] : representativesAndRoots)
                {
                    data->reebSpace.addElements({twinFaceID, root});
                }

                // If we want fiber computation, cache the seeds for the fiber components
                if (false == discardFiberSeedsSets)
                {
                    data->fiberSeeds[twinFaceID] = representativesAndRoots;
                }
            }

            // Sanity check, all graphs should have either been computed before or now
            //assert(false == data->preimageGraphs[twinFaceID].isEmpty());

            // Compute the correspondence with the neighbours, but only if they are at a higher level, or we are at the same level, currentFaceID < twinFaceID is used to avoid double work, we only need it once
            if (order[currentFaceID] < order[twinFaceID])
            {
                //std::cout << "Determining correspondence with neighbour " << twinFaceID << std::endl;
                ReebSpace::determineCorrespondence(data, curr);
            }

            ++curr;
        } while (curr != start);

        averageAraphsInMemory = averageAraphsInMemory + ((float)graphsInMemory - (float)averageAraphsInMemory) / (float)orderIndex;

        //printf("There are %d active preimage graphs with average %f at index %d/%ld.\n", graphsInMemory, averageAraphsInMemory, orderIndex, data->preimageGraphs.size());

        // If the threshold is zero the ticks bugs out, so we don't do it
        if (barTickThreshold > 0 && order[currentFaceID] % barTickThreshold == 0)
        {
            // Update bar state
            if (false == bar.is_completed()) {  bar.tick(); }
            if (false == bar.is_completed()) {  bar.tick(); }
        }


        // Dispose of the preimage graph we will no longer need it
        data->preimageGraphs[currentFaceID].clear();
        graphsInMemory--;
    }

    bar.set_progress(100); // all done
    printf("\n\nThere is an average of %f / %ld active preimage graphs.\n", averageAraphsInMemory, data->preimageGraphs.size());
    printf("The correspondence graphs has %ld vertices and the Reeb space has %ld sheets.\n\n", data->reebSpace.data.size(), data->reebSpace.getUniqueRepresentativesAndRoots().size());
}


void ReebSpace::determineCorrespondence(Data *data, Arrangement_2::Halfedge_const_handle &halfEdge)
{
    // Get the face IDs of the current face and its twin
    Face_const_handle face = halfEdge->face();
    Face_const_handle twinFace = halfEdge->twin()->face();

    const int faceID = data->arrangementFacesIdices[face];
    const int twinFaceID = data->arrangementFacesIdices[twinFace];

    // Get the originating edge
    const Segment_2 &segment = *data->arr.originating_curves_begin(halfEdge);
    const std::pair<int, int> originatingEdge = {data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]};

    // The triangls that are added/removd from face -> twinFace
    const auto& [minusTriangles, plusTriangles] = ReebSpace::getMinusPlusTrianglesIndex(halfEdge, data);


    // See which of the roots (connected components) are active (contain an active triangle)
    std::set<int> activeRootsFace;
    for (int triangle : minusTriangles)
    {
        activeRootsFace.insert(data->preimageGraphs[faceID].findTriangle(triangle));
    }

    std::set<int> activeRootsTwinFace;
    for (int triangle : plusTriangles)
    {
        activeRootsTwinFace.insert(data->preimageGraphs[twinFaceID].findTriangle(triangle));

    }

    // Definite edge
    if (activeRootsFace.size() == 0 || activeRootsTwinFace.size() == 0)
    {

        data->jacobiType[originatingEdge] = 0;
    }
    // Reeb-regular
    else if (activeRootsFace.size() == 1 && activeRootsTwinFace.size() == 1)
    {

        data->jacobiType[originatingEdge] = 1;
    }
    // Indefinite edge
    else
    {
        data->jacobiType[originatingEdge] = 2;

    }

    // This is a regular edge then connect the two active components, for singular edge we skip this, they are considered new
    if (activeRootsFace.size() == 1 && activeRootsTwinFace.size() == 1)
    {
        // Active roots point to each other
        // Add edge [faceId, activeRootsFace.begin()->second()] -> [twinFaceId, activeRootsTwinFace.begin()->second()]
        //connectedFacesAndRoots
        //std::set<int> faceRoot({faceID, *activeRootsFace.begin()->second(});
        //std::pair<std::set<int>, std::set<int>> reebSpaceConnection({t1, t2});


        // @TODO Ulgy, but hard to get the first element out otherwise
        for (const int &rFace :activeRootsFace)
        {
            for (const int &rTwinFace :activeRootsTwinFace)
            {
                //connectedFacesAndRoots.insert(reebSpaceConnection);
                //data->edgesH.push_back({
                        //{faceID, rFace},
                        //{twinFaceID, rTwinFace}
                        //});


                data->reebSpace.union_setsTriangle({faceID, rFace}, {twinFaceID, rTwinFace});



                //int faceVertexHindex = data->vertexHtoIndex[{faceID, rFace}];
                //int twinFaceVertexHindex = data->vertexHtoIndex[{twinFaceID, rTwinFace}];

                    //data->indexToVertexH.push_back(vertexH);
                    //data->vertexHtoIndex[vertexH] = verexHIndex;


                //data->reebSpace.union_setsTriangle(faceVertexHindex, twinFaceVertexHindex);

            }
        }
    }

    //
    // Link together all other connected components
    //
    for (const auto &[t, id] : data->preimageGraphs[faceID].data)
    {
        // The root of the triangle in the face
        const int triangleRootFace = data->preimageGraphs[faceID].findTriangle(t);

        // We have already deal with the active fiber
        if (activeRootsFace.contains(triangleRootFace)) { continue; }

        // The root of the triangle in the twin face
        const int triangleRootTwinFace = data->preimageGraphs[twinFaceID].findTriangle(t);

        //connectedFacesAndRoots.insert(reebSpaceConnection);
        //data->edgesH.push_back({
                //{faceID, triangleRootFace}, 
                //{twinFaceID, triangleRootTwinFace}
                //});
        data->reebSpace.union_setsTriangle({faceID, triangleRootFace}, {twinFaceID, triangleRootTwinFace});

        //int faceVertexHindex = data->vertexHtoIndex[{faceID, triangleRootFace}];
        //int twinFaceVertexHindex = data->vertexHtoIndex[{twinFaceID, triangleRootTwinFace}];
        //data->reebSpace.union_setsTriangle(faceVertexHindex, twinFaceVertexHindex);
    }
}




void ReebSpace::computeCorrespondenceGraph(Data *data)
{

    // For evey face
    for (auto face = data->arr.faces_begin(); face != data->arr.faces_end(); ++face) 
    {
        // Skip the outer face
        if (face->is_unbounded()) { continue; }


        // Walk around the boundary of the face
        Arrangement_2::Ccb_halfedge_const_circulator start = face->outer_ccb();
        Arrangement_2::Ccb_halfedge_const_circulator curr = start;
        do {

            ReebSpace::determineCorrespondence(data, curr);
            ++curr;
        } while (curr != start);

    }
}



void ReebSpace::computeReebSpacePostprocess(Data *data)
{

    Timer::start();
    // Path compression to make sure every element of H is poiting to its root
    data->reebSpace.update();
    Timer::stop("Updating RS disjoint set               :");


    // For faster lookups cache the sheetd IDs for each face
    std::vector<std::unordered_set<int>> faceSheets(data->fiberSeeds.size());
    for (int i = 0 ; i < data->fiberSeeds.size() ; i++)
    {
        for (const auto &[triangleId, fiberComponentId] : data->fiberSeeds[i])
        {
            const int sheetId = data->reebSpace.findTriangle({i, fiberComponentId});
            faceSheets[i].insert(sheetId);
        }

    }

    Timer::start();
    //
    // In order to compute the polygon of each sheet, first obtain a halfEdge of the arrangement that is on the boundary of the sheet
    //
    std::unordered_map<int, Arrangement_2::Halfedge_const_handle> sheetSeeds;
    for (auto currentFaceIterator = data->arr.faces_begin(); currentFaceIterator != data->arr.faces_end(); ++currentFaceIterator) 
    {
        Arrangement_2::Face_const_handle currentFace = currentFaceIterator;
        if (currentFace->is_unbounded()) { continue; }

        int currentFaceID = data->arrangementFacesIdices[currentFace];

        //printf("\n\nFace %d has these sheets - ", currentFaceID);

        // Which sheets contain the current face
        std::unordered_set<int> currentFaceSheetIds = faceSheets[currentFaceID];
        //for (const auto &[triangleId, fiberComponentId] : data->fiberSeeds[currentFaceID])
        //{
            //const int sheetId = data->reebSpace.findTriangle({currentFaceID, fiberComponentId});
            //currentFaceSheetIds.push_back(sheetId);
            ////printf("%d ", sheetId);
        //}
        //printf("\n");


        Arrangement_2::Ccb_halfedge_const_circulator start = currentFace->outer_ccb();
        Arrangement_2::Ccb_halfedge_const_circulator curr = start;
        do {
            Face_const_handle twinFace = curr->twin()->face();
            if (twinFace->is_unbounded()) { curr++; continue; }

            const int twinFaceID = data->arrangementFacesIdices[twinFace];
            std::unordered_set<int> twinFaceSheetIds = faceSheets[twinFaceID];

            // Which sheets are in the currentFace, but NOT in the twin face
            std::vector<int> diff;
            for (const int &sheetId : currentFaceSheetIds) 
            {
                if (false == twinFaceSheetIds.contains(sheetId))
                {
                    diff.push_back(sheetId);
                }
            }

            for (const int &sheetId : diff)
            {
                //std::cout << "Half-edge   from: " << curr->source()->point() << " to " << curr->target()->point() << " is a seed for sheet " << sheetId << std::endl;

                if (false == sheetSeeds.contains(sheetId))
                {
                    sheetSeeds[sheetId] = curr;
                }
            }

            ++curr;
        } while (curr != start);
    }
    Timer::stop("Computing sheet seeds                  :");





    //
    // At this point we have a half edge seed for every sheet, we use it to traverse the boundary of the sheet to obtain its geometry
    // We rely on the half edge data structure provided by the arrangement
    // if ab is an edge between vertices a and b in A and, then we can get the next edge bc as the one of the neighours of b that has
    // the sheet in its face and it doesn't have the sheet in the twin face
    //


    //printf("Done with sheet seeds\n");


    Timer::start();
    for (const auto &[sheetId, halfEdge]: sheetSeeds)
    {
        // Sanity check, make sure the seed we have picked actually is on the boundary of the sheet
        Face_const_handle faceSeedEdge = halfEdge->face();
        const int faceSeedEdgeId = data->arrangementFacesIdices[faceSeedEdge];
        if (false == faceSheets[faceSeedEdgeId].contains(sheetId))
        {
            throw std::runtime_error("The face of the seed edge does is not in the sheet.");
        }

        Face_const_handle faceSeedTwinEdge = halfEdge->twin()->face();
        const int faceSeedTwinEdgeId = data->arrangementFacesIdices[faceSeedTwinEdge];
        if (true == faceSheets[faceSeedTwinEdgeId].contains(sheetId))
        {
            throw std::runtime_error("The face of the twin of the seed edge does is no in the sheet.");
        }


        //printf("---------------------------------------------------------------------------------------------       At sheet %d \n", sheetId);

        // Starting location, we don't make it visited on purpose, we want to discover it latex to finish the loop
        Vertex_const_handle startVertex = halfEdge->source();

        // Time to find the next vertex
        Vertex_const_handle currentVertex = startVertex;
        std::unordered_set<Vertex_const_handle> visited;
        do 
        {

            // Add the current vertex to the polygon
            data->sheetPolygon[sheetId].push_back({
                    CGAL::to_double(currentVertex->point().x()), 
                    CGAL::to_double(currentVertex->point().y())
                    });


            //printf("-------------------------- The current vertex (%f, %f) \n ", CGAL::to_double(currentVertex->point().x()), CGAL::to_double(currentVertex->point().y()));



            //
            // Find the next vertex on the boundary of the sheet, it will be adjacent to the current vertex
            //
            int sheetBoundaryNeighbours = 0;
            int allNeighbours = 0;
            visited.insert(currentVertex);

            // This loop should always finish
            const auto begin = currentVertex->incident_halfedges();
            auto circ = begin;
            do {

                allNeighbours++;

                // Too many neighourss
                if (allNeighbours > data->arr.number_of_vertices())
                {
                    break;
                }

                const auto twinHalfEdge = circ->twin();
                const auto nextVertex = twinHalfEdge->target();
                //printf("Looking at neighbour (%f, %f) \n ", CGAL::to_double(twinHalfEdge->target()->point().x()), CGAL::to_double(twinHalfEdge->target()->point().y()));

                Face_const_handle faceA = circ->face();
                const int faceAId = data->arrangementFacesIdices[faceA];

                Face_const_handle faceB = twinHalfEdge->face();
                const int faceBId = data->arrangementFacesIdices[faceB];

                // If one of the face contains the sheet
                const bool faceHalfEdgeContainsSheet = faceSheets[faceAId].contains(sheetId);
                const bool faceTwinHalfEdgeContainsEdge = faceSheets[faceBId].contains(sheetId);

                if (faceHalfEdgeContainsSheet == true && faceTwinHalfEdgeContainsEdge == false)
                {
                    currentVertex = nextVertex;
                    sheetBoundaryNeighbours++;
                }

                ++circ;
            } while (circ != begin);

            // A the polygon of a sheet sheet cannot have more vertices than the number of the vertices in the arrangement, something went wrong.
            if (data->sheetPolygon[sheetId].size() > data->arr.number_of_vertices())
            {
                printf("The boundary of sheet %d is degenerate, it has more vertices %ld than the arrangement.\n", sheetId, data->sheetPolygon[sheetId].size());
                data->incompleteSheets.insert(sheetId);
                break;
            }
            
            if (sheetBoundaryNeighbours !=1)
            {
                printf("The boundary of sheet %d is degenerate, a vertex has %d boundary neighours, should only be one.\n", sheetId, sheetBoundaryNeighbours);
                data->incompleteSheets.insert(sheetId);
                break;
            }

            // I'm not sure if this will ever get triggered, since if we loop back to someone else than the start 
            // that should mean 
            if (currentVertex != startVertex && visited.contains(currentVertex))
            {
                printf("The boundary of sheet %d is degenerate, a vertex loops back to an already visited one.\n", sheetId);
                data->incompleteSheets.insert(sheetId);
                break;

            }

            if (allNeighbours > data->arr.number_of_vertices())
            {
                printf("The boundary of sheet %d is degenerate, the inner loop kept going for too long.\n", sheetId);
                data->incompleteSheets.insert(sheetId);
                break;
            }

        } while (currentVertex != startVertex);
    }
    Timer::stop("Computing sheet boundary polygons      :");

    if (data->incompleteSheets.size() == 0)
    {
        std::cout << "\nThe boundaries of all sheets are okay! No degeneracy.\n" << std::endl;
    }


    // Compute the areas of sheets


    Timer::start();
    for (const auto &[sheetId, polygon] : data->sheetPolygon)
    {
        float area = 0.0;

        // If the sheet is incomplete sum up the faces that make it
        if (data->incompleteSheets.contains(sheetId))
        {
            // Loop through all faces to see which ones are in the sheet
            for (auto f = data->arr.faces_begin(); f != data->arr.faces_end(); ++f) 
            {
                const int currentFaceID = data->arrangementFacesIdices[f];

                // For each fiber component in the face, see if one of those is in our sheet
                for (const auto &[triangleId, fiberComponentId] : data->fiberSeeds[currentFaceID])
                {
                    const int componentSheetId = data->reebSpace.findTriangle({currentFaceID, fiberComponentId});

                    // Now we can add the polygon
                    if (componentSheetId == sheetId)
                    {
                        CartesianPolygon_2 poly;

                        typename Arrangement_2::Ccb_halfedge_const_circulator circ = f->outer_ccb();
                        typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
                        do {
                            typename Arrangement_2::Halfedge_const_handle e = curr;

                            // Get point from CGAL (and convert to double )
                            const float u = CGAL::to_double(e->source()->point().x());
                            const float v = CGAL::to_double(e->source()->point().y());

                            poly.push_back({u, v});
                        } while (++curr != circ);

                        area += poly.area();
                    }
                }
            }
            printf("Incomplete sheet %d has area %f\n", sheetId, area);
        }
        else
        {
            area = abs(polygon.area());
        }

        data->sheetArea[sheetId] = area;
        
        //printf("The polygon of sheet %d has size %ld and area %f\n", sheetId, polygon.size(), area);
        //printf("The area of sheet %d is %f \n", sheetId, area);
    }
    Timer::stop("Computing sheet areas                  :");
    //Timer::stop("Computing sheet boundary polygons      :");

    Timer::start();
    // Transfer the map entries into a vector of pairs so that we can sort
    std::vector<std::pair<int, double>> sheetAreaSortVector(data->sheetArea.begin(), data->sheetArea.end());

    //printf("The sort vector has size %ld\n", sheetAreaSortVector.size());

    // Sort by area
    std::sort(sheetAreaSortVector.begin(), sheetAreaSortVector.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
            return a.second > b.second; 
            });

    // Assign a sorted order to each sheet
    for (int i = 0 ; i < sheetAreaSortVector.size() ; i++)
    {
        const auto &[sheetId, area] = sheetAreaSortVector[i];
        data->sheetToColour[sheetId] = i;
    }
    Timer::stop("Sorting sheets and labeling them       :");

    std::cout << "\nHere are the top 20 sheets sorted by range area.\n";
    // Print to debug, at least the first new
    //for (int i = 0 ; i < sheetAreaSortVector.size() ; i++)
    const int maxSheets = std::min((size_t)30, sheetAreaSortVector.size());
    for (int i = 0 ; i < maxSheets ; i++)
    {
        std::cout << i << " -- sheet " << sheetAreaSortVector[i].first << " has area " << sheetAreaSortVector[i].second << std::endl;
    }



}

void ReebSpace::countIntersectionsTypes(Data *data)
{
    int regularRegularIntersections = 0;
    int singularRegularIntersections = 0;
    int singularSingularIntersections = 0;
    int degenerateIntersections = 0;
    int maxNeihbours = 0;

    // Now loop over all vertices
    for (auto currentVertex = data->arr.vertices_begin(); currentVertex != data->arr.vertices_end(); ++currentVertex)
    {
        const Point_2& p = currentVertex->point();
        //std::cout << "Vertex at: " << p << std::endl;

        // Skip the vertices, we just want to count the inter intersections.
        if (data->arrangementPointsIdices.contains(p))
        {
            continue;
        }


        // This loop should always finish
        const auto begin = currentVertex->incident_halfedges();
        auto circ = begin;

        int regularNeighbours = 0;
        int singularNeighbours = 0;

        do {

            const Segment_2 &segment = *data->arr.originating_curves_begin(circ);

            // These will always be sorted, it's how we created the segments
            const int aIndex = data->arrangementPointsIdices[segment.source()];
            const int bIndex = data->arrangementPointsIdices[segment.target()];

            if (data->jacobiType[{aIndex, bIndex}] == 1)
            {
                regularNeighbours++;
            }
            else
            {
                singularNeighbours++;
            }

            ++circ;
        } while (circ != begin);

        if (regularNeighbours + singularNeighbours > 4)
        {
            degenerateIntersections++;
        }

        maxNeihbours = std::max(maxNeihbours, regularNeighbours + singularNeighbours);



        if (regularNeighbours == 0)
        {
            singularSingularIntersections++;
        }
        else if (singularNeighbours == 0)
        {
            regularRegularIntersections++;
        }
        else
        {
            singularRegularIntersections++;
        }
    }


    const int totalIntersctionPoints = regularRegularIntersections + singularRegularIntersections + singularSingularIntersections;

    printf("Here is a summary of the intersections types:\n");
    printf("regular   - regular  intersections: %d the ratio is %.2f%%.\n", regularRegularIntersections, 100.0 * (float)regularRegularIntersections / (float)totalIntersctionPoints);
    printf("singular  - regular  intersections: %d the ratio is %.2f%%.\n", singularRegularIntersections, 100.0 * (float)singularRegularIntersections / (float)totalIntersctionPoints);
    printf("singular  - sigular  intersections: %d the ratio is %.2f%%.\n", singularSingularIntersections, 100.0 * (float)singularSingularIntersections / (float)totalIntersctionPoints);
    printf("The number of degenerate intersections is %d with maximum being %d.\n", degenerateIntersections, maxNeihbours);

}
