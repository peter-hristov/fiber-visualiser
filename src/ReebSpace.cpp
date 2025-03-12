#include "./ReebSpace.h"
#include "./DisjointSet.h"
#include "./CGALTypedefs.h"
#include "./Timer.h"

#include <cstddef>
#include <map>
#include <queue>
#include <set>


std::pair<std::vector<std::set<int>>, std::vector<std::set<int>>> ReebSpace::getMinusPlusTriangles(Arrangement_2::Halfedge_const_handle currentHalfEdge, Data *data)
{
    // Step 1. Initialize lists
    std::vector<std::set<int>> plusTriangles;
    std::vector<std::set<int>> minusTriangles;

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
            minusTriangles.push_back({aIndex, bIndex, v});
            //preimageGraphs[twinFaceID].erase({aIndex, bIndex, v});
            //printf("%d ", v);
        }
        //printf("\n");

        //printf("The lower link of (%d, %d) : ", aIndex, bIndex);
        for (const int v: data->lowerLink[std::pair<int, int>({aIndex, bIndex})]) 
        {
            plusTriangles.push_back({aIndex, bIndex, v});
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
            plusTriangles.push_back({aIndex, bIndex, v});
            //preimageGraphs[twinFaceID].insert({aIndex, bIndex, v});
            //printf("%d ", v);
        }
        //printf("\n");

        //printf("The lower link of (%d, %d) : ", aIndex, bIndex);
        for (const int v: data->lowerLink[std::pair<int, int>({aIndex, bIndex})]) 
        {
            //preimageGraphs[twinFaceID].erase({aIndex, bIndex, v});
            minusTriangles.push_back({aIndex, bIndex, v});
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
        assert(false);
    }

    // Paranoid assert
    assert(false);
}


void ReebSpace::computeUpperLowerLink(Data *data)
{
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

                // Search though the other two unused vertices
                for (int v = 0 ; v < 4 ; v++)
                {
                    int vIndex = tet[v];

                    // If the currect vertex is not one of two used to define the edge
                    if (vIndex != aIndex && vIndex != bIndex) 
                    {
                        bool isUpperLink = isUpperLinkEdgeVertex(aIndex, bIndex, vIndex, data);

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


    // If the upper link is empty we have a definite edge
    for (const auto &[edge, linkVertices] : data->upperLink)
    {
        //printf("Upper link of (%d, %d) is:\n", edge.first, edge.second);
        //for (const auto &v : linkVertices)
        //{
            //printf("(%d, %d, %d)\n", edge.first, edge.second, v);
        //}

        // If the upper link is full, but the lower link is empty
        if (false == data->lowerLink.contains(edge))
        {
            data->jacobiType[edge] = 0;
        }

    }

    for (const auto &[edge, linkVertices] : data->lowerLink)
    {
        //printf("Lower link of (%d, %d) is:\n", edge.first, edge.second);

        //for (const auto &v : linkVertices)
        //{
            //printf("(%d, %d, %d)\n", edge.first, edge.second, v);
        //}

        // If the lower link is full but the upper link is empty we have a definite edge
        if (false == data->upperLink.contains(edge))
        {
            data->jacobiType[edge] = 0;
            continue;
        }

        // If both the upper and lower link are not empty, we can have a regular or indefinite edge
        assert(data->lowerLink.contains(edge) && data->upperLink.contains(edge));

        // The list of the triangles adjacent to the edge in the star
        std::set<std::set<int>> adjacentTriangles;

        // Assemble the triangles
        for (const auto &v : linkVertices)
        {
            adjacentTriangles.insert({edge.first, edge.second, v});
        }

        // Compute the number of connected components in the link
        DisjointSet<std::set<int>> linkConnectedComponents;
        linkConnectedComponents.initialize(adjacentTriangles);

        for (const auto &t1 : adjacentTriangles)
        {
            for (const auto &t2 : adjacentTriangles)
            {
                if (data->connectedTriangles.contains({t1, t2}))
                {
                    linkConnectedComponents.union_setsTriangle(t1, t2);
                }

            }
        }
 
        data->jacobiType[edge] = linkConnectedComponents.countConnectedComponents();
    }

    std::map<int, int> typeCount;
    for (auto &[edge, jacobiType]: data->jacobiType)
    {
        if (false == typeCount.contains(jacobiType))
        {
            typeCount[jacobiType] = 0;
        }

        typeCount[jacobiType] += 1;

        //printf("The Jacobi type of edge (%d, %d) is %d.\n", edge.first, edge.second, jacobiType);
    }

    for (auto &[jacobiType, count]: typeCount)
    {
        //printf("There are %d edges of type %d.\n", count, jacobiType);
    }



    // Print upper/lower links to debug

    // otherwise the lower and upper link flip around
    //for (const std::vector<size_t> tet : data->tetrahedra)
    //{
        //// All pairs give you all six edges
        //for (int a = 0 ; a < 4 ; a++)
        //{
            //for (int b = a + 1 ; b < 4 ; b++)
            //{
                //int aIndex = tet[a];
                //int bIndex = tet[b];
                
                //// Make sure the vertices of the edge are in sorted order to have consistent orientation
                //if (aIndex > bIndex)
                //{
                    //std::swap(aIndex, bIndex);
                //}

                ////printf("The upper link of (%d, %d) : ", aIndex, bIndex);
                ////for (const int v: data->upperLink[std::pair<int, int>({aIndex, bIndex})]) 
                ////{
                    ////printf("%d ", v);
                ////}
                ////printf("\n");

                ////printf("The lower link of (%d, %d) : ", aIndex, bIndex);
                ////for (const int v: data->lowerLink[std::pair<int, int>({aIndex, bIndex})]) 
                ////{
                    ////printf("%d ", v);
                ////}
                ////printf("\n");
            //}
        //}
    //}
}


void ReebSpace::computeTriangleAdjacency(Data *data)
{
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

    Timer::start();
    // Insert all the segments into the arrangement using the range insert method
    CGAL::insert(data->arr, segments.begin(), segments.end());
    Timer::stop("Computed arrangement                   :");

    std::cout << std::endl << std::endl << "The arrangement size:\n"
        << "   |V| = " << data->arr.number_of_vertices()
        << ",  |E| = " << data->arr.number_of_edges()
        << ",  |F| = " << data->arr.number_of_faces() << std::endl;


    // Print out all the faces in the arrangement
    //std::cout << "Faces in the arrangement:" << std::endl;

    int counter = 0;
    for (auto f = data->arr.faces_begin(); f != data->arr.faces_end(); ++f) 
    {

        Arrangement_2::Face_const_handle a = f;
        data->arrangementFacesIdices[a] = counter;
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

void ReebSpace::computePreimageGraphs(Data *data)
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
    Halfedge_const_handle outerHalfEdge = *outerFace->holes_begin();

    // Make sure there is only one originating curve, something has gone wrong otherwise (edge overlap)
    assert(std::distance(data->arr.originating_curves_begin(outerHalfEdge), data->arr.originating_curves_end(outerHalfEdge)) == 1);

    // Extract the original curve
    //const Segment_2& segment = *data->arr.originating_curves_begin(outerHalfEdge);

    //printf("The half edge came from this edge:\n");
    //std::cout << segment.source() << " -> " << segment.target() << std::endl;

    //printf("The original indices are %d and %d", data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]);
    //printf("\n\n");

    // Starting from an outer half edge
    std::queue<Arrangement_2::Halfedge_const_handle> traversalQueue;
    traversalQueue.push(outerHalfEdge);

    // Which faces are visited
    std::set<Arrangement_2::Face_const_handle> visited;
    visited.insert(outerHalfEdge->face());

    // The disjoint set to track the connected components of the preimage graph
    data->preimageGraphs.resize(data->arrangementFacesIdices.size());

    // The number of connected components for each preimage graph (computed from the disjoint set)
    data->arrangementFiberComponents.resize(data->arrangementFacesIdices.size(), -1);

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

            // Type is [std::vector<std::set<int>>, std::vector<std::set<int>>]
            auto [minusTriangles, plusTriangles] = ReebSpace::getMinusPlusTriangles(currentHalfEdge, data);

            // Set the current preimage graph to be the preimage graph of the parent
            std::set<std::set<int>> preimageGraph;
            for (const auto &[t, id] : data->preimageGraphs[currentFaceID].data)
            {
                preimageGraph.insert(t);
            }


            // Add and remove triangles of the upper/lower link of the crossed edge
            for (const auto triangle: minusTriangles)
            {
                preimageGraph.erase(triangle);
            }

            for (const auto triangle: plusTriangles)
            {
                preimageGraph.insert(triangle);
            }


            //
            // Step 3. Compute disjointSets[twinFaceID]
            //

            data->preimageGraphs[twinFaceID].initialize(preimageGraph);

            // @TODO very inefficient squared running time in the size of the contour
            for (const auto &[t1, id1] : data->preimageGraphs[twinFaceID].data)
            {
                for (const auto &[t2, id2] : data->preimageGraphs[twinFaceID].data)
                {
                    if (t1 == t2)
                    {
                        continue;
                    }

                    //printf("First triangle: ");
                    //for (const auto v: t1)
                    //{
                        //printf("%d ", v);
                    //}

                    //printf("\nSecond triangle: ");
                    //for (const auto v: t2)
                    //{
                        //printf("%d ", v);
                    //}


                    std::pair<std::set<int>, std::set<int>> trianglePair({t1, t2});

                    if (data->connectedTriangles.find(trianglePair) != data->connectedTriangles.end()) 
                    {
                        //printf("Unioning\n");
                        data->preimageGraphs[twinFaceID].union_setsTriangle(t1, t2);
                    }

                    //printf("-----------\n");
                }
            }
            
            // Finaly make sure everyon points to their root
            data->preimageGraphs[twinFaceID].update();

            //printf("That preimage graph has %d connected components.\n", data->faceDisjointSets[twinFaceID].countConnectedComponents());
            data->arrangementFiberComponents[twinFaceID] = data->preimageGraphs[twinFaceID].countConnectedComponents();

            //printf("------------------------------------------------------------------------ \n");



            //
            // Push all the half edge around the twinFace, that lead to other faces
            //

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


void ReebSpace::computeReebSpace(Data *data)
{
    // Assemble all faces and their roots, the face is given as an int ID and same for the root (in it's respective DS)
    std::set<std::pair<int, int>> facesAndRoots;
    for (auto face = data->arr.faces_begin(); face != data->arr.faces_end(); ++face) 
    {
        // Skip the outer face
        if (face->is_unbounded()) { continue; }

        const int faceID = data->arrangementFacesIdices[face];

        // Make ure everyone is pointing to their root
        data->preimageGraphs[faceID].update();
        for (const int &r : data->preimageGraphs[faceID].getUniqueRoots())
        {
            facesAndRoots.insert({faceID, r});
        }
    }

    // Talls us which faces, root pairs are connected to each other
    std::set<std::pair<std::pair<int, int>, std::pair<int, int>>> connectedFacesAndRoots;



    // For evey face
    for (auto face = data->arr.faces_begin(); face != data->arr.faces_end(); ++face) 
    {
        // Skip the outer face
        if (face->is_unbounded()) { continue; }

        //
        // Traverse the neighbouring faces to collect them and the minusPlusTriangles
        //

        // List of adjacent faces
        std::vector<Face_const_handle> adjacentFaces;

        // Map from a face to it's plus and minus triangles as we go from face -> twinFace
        std::map<Face_const_handle, std::pair<std::vector<std::set<int>>, std::vector<std::set<int>>>> minusPlusTriangles;

        // Walk around the boundary of the face
        Arrangement_2::Ccb_halfedge_const_circulator start = face->outer_ccb();
        Arrangement_2::Ccb_halfedge_const_circulator curr = start;

        do {
            // Only for debug purposes
            //const Segment_2 &segment = *data->arr.originating_curves_begin(curr);
            //std::cout << "Half-edge   from: " << curr->source()->point() << " to " << curr->target()->point() << std::endl;
            //std::cout << "Source-edge from: " << segment.source() << " to " << segment.target() << std::endl;
            //printf("The original indices are %d and %d", data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]);
            //printf("\n\n");

            // Pull out the face and the minusPlusTriangles
            Arrangement_2::Halfedge_const_handle twinHalfEdge = curr->twin();
            Arrangement_2::Face_const_handle twinFace = twinHalfEdge->face();
            minusPlusTriangles[twinFace] = ReebSpace::getMinusPlusTriangles(curr, data);
            ++curr;
        } while (curr != start);


        //
        // Go over all adjacent twinFaces
        //
        for (const auto& [twinFace, mpt] : minusPlusTriangles) 
        { 
            // Get the plusMinus triangles for the current face
            const auto& [minusTriangles, plusTriangles] = mpt;

            // Face Graph - minusTriangles + plusTriangles = twinFace Graph
            const int faceID = data->arrangementFacesIdices[face];
            const int twinFaceID = data->arrangementFacesIdices[twinFace];

            // See which of he roots (connected components) are active (contain an active triangle)
            std::set<int> activeRootsFace;
            for (std::set<int> triangle : minusTriangles)
            {
                activeRootsFace.insert(data->preimageGraphs[faceID].findTriangle(triangle));
            }

            std::set<int> activeRootsTwinFace;
            for (std::set<int> triangle : plusTriangles)
            {
                activeRootsTwinFace.insert(data->preimageGraphs[twinFaceID].findTriangle(triangle));

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
                        connectedFacesAndRoots.insert({
                                {faceID, rFace},
                                {twinFaceID, rTwinFace}
                                });
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
                connectedFacesAndRoots.insert({
                        {faceID, triangleRootFace}, 
                        {twinFaceID, triangleRootTwinFace}
                        });
            }

        }
    }


    //
    // Compute the Reeb space
    //
    data->reebSpace.initialize(facesAndRoots);
    for (const auto &[faceRoot1, faceRoot2] : connectedFacesAndRoots)
    {
        data->reebSpace.union_setsTriangle(faceRoot1, faceRoot2);
    }

    // Path compression to make sure everyone is poiting to their root
    data->reebSpace.update();


    //
    // Set up colourIDs for each sheet
    //
    std::set<int> uniqueSheetIDs;

    for (const auto &[key, value] : data->reebSpace.data)
    {
        uniqueSheetIDs.insert(data->reebSpace.findTriangle(key));
    }
    int counter = 0;
    for (const auto &sheetID : uniqueSheetIDs)
    {
        data->sheetToColour[sheetID] = counter++;
    }
}
