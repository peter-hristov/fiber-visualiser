#include "./ReebSpace.h"
#include "src/CGALTypedefs.h"

#include <cstddef>
#include <map>
#include <queue>
#include <set>

template <typename Arrangement>
void print_ccb(typename Arrangement::Ccb_halfedge_const_circulator circ) 
{

    std::cout << "Start = (" << circ->source()->point() << ")" << std::endl;
    typename Arrangement::Ccb_halfedge_const_circulator curr = circ;
    do {
        typename Arrangement::Halfedge_const_handle e = curr;
        // Get the twin half-edge and its adjacent face
        typename Arrangement::Halfedge_const_handle twin = e->twin();

        std::cout << "   [" << e->curve() << "]   " << "(" << e->target()->point() << ")" << std::endl;
    } while (++curr != circ);
    std::cout << std::endl;
}

bool ReebSpace::isUpperLinkEdgeVertex(int aIndex, int bIndex, int vIndex, Data *data)
{
    // Make sure the vertices of the edge are in sorted order to have consistent orientation
    if (aIndex > bIndex)
    {
        std::swap(aIndex, bIndex);
    }

    // Define the two points that form the line
    Point_2 a(data->vertexCoordinatesF[aIndex], data->vertexCoordinatesG[aIndex]);
    Point_2 b(data->vertexCoordinatesF[bIndex], data->vertexCoordinatesG[bIndex]);

    // Define the test point
    Point_2 v(data->vertexCoordinatesF[vIndex], data->vertexCoordinatesG[vIndex]);  // Change this to test different locations

    // Determine which half-plane r is in
    CGAL::Orientation result = CGAL::orientation(a, b, v);

    printf("Checking line (%d, %d) against vertex %d\n", aIndex, bIndex, vIndex);

    std::cout << "a coords = " << a << std::endl;
    std::cout << "b coords = " << b << std::endl;
    std::cout << "v coords = " << v << std::endl;


    // Upper link = left
    if (result == CGAL::LEFT_TURN) {
        return true;
        std::cout << "Point r is in the LEFT half-plane.\n";
    // Lower link = right
    } else if (result == CGAL::RIGHT_TURN) {
        return false;
        std::cout << "Point r is in the RIGHT half-plane.\n";
    // This should not happen for generic maps
    } else {
        std::cout << "Point r is on the line.\n";
        assert(false);
    }

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

                // Search though the other two unused vertices
                for (int v = 0 ; v < 4 ; v++)
                {
                    int vIndex = tet[v];

                    // If the currect vertex is not one of two used to define the edge
                    if (vIndex != aIndex && vIndex != bIndex) 
                    {
                        bool isUpperLink = isUpperLinkEdgeVertex(aIndex, bIndex, vIndex, data);

                        if (true == isUpperLink ) {
                            data->upperLink[std::pair<int, int>({aIndex, bIndex})].insert(vIndex);
                        } else {
                            data->lowerLink[std::pair<int, int>({aIndex, bIndex})].insert(vIndex);
                        }
                    }
                }
            }
        }
    }

    // Print upper/lower links to debug

    // otherwise the lower and upper link flip around
    for (const std::vector<size_t> tet : data->tetrahedra)
    {
        // All pairs give you all six edges
        for (int a = 0 ; a < 4 ; a++)
        {
            for (int b = a + 1 ; b < 4 ; b++)
            {
                int aIndex = tet[a];
                int bIndex = tet[b];
                
                // Make sure the vertices of the edge are in sorted order to have consistent orientation
                if (aIndex > bIndex)
                {
                    std::swap(aIndex, bIndex);
                }

                printf("The upper link of (%d, %d) : ", aIndex, bIndex);
                for (const int v: data->upperLink[std::pair<int, int>({aIndex, bIndex})]) 
                {
                    printf("%d ", v);
                }
                printf("\n");

                printf("The lower link of (%d, %d) : ", aIndex, bIndex);
                for (const int v: data->lowerLink[std::pair<int, int>({aIndex, bIndex})]) 
                {
                    printf("%d ", v);
                }
                printf("\n");
            }
        }
    }
}


void ReebSpace::computeArrangement(Data *data) 
{


    // Example simple arrangement
    // |3--2|
    // |\  /|
    // | \/ |
    // | /\ |
    // |/  \|
    // |0--1|

    // Add in the points
    for (int i = 0 ; i < data->vertexCoordinatesF.size() ; i++)
    {
        const float u = data->vertexCoordinatesF[i];
        const float v = data->vertexCoordinatesG[i];
        Point_2 point(u, v);

        data->arrangementPoints.push_back(point);
        data->arrangementPointsIdices[point] = i;
    };

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

    // Add the unique edges as setments to the arrangement
    std::vector<Segment_2> segments;
    for (const auto& edge : uniqueEdges) 
    {
        // Put in a vector for easy access
        std::vector<int> edgeVector(edge.first.begin(), edge.first.end());
        assert(edgeVector.size() == 2);

        segments.push_back(Segment_2(data->arrangementPoints[edgeVector[0]], data->arrangementPoints[edgeVector[1]]));
        std::cout << "Adding edge " << edgeVector[0] << " - " << edgeVector[1] << std::endl;
    }


    //std::vector<Segment_2> segments = {
        //// Edges of the square
        //Segment_2(points[0], points[1]),
        //Segment_2(points[1], points[2]),
        //Segment_2(points[2], points[3]),
        //Segment_2(points[3], points[0]),

        //// Diaganals
        //Segment_2(points[0], points[2]),
        //Segment_2(points[1], points[3]),

    //};

    // Create the arrangement to store the lines

    // Insert all the segments into the arrangement using the range insert method
    CGAL::insert(data->arr, segments.begin(), segments.end());


    std::cout << "The arrangement size:\n"
        << "   |V| = " << data->arr.number_of_vertices()
        << ",  |E| = " << data->arr.number_of_edges()
        << ",  |F| = " << data->arr.number_of_faces() << std::endl;



    // Print out all the faces in the arrangement
    std::cout << "Faces in the arrangement:" << std::endl;

    int counter = 0;
    for (auto f = data->arr.faces_begin(); f != data->arr.faces_end(); ++f) 
    {

        Arrangement_2::Face_const_handle a = f;
        data->arrangementFacesIdices[a] = counter;
        counter++;

        std::cout << data->arrangementFacesIdices[a] << std::endl;

        if (f->is_unbounded()) {
            std::cout << "Unbounded face" << std::endl;
            continue;
        }

        std::cout << "Bounded face with " << f->number_of_holes() << " holes " << std::endl;


        std::cout << "inner     = " << f->number_of_inner_ccbs() << std::endl;
        std::cout << "outer     = " << f->number_of_outer_ccbs() << std::endl;
        std::cout << "holes     = " << f->number_of_holes() << std::endl;
        std::cout << "isolated  = " << f->number_of_isolated_vertices() << std::endl;

        print_ccb<Arrangement_2>(f->outer_ccb());


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

void ReebSpace::BFS(Data *data)
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

    // Sanity check
    assert(outerFace->number_of_inner_ccbs() == 1);
    assert(outerFace->number_of_outer_ccbs() == 0);
    assert(outerFace->number_of_holes() == 1);
    assert(outerFace->number_of_isolated_vertices() == 0);


    // Loop the edges of the boundary
    Arrangement_2::Ccb_halfedge_const_circulator inner_ccb = *outerFace->holes_begin();
    Arrangement_2::Ccb_halfedge_const_circulator he = inner_ccb;

    // Iterate over the halfedges forming the inner CCB
    do {
        std::cout << he->source()->point() << " -> ";
        std::cout << he->target()->point() << std::endl;

        // Get the originating curve of the half edge
        // Maybe need a function for this
        const Segment_2& segment = *data->arr.originating_curves_begin(he);
        std::cout << data->arrangementPointsIdices[segment.source()] << ", ";
        std::cout << data->arrangementPointsIdices[segment.target()] << std::endl;
    } while (++he != inner_ccb);

    printf("\n\n");

    // Get the first edge
    printf("This is the first half edge:\n");
    Halfedge_const_handle outerHalfEdge = *outerFace->holes_begin();
    std::cout << outerHalfEdge->source()->point() << " -> " << outerHalfEdge->target()->point() << std::endl;

    printf("\n\n");


    // Make sure there is only one originating curve, something has gone wrong otherwise (edge overlap)
    assert(std::distance(data->arr.originating_curves_begin(outerHalfEdge), data->arr.originating_curves_end(outerHalfEdge)) == 1);

    // Extract the original curve
    const Segment_2& segment = *data->arr.originating_curves_begin(outerHalfEdge);

    printf("The half edge came from this edge:\n");
    std::cout << segment.source() << " -> " << segment.target() << std::endl;

    printf("The original indices are %d and %d", data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]);
    printf("\n\n");

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


        // Get the twin
        Arrangement_2::Halfedge_const_handle twin = currentHalfEdge->twin();
        Arrangement_2::Face_const_handle twinFace = twin->face();

        // If we have never visited this face, then we have never visited any of the half edges.
        if (visited.find(twinFace) == visited.end())
        {
            printf("NEW FACE ------------------------------------------ %d \n", data->arrangementFacesIdices[twinFace]);
            visited.insert(twinFace);

            // We should only have onbounded face, the outside face.
            assert(false == twinFace->is_unbounded());

            Arrangement_2::Ccb_halfedge_const_circulator start = twinFace->outer_ccb();
            Arrangement_2::Ccb_halfedge_const_circulator curr = start;

            do {
                traversalQueue.push(curr);


                // Make sure there is only one originating curve (sanity check)
                assert(std::distance(data->arr.originating_curves_begin(outerHalfEdge), data->arr.originating_curves_end(outerHalfEdge)) == 1);
                const Segment_2 &segment = *data->arr.originating_curves_begin(curr);

                std::cout << "Half-edge   from: " << curr->source()->point() << " to " << curr->target()->point() << std::endl;
                std::cout << "Source-edge from: " << segment.source() << " to " << segment.target() << std::endl;

                printf("The original indices are %d and %d", data->arrangementPointsIdices[segment.source()], data->arrangementPointsIdices[segment.target()]);
                printf("\n\n");

                ++curr;
            } while (curr != start);
        }
    }

    //do {
        //typename Arrangement::Halfedge_const_handle e = curr;
        //// Get the twin half-edge and its adjacent face
        //typename Arrangement::Halfedge_const_handle twin = e->twin();

        //std::cout << "   [" << e->curve() << "]   " << "(" << e->target()->point() << ")" << std::endl;
    //} while (++curr != circ);
}
