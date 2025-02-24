#include "./ReebSpace.h"
#include "src/CGALTypedefs.h"

#include <cstddef>
#include <map>
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

Arrangement_2 ReebSpace::computeArrangement(Data *data) 
{

    // Map from EdgeSet -> Vertex Subset (Upper/Lower Link)
    std::map<std::pair<int, int>, std::set<int>> upperLink;
    std::map<std::pair<int, int>, std::set<int>> lowerLink;

    // Ideally we do want this with a proper data structure with a STAR, then there's no write conflits
    // Make sure to always keep the edge (u, v) such that u < v in index value,
    // otherwise the lower and upper link flip around
    for (int i = 0 ; i < data->tetrahedra.size() ; i++)
    {
        std::vector<size_t> tet = data->tetrahedra[i];

        // All pairs give you all six edges
        for (int j = 0 ; j < 4 ; j++)
        {
            for (int k = j + 1 ; k < 4 ; k++)
            {
                // Get the indices of the vertices for the edge
                int aIndex = tet[j];
                int bIndex = tet[k];

                // Make sure they are in a sorted order to have consistent orientation
                if (aIndex > bIndex)
                {
                    std::swap(aIndex, bIndex);
                }

                // Now we need the ohere two vertices to see which halflpane they are on
                std::set<int> usedVertices{aIndex, bIndex};

                // Search though the other two unused vertices
                for (int z = 0 ; z < 4 ; z++)
                {
                    // If the currect vertex is not one of two used to define the edge
                    if (usedVertices.find(tet[z]) == usedVertices.end()) 
                    {
                        // Define the two points that form the line
                        Point_2 a(data->vertexCoordinatesF[aIndex], data->vertexCoordinatesG[aIndex]);
                        Point_2 b(data->vertexCoordinatesF[bIndex], data->vertexCoordinatesG[bIndex]);

                        // Define the test point
                        Point_2 v(data->vertexCoordinatesF[tet[z]], data->vertexCoordinatesG[tet[z]]);  // Change this to test different locations

                        // Determine which half-plane r is in
                        CGAL::Orientation result = CGAL::orientation(a, b, v);

                        printf("Checking line (%d, %d) against vertex %ld\n", aIndex, bIndex, tet[z]);

                        std::cout << "a coords = " << a << std::endl;
                        std::cout << "b coords = " << b << std::endl;
                        std::cout << "v coords = " << v << std::endl;

                        if (result == CGAL::LEFT_TURN) {
                            //upperLink[]
                            upperLink[std::pair<int, int>({aIndex, bIndex})].insert(tet[z]);
                            std::cout << "Point r is in the LEFT half-plane.\n";
                        } else if (result == CGAL::RIGHT_TURN) {
                            lowerLink[std::pair<int, int>({aIndex, bIndex})].insert(tet[z]);
                            std::cout << "Point r is in the RIGHT half-plane.\n";
                        } else {
                            std::cout << "Point r is on the line.\n";
                            assert(false);
                        }
                    }
                }
            }
        }
    }

    // Example simple arrangement
    // |3--2|
    // |\  /|
    // | \/ |
    // | /\ |
    // |/  \|
    // |0--1|

    // Add in the points
    std::vector<Point_2> points;
    for (int i = 0 ; i < data->vertexCoordinatesF.size() ; i++)
    {
        const float u = data->vertexCoordinatesF[i];
        const float v = data->vertexCoordinatesG[i];
        points.push_back(Point_2(u, v));
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

        segments.push_back(Segment_2(points[edgeVector[0]], points[edgeVector[1]]));
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
    Arrangement_2 arr;

    // Insert all the segments into the arrangement using the range insert method
    CGAL::insert(arr, segments.begin(), segments.end());


    std::cout << "The arrangement size:\n"
        << "   |V| = " << arr.number_of_vertices()
        << ",  |E| = " << arr.number_of_edges()
        << ",  |F| = " << arr.number_of_faces() << std::endl;



    // Print out all the faces in the arrangement
    std::cout << "Faces in the arrangement:" << std::endl;
    for (auto f = arr.faces_begin(); f != arr.faces_end(); ++f) {
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


return arr;
}
