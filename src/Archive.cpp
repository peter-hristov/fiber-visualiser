
//
// Good test for the arrangement
//
#include "./ReebSpace.h"

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

Arrangement_2 ReebSpace::computeArrangement(Data *data) {


    // Example simple arrangement
    // |3--2|
    // |\  /|
    // | \/ |
    // | /\ |
    // |/  \|
    // |0--1|

    std::vector<Point_2> points = {
        Point_2(0,0),
        Point_2(1,0),
        Point_2(1,1),
        Point_2(0,1)
    };

    std::vector<Segment_2> segments = {
        // Edges of the square
        Segment_2(points[0], points[1]),
        Segment_2(points[1], points[2]),
        Segment_2(points[2], points[3]),
        Segment_2(points[3], points[0]),

        // Diaganals
        Segment_2(points[0], points[2]),
        Segment_2(points[1], points[3]),

    };

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


    for(auto e = arr.edges_begin(); e != arr.edges_end(); ++e) {
        std::cout << "Edge: " << e->curve() << "\n";
        std::cout << "  Originating curve(s):\n";
        for(auto oc = arr.originating_curves_begin(e);
                oc != arr.originating_curves_end(e); ++oc)
        {
            std::cout << "    " << *oc << "\n";
        }
        std::cout << std::endl;
    }



    return arr;





}
