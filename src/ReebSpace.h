#pragma once

#include "./Data.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Arr_segment_traits_2<K> Traits_2;
typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef K::Point_2 Point_2;
typedef K::Line_2 Line_2;
typedef K::Segment_2 Segment_2;

class ReebSpace
{
    public:
        ReebSpace(Data *data) {


            // Example simlpe arrangement
            // |----|
            // |\  /|
            // | \/ |
            // | /\ |
            // |/  \|
            // |----|

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

            for (int i = 0 ; i < segments.size(); i++)
            {
                CGAL::insert(arr, segments[i]);

            }


            std::cout << "The arrangement size:\n"
            << "   |V| = " << arr.number_of_vertices()
            << ",  |E| = " << arr.number_of_edges()
            << ",  |F| = " << arr.number_of_faces() << std::endl;



            // Print out all the faces in the arrangement
            std::cout << "Faces in the arrangement:" << std::endl;
            for (auto f = arr.faces_begin(); f != arr.faces_end(); ++f) {
                if (f->is_unbounded()) {
                    std::cout << "Unbounded face" << std::endl;
                } else {
                    std::cout << "Bounded face" << std::endl;
                }


                for (auto e = f->outer_ccbs_begin(); e != f->outer_ccbs_end(); ++e) {
                    std::cout << "E" << std::endl;

                    
                    
                }
            }
        }
};
