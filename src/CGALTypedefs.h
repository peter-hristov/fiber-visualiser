#pragma once

// Correct and slow
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//typedef CGAL::Exact_predicates_exact_constructions_kernel K;


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Lazy_exact_nt.h>
typedef CGAL::Lazy_exact_nt<double> NT;
typedef CGAL::Simple_cartesian<NT> K;

//#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<double> K;

//#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_segment_traits_2.h>

#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>


typedef CGAL::Arr_segment_traits_2<K> Traits_2;
typedef CGAL::Arrangement_with_history_2<Traits_2> Arrangement_2;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;

typedef CGAL::Arr_segment_2<K> Curve_2;

typedef Arrangement_2::Halfedge_handle Halfedge_handle;
typedef Arrangement_2::Face_handle Face_handle;
typedef Arrangement_2::Vertex_iterator Vertex_iterator;

// Const iterators and handles
typedef Arrangement_2::Face_const_iterator Face_const_iterator;

//typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
typedef Arrangement_2::Face_const_handle Face_const_handle;

// Used to find which face a point is in 
#include <CGAL/Arr_trapezoid_ric_point_location.h>
typedef CGAL::Arr_trapezoid_ric_point_location<Arrangement_2> Point_location;
