#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_segment_traits_2.h>

#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Arr_segment_traits_2<K> Traits_2;
typedef CGAL::Arrangement_with_history_2<Traits_2> Arrangement_2;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;

//typedef CGAL::Arrangement_2<Traits_2> Arrangement_2;
typedef Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
typedef Arrangement_2::Face_const_handle Face_const_handle;

typedef Arrangement_2::Halfedge_handle Halfedge_handle;
typedef Arrangement_2::Face_handle Face_handle;
typedef Arrangement_2::Vertex_iterator Vertex_iterator;
