//
// CGAL includes 
//
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel K;

#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>

typedef CGAL::Arr_segment_traits_2<K> Traits_2;
typedef CGAL::Arrangement_with_history_2<Traits_2> Arrangement_2;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;

typedef CGAL::Arr_segment_2<K> Curve_2;

typedef Arrangement_2::Halfedge_handle Halfedge_handle;
typedef Arrangement_2::Face_handle Face_handle;
typedef Arrangement_2::Vertex_const_handle Vertex_const_handle;
typedef Arrangement_2::Vertex_iterator Vertex_iterator;
typedef Arrangement_2::Face_const_iterator Face_const_iterator;
typedef Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
typedef Arrangement_2::Face_const_handle Face_const_handle;




//
// Data structure - tetrahedral mesh with two float values defined at each vertex, each edge of the mesh is mapped to a segment in the plane. 
// We want the arrangement of those mapped segments. They have some structure (no endpoints of a segment has degree one).
//

class Data
{
  public:
    std::vector<std::vector<size_t>> tetrahedra;
    std::vector<float> vertexCoordinatesF;
    std::vector<float> vertexCoordinatesG;

    Arrangement_2 arr;

    // Question - is there a way to put a integer index in each face, so I don't have to use maps?

    // The points of the line segments that define the arrangement, does not include new intersection points
    std::vector<Point_2> arrangementPoints;
    // The inverse map of arrangementPoints, returns the index of a point
    std::map<Point_2, int> arrangementPointsIdices;

    // List of the faces in the arrangement
    std::vector<Arrangement_2::Face_const_handle> arrangementIndexToFace;
    // The inverse map of arrangementIndexToFace, returns the index of a point
    std::map<Arrangement_2::Face_const_handle, int> arrangementFacesIdices;
};




//
// Compute the arrangement
//

data->arrangementPoints.resize(data->vertexCoordinatesF.size());
for (int i = 0 ; i < data->vertexCoordinatesF.size() ; i++)
{
    const float u = data->vertexCoordinatesF[i];
    const float v = data->vertexCoordinatesG[i];
    const Point_2 point(u, v);

    data->arrangementPoints[i] = point;
    data->arrangementPointsIdices[point] = i;
};

// Make sure you don't add duplicate edge to the arrangement, each tet has 6 edges
std::set<std::set<size_t>> uniqueSegments;
for (int i = 0 ; i < data->tetrahedra.size() ; i++)
{
    // Bottom face
    uniqueSegments.insert({data->tetrahedra[i][0], data->tetrahedra[i][1]});
    uniqueSegments.insert({data->tetrahedra[i][1], data->tetrahedra[i][2]});
    uniqueSegments.insert({data->tetrahedra[i][2], data->tetrahedra[i][0]});

    // Connect bottom face to top vertex
    uniqueSegments.insert({data->tetrahedra[i][0], data->tetrahedra[i][3]});
    uniqueSegments.insert({data->tetrahedra[i][1], data->tetrahedra[i][3]});
    uniqueSegments.insert({data->tetrahedra[i][2], data->tetrahedra[i][3]});
}

// Add the unique edges as setments to the arrangement
std::vector<Segment_2> segments;
for (const std::set<size_t> &segment : uniqueSegments) 
{
    const size_t endpointA = *segment.begin();
    const size_t endpointB = *std::prev(segment.end());
    segments.push_back(Segment_2(data->arrangementPoints[endpointA], data->arrangementPoints[endpointB]));
}

CGAL::insert(data->arr, segments.begin(), segments.end());





//
// Using the arrangement - BFS over the incidence graph (traverse the faces)
//


// Find the unbounded face to use a starting point for the graph search (BFS)
Face_const_handle outerFace;

// Iterate over all faces and find the unbounded one
for (Face_const_iterator fit = data->arr.faces_begin(); fit != data->arr.faces_end(); ++fit) 
{
    if (fit->is_unbounded()) 
    {
        outerFace = fit;
        break;
    }
}

// Sanity check, make sure the outer face is simple
assert(outerFace->number_of_inner_ccbs() == 1);
assert(outerFace->number_of_outer_ccbs() == 0);
assert(outerFace->number_of_holes() == 1);
assert(outerFace->number_of_isolated_vertices() == 0);

std::queue<Face_const_handle> traversalQueue;
std::vector<float> visited(data->arrangementFacesIdices.size(), false);

traversalQueue.push(outerFace);
visited[data->arrangementFacesIdices[outerFace]] = true;

while (false == traversalQueue.empty())
{
    // Pop an half edge out
    Face_const_handle currentFace = traversalQueue.front();
    traversalQueue.pop();

    // Get ids of the current face and the twin face
    int currentFaceID = data->arrangementFacesIdices[currentFace];

    // Iterate over all the neighbours of the current face
    Halfedge_const_handle start;

    // The unbounded face (the starting one) has a different way of addressing its neighbours
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
        if (false == visited[twinFaceID])
        {
            traversalQueue.push(twinFace);
            visited[twinFaceID] = true;

            // Compute the preimage graph of this unvisited face
            ReebSpace::computeTwinFacePreimageGraph(data, curr);
        }


        ++curr;
    } while (curr != start);
}

// Used in function ReebSpace::computeTwinFacePreimageGraph
const Segment_2 &segment = *data->arr.originating_curves_begin(currentHalfEdge);
const int aIndex = data->arrangementPointsIdices[segment.source()];
const int bIndex = data->arrangementPointsIdices[segment.target()];
