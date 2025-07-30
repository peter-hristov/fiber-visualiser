#pragma once

#include "./CGALTypedefs.h"
#include "./TetMesh.h"

#include <map>
#include <vector>

// This is mostly a wrapper for the CGAL arrangement class
class Arrangement
{
    public:

    // Singular stuff
    enum class SegmentMode
    {
        UseSingularSegments,
        UseAllSegments
    };

    // Regular points along each half-edge
    std::map<Halfedge_const_handle, std::vector<Point_2>> halfEdgePoints;

    // Previous stuff

    Arrangement_2 arr;

    // The points of the line segments that define the arrangement, does not include new intersection points
    std::vector<Point_2> arrangementPoints;

    Face_const_handle getActiveFace(const std::array<float, 2>);

    // The inverse map of arrangementPoints, returns the index of a point
    std::map<Point_2, int> arrangementPointIndices;

    std::map<Point_2, Vertex_const_handle> arrangementPointHandles;

    // The integer ID of every face in the arrangement
    std::map<Arrangement_2::Face_const_handle, int> arrangementFacesIdices;
    std::vector<Arrangement_2::Face_const_handle> arrangementIndexToFace;

    // Which faces are adjacent to a singular edge
    std::set<Arrangement_2::Face_const_handle> singularFaces;


    std::unique_ptr<Point_location> pl;  // nullptr by default

    // We only add data->arr, the rest of data is unchanged
    void computeArrangement(const TetMesh&, const SegmentMode&);
    void computePointLocationDataStructure();


};
