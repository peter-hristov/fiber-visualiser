#pragma once

#include <map>
#include <vector>

#include "./CGALTypedefs.h"
#include "./TetMesh.h"

// This is mostly a wrapper for the CGAL arrangement class
class Arrangement
{
    public:

    Arrangement_2 arr;

    // The points of the line segments that define the arrangement, does not include new intersection points
    std::vector<Point_2> arrangementPoints;

    // The inverse map of arrangementPoints, returns the index of a point
    std::map<Point_2, int> arrangementPointsIdices;

    // The integer ID of every face in the arrangement
    std::map<Arrangement_2::Face_const_handle, int> arrangementFacesIdices;
    std::vector<Arrangement_2::Face_const_handle> arrangementIndexToFace;

    // We only add data->arr, the rest of data is unchanged
    void computeArrangement(const TetMesh&);
};
