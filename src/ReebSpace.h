#pragma once

#include "./Data.h"

#include <vector>

#include "CGALTypedefs.h"



// What do I need from the arrangement?
// Done - Compute the arrangement of a number of line segments
// Done - Traverse the half-edge arrangement data struture (use ccb and twin)
// Done - How to see where the segment came from? Use with_history and originating curve
// Set up IDs for the faces (use a map for now, store ID field in data later)
//
// Get a running BFS
// Lower and upper star for each edge (loop all tets, maybe we don't need a big datastruture)
//
//
//
// Later
// - Can you tell the algo which segments have the same point? Does that help?
// - How do we use rational numbers for precision? Robustness.


namespace ReebSpace
{
        bool isUpperLinkEdgeVertex(int aIndex, int bIndex, int vIndex, Data *data);

        // We only add upperLink/lowerLink to data, the rest of data is unchagned
        void computeUpperLowerLink(Data *);
        // We only add data->arr, the rest of data is unchanged
        void computeArrangement(Data *data);
        void BFS(Data *);
};
