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

        // We only add upperLink/lowerLink to data, the rest of data is unchagned
        void computeUpperLowerLink(Data *);
        
        // We only add upperLink/lowerLink to data, the rest of data is unchagned
        void computeTriangleAdjacency(Data *);

        // We only add data->arr, the rest of data is unchanged
        void computeArrangement(Data *data);

        // Computes the correspondence between two faces given a half edge between them
        void computeTwinFacePreimageGraph(Data *data, Arrangement_2::Halfedge_const_handle &);

        // Compute the preimage graphs Gi for each cell in the arrangement
        void computePreimageGraphs(Data *, const bool);
        void computePreimageGraphsSingular(Data *, const bool);

        // Computes the correspondence between two faces given a half edge between them
        void determineCorrespondence(Data *data, Arrangement_2::Halfedge_const_handle &);

        // Computing corresponded graph H
        void computeCorrespondenceGraph(Data *);

        // Compute the Reeb space from all the preimage graphs
        void computeReebSpacePostprocess(Data *);

        void countIntersectionsTypes(Data *);


        //
        // Helper functions
        //

        // Give the edge (aIndex, bIndex), is the vertex vIndex from its link in the upper and lower link of the edge
        // We assume that aIndex < bIndex for consistent orientation.
        bool isUpperLinkEdgeVertex(int aIndex, int bIndex, int vIndex, Data *data);

        // Get the upper/lower link with the orientation of the half edge with respect to the original edge.
        //std::pair<std::vector<std::set<int>>, std::vector<std::set<int>>> getMinusPlusTriangles(Arrangement_2::Halfedge_const_handle currentHalfEdge, Data *data);

        std::pair<std::set<int>, std::set<int>> getMinusPlusTrianglesIndex(Arrangement_2::Halfedge_const_handle currentHalfEdge, Data *data);

        void testTraverseArrangement(Data *data);



};
