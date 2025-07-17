#pragma once

#include "./TetMesh.h"
#include "./Arrangement.h"

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


class ReebSpace
{
    public:
        // The connected components of the preimage graph for each face in the arrangement
        std::vector<DisjointSet<int>> preimageGraphs;

        // The actual Reeb space, map from a connected component of a preimage graph to a sheet.
        DisjointSet<std::pair<int, int>> reebSpace;

        // For each faceID we have a number of triangle seeds, which are given in fiberComponentId
        std::vector<std::vector<std::pair<int, int>>> fiberSeeds;

        std::unordered_map<int, CartesianPolygon_2> sheetPolygon;

        std::map<int, double> sheetArea;

        std::unordered_set<int> incompleteSheets;

        // Maps the IDs of the Reeb space sheets to consequitive integers, useful for colouring things
        std::unordered_map<int, int> sheetToColour;

        // Jacobi type computing from the Reeb space algorithm, for comparison with the Jacobi type computed from the mesh
        std::unordered_map<std::pair<int, int>, int, MyHash<std::pair<int, int>>> jacobiType;


        // Computes the correspondence between two faces given a half edge between them
        void computeTwinFacePreimageGraph(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle &);

        // Compute the preimage graphs Gi for each cell in the arrangement
        void computePreimageGraphs(const TetMesh &tetMesh, const Arrangement &arrangement, const bool);

        // This function is not used now, it's incorporated in computePreimageGraphs, otherwise we can't just keep seeds.
        // Computes the correspondence between two faces given a half edge between them
        void determineCorrespondence(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle &);

        // Computing corresponded graph H
        void computeCorrespondenceGraph(const TetMesh &tetMesh, const Arrangement &arrangement);

        // Compute the Reeb space from all the preimage graphs
        void computeReebSpacePostprocess(const TetMesh &tetMesh, const Arrangement &arrangement);

        //void countIntersectionsTypes(TetMesh &tetMesh, Arrangement &arrangement);


        //
        // Helper functions
        //
        std::pair<std::set<int>, std::set<int>> getMinusPlusTrianglesIndex(const TetMesh &tetMesh, const Arrangement &arrangement, const Arrangement_2::Halfedge_const_handle currentHalfEdge);
        void testTraverseArrangement(const Arrangement &arrangement);
};
