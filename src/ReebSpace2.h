#pragma once

#include "./TetMesh.h"
#include "./Arrangement.h"
#include "./CGALTypedefs.h"

class ReebSpace2
{
    public:


        //
        // Geometric computation
        //
        std::map<Halfedge_const_handle, std::set<int>> edgeRegionSegments;
        std::map<Halfedge_const_handle, std::set<int>> vertexRegionSegments;

        void computeEdgeRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeVertexRegionSegments(const TetMesh &tetMesh, Arrangement &singularArrangement);

        bool doSegmentEndpointsOverlap(const Segment_2 &s1, const Segment_2 &s2);
        bool ifSegmentInHalfEdgeRegion(Arrangement_2::Halfedge_around_vertex_const_circulator &halfEdgeCirculator, const Segment_2 &segment);
        Halfedge_const_handle getSegmentRegion(Vertex_const_handle &vertexHandle, const Segment_2 &segment);

        void unitTest(const TetMesh &tetMesh, Arrangement &singularArrangement, Arrangement &regularArrangement);




        std::map<Halfedge_const_handle, std::set<int>> edgeRegionMinusTriangles;
        std::map<Halfedge_const_handle, std::set<int>> edgeRegionPlusTriangles;
        std::map<Halfedge_const_handle, std::set<int>> vertexRegionMinusTriangles;
        std::map<Halfedge_const_handle, std::set<int>> vertexRegionPlusTriangles;

        void computeEdgeRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement);
        void computeVertexRegionMinusPlusTriangles(const TetMesh &tetMesh, Arrangement &singularArrangement);

        bool areHalfEdgeRegionMapsEqual(const std::map<Halfedge_const_handle, std::set<int>>& a, const std::map<Halfedge_const_handle, std::set<int>>& b);


};
