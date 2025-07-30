#include "./CGALTypedefs.h"

#include "./Timer.h"
#include "./ReebSpace2.h"
#include <CGAL/enum.h>



bool ifSegmentInHalfEdgeRegion(Arrangement_2::Halfedge_around_vertex_const_circulator &halfEdgeCirculator, const Segment_2 &segment)
{
    assert(halfEdgeCirculator != nullptr); 

    const auto halfEdge = halfEdgeCirculator;

    auto halfEdgeCirculatorPrevious = halfEdgeCirculator;
    const auto halfEdgeNext = --halfEdgeCirculatorPrevious;


    // The image of the vertex is o, the endpoints of the half-edges are a and b and the other endpoints of the segment is b
    // See bellow for pictures
    const Point_2 &o = halfEdge->target()->point();
    const Point_2 &a = halfEdge->source()->point();
    const Point_2 &b = segment.source() == o ? segment.target() : segment.source();
    const Point_2 &c = halfEdgeNext->source()->point();

    //std::cout << "o = " << o << "\na = " << a << "\nb = " << b << "\nc = " << c << std::endl;

    // Case 1. 
    //
    //      a
    //  b   |
    //   \  |
    //    \ |
    //     \|
    // c----o
    if (CGAL::orientation(o, a, c) == CGAL::LEFT_TURN)
    {
        //printf("Returning left, left");
        //std::cout << "Return left, left " << CGAL::orientation(o, a, b) << ", " << CGAL::orientation(o, b, c) << std::endl;

        return CGAL::orientation(o, a, b) == CGAL::LEFT_TURN && CGAL::orientation(o, b, c) == CGAL::LEFT_TURN;
    }
    // Case 2. This is the only case where c \notin the region a-b, we need the negation of this case.
    //
    // a
    // |   b
    // |  /
    // | /
    // |/
    // o----c
    else if (CGAL::orientation(o, a, c) == CGAL::RIGHT_TURN)
    {
        //printf("Returning not right, right");
        return ! (CGAL::orientation(o, a, b) == CGAL::RIGHT_TURN && CGAL::orientation(o, b, c) == CGAL::RIGHT_TURN);
    }
    else
    {
        assert(false);
    }
}

Halfedge_const_handle getSegmentRegion(Vertex_const_handle &vertexHandle, const Segment_2 &segment)
{
    assert(segment.source() == vertexHandle->point() || segment.target() == vertexHandle->point());
    assert(false == vertexHandle->is_isolated());

    const auto first = vertexHandle->incident_halfedges();
    auto current = first;

    // Iterate in a CCW fashion to keep things consistent (faces are iterated in a CCW fashion too)
    do 
    {
        std::cout << "  Halfedge from " << current->source()->point() << " to " << current->target()->point() << "\n";

        if (ifSegmentInHalfEdgeRegion(current, segment))
        {
            std::cout << "Found its rightful place.\n";
            return current;
        }

    } while (--current != first);

    // Something has gone terribly wrong if we are here.
    assert(false);
}



void ReebSpace2::compute(const TetMesh &tetMesh, Arrangement &regularArrangement, Arrangement &singularArrangement)
{
    Timer::start();





    //std::vector<Segment_2> regularSegments;
    //regularSegments.reserve(tetMesh.regularEdgesNumber);

    //std::vector<Segment_2> singularSegments;
    //singularSegments.reserve(tetMesh.singularEdgesNumber);

    //std::vector<Box> regularBoxes; 
    //regularBoxes.reserve(tetMesh.regularEdgesNumber);

    //std::vector<Box> singularBoxes;
    //singularBoxes.reserve(tetMesh.singularEdgesNumber);

    //for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    //{
        //if (type == 1)
        //{
            //regularSegments.emplace_back(arrangement.arrangementPoints[edge[0]], arrangement.arrangementPoints[edge[1]]);
            //regularBoxes.emplace_back(regularSegments.back().bbox(), &regularSegments.back());
        //}
        //else
        //{
            //singularSegments.emplace_back(arrangement.arrangementPoints[edge[0]], arrangement.arrangementPoints[edge[1]]);
            //singularBoxes.emplace_back(singularSegments.back().bbox(), &singularSegments.back());
        //}
    //}













    std::vector<MySegment_2> singularSegments;
    singularSegments.reserve(singularArrangement.arr.number_of_edges());

    std::vector<Box> singularBoxes;
    singularBoxes.reserve(singularArrangement.arr.number_of_edges());

    for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
    {
        if (he < he->twin())
        {
            singularSegments.emplace_back(Segment_2(he->source()->point(), he->target()->point()), he);
            singularBoxes.emplace_back(singularSegments.back().seg.bbox(), &singularSegments.back());
        }
    }

    std::vector<MySegment_2> regularSegments;
    regularSegments.reserve(tetMesh.regularEdgesNumber);

    std::vector<Box> regularBoxes; 
    regularBoxes.reserve(tetMesh.regularEdgesNumber);

    for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    {
        if (type == 1)
        {
            regularSegments.emplace_back(Segment_2(singularArrangement.arrangementPoints[edge[0]], singularArrangement.arrangementPoints[edge[1]]), std::nullopt);
            regularBoxes.emplace_back(regularSegments.back().seg.bbox(), &regularSegments.back());
        }
    }


    int intersections = 0;
    int falsePositives = 0;

    auto cb = [&](const Box& regularSegmentBox, const Box& singularSegmentBox) {
        
        const auto &regularSegmentHandle = *regularSegmentBox.handle();
        const auto &singularSegmentHandle  = *singularSegmentBox.handle();


        std::optional<std::variant<Point_2, Segment_2>> result = CGAL::intersection(regularSegmentHandle.seg, singularSegmentHandle.seg);


        if (result)
        {
            if (std::holds_alternative<Point_2>(*result))
            {
                Point_2 ip = std::get<Point_2>(*result);
                intersections++;

                if (ip != singularSegmentHandle.seg.source() && ip != singularSegmentHandle.seg.target())
                {
                    Halfedge_const_handle he = *(singularSegmentHandle.originatingHalfEdge); 
                    singularArrangement.halfEdgePoints[he].push_back(ip);
                }

                //singularArrangement.arr.split_edge(*(singularSegmentHandle.originatingHalfEdge), ip);

                //const int aIndex = singularArrangement.arrangementPointIndices.at(singularSegmentHandle.seg.source());
                //const int bIndex = singularArrangement.arrangementPointIndices.at(singularSegmentHandle.seg.target());


                //intersectionPoints[{aIndex, bIndex}].push_back(ip);

                //std::cout << "Intersection point: " << ip << "\n";
            }
            else
            {
                std::cout << "Segments overlap (result is a segment).\n";
            }
        }
        else
        {
            //std::cout << "No intersection.\n";
            falsePositives++;
        }


        //if (CGAL::do_intersect(regularSegmentHandle, singularSegmentHandle)) 
        //{
            //intersections++;
        //}
        //else
        //{
            //falsePositives++;
        //}
    };

    CGAL::box_intersection_d(regularBoxes.begin(), regularBoxes.end(),
                              singularBoxes.begin(), singularBoxes.end(),
                              cb);




    Timer::stop("Computed red/blud intersetions         :");
    printf("There were %d intersections with %d false positives.\n", intersections, falsePositives);




    std::map<Halfedge_const_handle, std::vector<const Segment_2*>> halfEdgeRegionSegments;

    for (const MySegment_2 &mySegment : regularSegments)
    {
        const Segment_2 &segment = mySegment.seg;

        const int aIndex = singularArrangement.arrangementPointIndices.at(segment.source());
        const int bIndex = singularArrangement.arrangementPointIndices.at(segment.target());
        printf("-----------------------------------------------------  At segment [%d, %d].\n", aIndex, bIndex);

        const auto itSource = singularArrangement.arrangementPointHandles.find(segment.source());
        if (itSource != singularArrangement.arrangementPointHandles.end())
        {
            std::cout << "The source point is : " << itSource->second->point() << " | " << singularArrangement.arrangementPointIndices.at(itSource->second->point()) << std::endl;
            Halfedge_const_handle regionHalfEdgeHandle = getSegmentRegion(itSource->second, segment);
            halfEdgeRegionSegments[regionHalfEdgeHandle].push_back(&segment);
        }


        const auto itTarget = singularArrangement.arrangementPointHandles.find(segment.target());
        if (itTarget != singularArrangement.arrangementPointHandles.end())
        {
            std::cout << "\n\nThe target point is : " << itTarget->second->point() << " | " <<  singularArrangement.arrangementPointIndices.at(itTarget->second->point()) << std::endl;
            Halfedge_const_handle regionHalfEdgeHandle = getSegmentRegion(itTarget->second, segment);
            halfEdgeRegionSegments[regionHalfEdgeHandle].push_back(&segment);
        }
    }
}
