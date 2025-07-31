#include "./CGALTypedefs.h"

#include "./Timer.h"
#include "./ReebSpace2.h"
#include <CGAL/enum.h>

bool areHalfEdgeRegionMapsEqual(const std::map<Halfedge_const_handle, std::set<const Segment_2*>>& a, const std::map<Halfedge_const_handle, std::set<const Segment_2*>>& b)
{
    if (a.size() != b.size())
    {
        std::cerr << "Size is not equal.\n";
        return false;
    }

    for (const auto& [key, setA] : a)
    {
        auto itB = b.find(key);
        if (itB == b.end())
        {
            std::cerr << "Second map does not have a key.\n";
            return false;
        }

        const std::set<const Segment_2*>& setB = itB->second;
        if (setA.size() != setB.size())
        {
            std::cerr << "Sets for the same key do not match in size.\n";
            return false;
        }

        if (setA != setB)
        {
            std::cerr << "Set elements for the same key do not match.\n";
            return false;
        }
    }

    return true;
}

bool doSegmentEndpointsOverlap(const Segment_2 &s1, const Segment_2 &s2)
{
    return 
        s1.source() == s2.source() || s1.source() == s2.target() ||
        s1.target() == s2.source() || s1.target() == s2.target();
}

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
        //std::cout << "Return not right, right " << CGAL::orientation(o, a, b) << ", " << CGAL::orientation(o, b, c) << std::endl;
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
        //std::cout << "  Halfedge from " << current->source()->point() << " to " << current->target()->point() << "\n";

        if (ifSegmentInHalfEdgeRegion(current, segment))
        {
            //std::cout << "Found its rightful place.\n";
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


    std::map<Halfedge_const_handle, std::set<const Segment_2*>> halfEdgeEdgeRegionSegments;

    int intersections = 0;
    int falsePositives = 0;

    auto cb = [&](const Box& regularSegmentBox, const Box& singularSegmentBox) {
        
        const auto &regularSegmentHandle = *regularSegmentBox.handle();
        const auto &singularSegmentHandle  = *singularSegmentBox.handle();

        if (doSegmentEndpointsOverlap(regularSegmentHandle.seg, singularSegmentHandle.seg)) { return; }

        std::optional<std::variant<Point_2, Segment_2>> result = CGAL::intersection(regularSegmentHandle.seg, singularSegmentHandle.seg);


        if (result && std::holds_alternative<Point_2>(*result))
        {
            Point_2 ip = std::get<Point_2>(*result);
            intersections++;

            Halfedge_const_handle he = *(singularSegmentHandle.originatingHalfEdge); 
            singularArrangement.halfEdgePoints[he].push_back(ip);

            halfEdgeEdgeRegionSegments[he].insert(&(regularSegmentHandle.seg));

            //singularArrangement.arr.split_edge(*(singularSegmentHandle.originatingHalfEdge), ip);

            //const int aIndex = singularArrangement.arrangementPointIndices.at(singularSegmentHandle.seg.source());
            //const int bIndex = singularArrangement.arrangementPointIndices.at(singularSegmentHandle.seg.target());


            //intersectionPoints[{aIndex, bIndex}].push_back(ip);

            //std::cout << "Intersection point: " << ip << "\n";
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




    std::map<Halfedge_const_handle, std::set<const Segment_2*>> halfEdgeVertexRegionSegments;

    for (const MySegment_2 &mySegment : regularSegments)
    {
        const Segment_2 &segment = mySegment.seg;

        const int aIndex = singularArrangement.arrangementPointIndices.at(segment.source());
        const int bIndex = singularArrangement.arrangementPointIndices.at(segment.target());
        //printf("-----------------------------------------------------  At segment [%d, %d].\n", aIndex, bIndex);

        const auto itSource = singularArrangement.arrangementPointHandles.find(segment.source());
        if (itSource != singularArrangement.arrangementPointHandles.end())
        {
            //std::cout << "The source point is : " << itSource->second->point() << " | " << singularArrangement.arrangementPointIndices.at(itSource->second->point()) << std::endl;
            Halfedge_const_handle regionHalfEdgeHandle = getSegmentRegion(itSource->second, segment);
            halfEdgeVertexRegionSegments[regionHalfEdgeHandle].insert(&segment);
        }


        const auto itTarget = singularArrangement.arrangementPointHandles.find(segment.target());
        if (itTarget != singularArrangement.arrangementPointHandles.end())
        {
            //std::cout << "\n\nThe target point is : " << itTarget->second->point() << " | " <<  singularArrangement.arrangementPointIndices.at(itTarget->second->point()) << std::endl;
            Halfedge_const_handle regionHalfEdgeHandle = getSegmentRegion(itTarget->second, segment);
            halfEdgeVertexRegionSegments[regionHalfEdgeHandle].insert(&segment);
        }
    }
















    // Make sure we have compute all intersections correctly
    Timer::start();
    std::map<Halfedge_const_handle, std::set<const Segment_2*>> halfEdgeEdgeRegionSegmentsUnitTest;

    for (const auto &myRegularSegment : regularSegments)
    {
        for (const auto &mySingularSegment : singularSegments)
        {
            if (doSegmentEndpointsOverlap(myRegularSegment.seg, mySingularSegment.seg)) { continue; }

            std::optional<std::variant<Point_2, Segment_2>> result = CGAL::intersection(myRegularSegment.seg, mySingularSegment.seg);

            if (result && std::holds_alternative<Point_2>(*result))
            {
                Point_2 ip = std::get<Point_2>(*result);
                intersections++;

                Halfedge_const_handle he = *(mySingularSegment.originatingHalfEdge); 

                halfEdgeEdgeRegionSegmentsUnitTest[he].insert(&(myRegularSegment.seg));
            }
        }
    }

    assert(areHalfEdgeRegionMapsEqual(halfEdgeEdgeRegionSegments, halfEdgeEdgeRegionSegmentsUnitTest));
    Timer::stop("Segment intersection unit test         :");

    Timer::start();


    // halfEdgeVertexRegionSegments;


    std::map<Halfedge_const_handle, std::set<const Segment_2*>> halfEdgeVertexRegionSegmentsUnitTest;

    // Iterate over the halfedges of every vertex
    for (auto v = regularArrangement.arr.vertices_begin(); v != regularArrangement.arr.vertices_end(); ++v)
    {
        if (v->is_isolated() || false == regularArrangement.arrangementPointIndices.contains(v->point())) 
        {
            continue;
        }

        //std::cout << "At vertex " << v->point() << std::endl;

        Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
        first = curr = v->incident_halfedges();

        bool singularEdgeFound = false;

        // Go to the first singular edge you find.
        do 
        {
            // If the edge is singular
            const Segment_2 &segment = *regularArrangement.arr.originating_curves_begin(curr);

            const int aIndex = regularArrangement.arrangementPointIndices.at(segment.source());
            const int bIndex = regularArrangement.arrangementPointIndices.at(segment.target());
            const int edgeType = tetMesh.edgeSingularTypes.at({aIndex, bIndex});
            if  (edgeType != 1)
            {
                first  = curr;
                singularEdgeFound = true;
                break;
            }

        } while (--curr != first);

        if (false == singularEdgeFound) { continue; }


        //std::cout << "Starting half-edge at " << curr->source()->point() << " to " << curr->target()->point() << "\n";

        Arrangement_2::Halfedge_const_handle currentHalfEdge;
        Arrangement_2::Halfedge_const_handle currentSingularHalfEdge;
        do 
        {
            // If the edge is singular
            const Segment_2 &segment = *regularArrangement.arr.originating_curves_begin(curr);


            
            const int aIndex = regularArrangement.arrangementPointIndices.at(segment.source());
            const int bIndex = regularArrangement.arrangementPointIndices.at(segment.target());
            const int edgeType = tetMesh.edgeSingularTypes.at({aIndex, bIndex});
            if  (edgeType != 1)
            {
                currentHalfEdge = curr;

                bool singularHalfEdgeFound = false;

                for (auto he = singularArrangement.arr.halfedges_begin(); he != singularArrangement.arr.halfedges_end(); ++he)
                {
                    const Segment_2 &singularSegment = *singularArrangement.arr.originating_curves_begin(he);

                    if (he->target()->point() == v->point() && segment == singularSegment)
                    {

                        currentSingularHalfEdge = he;
                        singularHalfEdgeFound = true;
                    }
                }

                assert(singularHalfEdgeFound);





                //std::cout << "---- Singular region start at " << curr->source()->point() << " to " << curr->target()->point() << "\n";
                //std::cout << "---- Singular region start at singular " << currentSingularHalfEdge->source()->point() << " to " << currentSingularHalfEdge->target()->point() << "\n";
                //std::cout << std::endl;
            }
            else
            {

                //std::cout << "---------- Adding regular edge " << curr->source()->point() << " to " << curr->target()->point() << "\n";

                const Segment_2 *nonTempSegment;

                bool nonTempSegmentFound = false;
                for (const MySegment_2 &mySegment : regularSegments)
                {
                    if (mySegment.seg.source() == segment.source() && mySegment.seg.target() == segment.target())
                    {
                        nonTempSegment = &mySegment.seg;
                        nonTempSegmentFound = true;
                    }
                }

                assert (nonTempSegmentFound);







                halfEdgeVertexRegionSegmentsUnitTest[currentSingularHalfEdge].insert(nonTempSegment);
                // The current half edge is equal to which half-edge in the singular arrangement?


            }


        } while (--curr != first);
    }


    assert(areHalfEdgeRegionMapsEqual(halfEdgeVertexRegionSegments, halfEdgeVertexRegionSegmentsUnitTest));
    //areHalfEdgeRegionMapsEqual(halfEdgeVertexRegionSegments, halfEdgeVertexRegionSegmentsUnitTest);



    //printf("Printing the  singular vertex regions ... \n");
    //for (const auto& [halfEdge, segments] : halfEdgeVertexRegionSegments)
    //{
        //std::cout << "At half edge " << halfEdge->source()->point() << " to " << halfEdge->target()->point() << "\n";
        //for (const auto &segment : segments)
        //{
            //std::cout << "----" << segment->source() << " to " << segment->target() << "\n";

        //}

    //}

    //printf("\n\nPrinting the singular vertex regions from the unit test... \n");
    //for (const auto& [halfEdge, segments] : halfEdgeVertexRegionSegmentsUnitTest)
    //{
        //std::cout << "At half edge " << halfEdge->source()->point() << " to " << halfEdge->target()->point() << "\n";
        //for (const auto &segment : segments)
        //{
            //std::cout << "----" << segment->source() << " to " << segment->target() << "\n";

        //}

    //}



    Timer::stop("Vertex regions unit test               :");
    Timer::stop("Segment intersection unit test         :");











}
