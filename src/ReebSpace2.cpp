#include "./Timer.h"
#include "./ReebSpace2.h"
#include "src/CGALTypedefs.h"

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



    //Timer::start();

    //intersections = 0;
    //for (const auto &seg1 : regularSegments)
    //{
        //for (const auto &seg2 : singularSegments)
        //{
            //if (CGAL::do_intersect(seg1, seg2))
            //{
                //intersections++;
            //}
        //}
    //}

    //printf("There were %d intersections.\n", intersections);

    //Timer::stop("Computed red/blud intersetions n^2     :");
}
