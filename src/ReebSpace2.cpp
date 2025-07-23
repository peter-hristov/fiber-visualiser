#include "./Timer.h"
#include "./ReebSpace2.h"

void ReebSpace2::compute(const TetMesh &tetMesh, const Arrangement &arrangement)
{
    std::vector<Segment_2> regularSegments;
    regularSegments.reserve(tetMesh.regularEdgesNumber);

    std::vector<Segment_2> singularSegments;
    singularSegments.reserve(tetMesh.singularEdgesNumber);

    std::vector<Box> regularBoxes; 
    regularBoxes.reserve(tetMesh.regularEdgesNumber);

    std::vector<Box> singularBoxes;
    singularBoxes.reserve(tetMesh.singularEdgesNumber);

    for (const auto &[edge, type] : tetMesh.edgeSingularTypes) 
    {
        if (type == 1)
        {
            regularSegments.emplace_back(arrangement.arrangementPoints[edge[0]], arrangement.arrangementPoints[edge[1]]);
            regularBoxes.emplace_back(regularSegments.back().bbox(), &regularSegments.back());
        }
        else
        {
            singularSegments.emplace_back(arrangement.arrangementPoints[edge[0]], arrangement.arrangementPoints[edge[1]]);
            singularBoxes.emplace_back(singularSegments.back().bbox(), &singularSegments.back());
        }
    }

    int intersections = 0;
    int falsePositives = 0;

    auto cb = [&](const Box& regularSegment, const Box& singularSegment) {
        
        const auto &regularSegmentHandle = *regularSegment.handle();
        const auto &singularSegmentHandle  = *singularSegment.handle();

        if (CGAL::do_intersect(regularSegmentHandle, singularSegmentHandle)) 
        {
            intersections++;
        }
        else
        {
            falsePositives++;
        }
    };

    CGAL::box_intersection_d(regularBoxes.begin(), regularBoxes.end(),
                              singularBoxes.begin(), singularBoxes.end(),
                              cb);


    Timer::stop("Computed red/blud intersetions         :");
    printf("There were %d intersections with %d false positives.\n", intersections, falsePositives);
}
