#include "./Arrangement.h"
#include "src/CustomArrangementTypes.h"

#include <algorithm>
#include <memory>
#include <vector>

void CustomArrangement::computeArrangementCustom(Data *data)
{
    //
    // Make an array of points,
    //

    // ASSUMPTION - these are already sorted lexicographically
    std::vector<PointHandler> points(data->vertexCoordinatesF.size());

    for (int i = 0 ; i < data->vertexCoordinatesF.size() ; i++)
    {
        const double u = data->vertexCoordinatesF[i];
        const double v = data->vertexCoordinatesG[i];

        points[i] = std::make_shared<Point>(Point(u, v, i));
    }


    //
    // Make an array of segments
    //

    // Get the unique edges
    std::map<std::set<size_t>, bool> uniqueEdges;
    for (int i = 0 ; i < data->tetrahedra.size() ; i++)
    {
        // Bottom face
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][0], data->tetrahedra[i][1]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][1], data->tetrahedra[i][2]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][2], data->tetrahedra[i][0]})] = true;

        // Connect bottom face to top vertex
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][0], data->tetrahedra[i][3]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][1], data->tetrahedra[i][3]})] = true;
        uniqueEdges[std::set<size_t>({data->tetrahedra[i][2], data->tetrahedra[i][3]})] = true;
    }

    std::vector<SegmentHandler> segments;
    for (const auto& edge : uniqueEdges) 
    {
        // Put in a vector for easy access
        std::vector<int> edgeVector(edge.first.begin(), edge.first.end());
        assert(edgeVector.size() == 2);

        // Swap if the first point is not smaller
        if (*points[edgeVector[1]] < *points[edgeVector[0]])
        {
            std::swap(edgeVector[0], edgeVector[1]);

        }

        segments.push_back(std::make_shared<Segment>(points[edgeVector[0]], points[edgeVector[1]], segments.size()));
    }

    //std::cout << "Printing all points...\n";

    for (int i = 0 ; i < points.size() ; i++)
    {
        //printf("Point %d (%.2f, %.2f).\n", points[i]->index, points[i]->x.get_d(), points[i]->y.get_d());
    }


    //std::cout << "\n\nPrinting all segments...\n";

    for (int i = 0 ; i < segments.size() ; i++)
    {
        auto pointA = segments[i]->a;
        auto pointB = segments[i]->b;
        //printf("Segment %d between point %d (%.2f, %.2f) and point %d (%.2f, %.2f).\n", segments[i]->index, pointA->index, pointA->x.get_d(), pointA->y.get_d(), pointB->index, pointB->x.get_d(), pointB->y.get_d());
    }














    data->customArrangementPoints = points;

    // All points in the segment, endpoints + intersection points
    // The pair is <pointId, \lambda \in [0,1]>
    std::vector<std::vector<std::pair<int, DataType>>> segmentPoints(segments.size());

    for (int i = 0 ; i < segments.size() ; i++)
    {
        // Push the both endpoints together, they will be sorted later anyway
        segmentPoints[i].push_back({segments[i]->a->index, 0});
        segmentPoints[i].push_back({segments[i]->b->index, 1});

        for (int j = i + 1 ; j < segments.size() ; j++)
        {
            // If the two segments have a common endpoint, skip
            if (segments[i]->a == segments[j]->a || segments[i]->a == segments[j]->b || segments[i]->b == segments[j]->b || segments[i]->b == segments[j]->a)
            {
                continue;
            }

            if (std::optional<std::tuple<PointHandler, DataType, DataType>> result = computeIntersection(segments[i], segments[j])) 
            {
                auto &[intersectionPoint, t, u] = *result;

                //std::cout << "Intersection at (" << intersectionPoint->x << ", " << intersectionPoint->y << ")\n";

                // Set the index of the point and add it to the rest
                intersectionPoint->index = data->customArrangementPoints.size();
                data->customArrangementPoints.push_back(intersectionPoint);

                // Add the index of the point to the segments it intersects
                segmentPoints[i].push_back({intersectionPoint->index, t});
                segmentPoints[j].push_back({intersectionPoint->index, u});
            }
        }
    }

    //printf("\n\n----------------------------\n\n");

    for (int i = 0 ; i < segments.size() ; i++)
    {
        // Sort intersection points along the segment using their lambda value
        std::sort(segmentPoints[i].begin(), segmentPoints[i].end(), [&](const std::pair<int, DataType> &a, std::pair<int, DataType> b) {
                return a.second < b.second;
                });

        //
        // Print points along the edge
        //
        //std::cout << "There are this many segment points " << segmentPoints[i].size() << " in segment " << i << std::endl;
        //for (int j = 0 ; j < segmentPoints[i].size() ; j++)
        //{
            //const int pointIndex = segmentPoints[i][j];
            ////printf("Point %d (%.2f, %.2f).\n", 
                    ////data->customArrangementPoints[pointIndex]->index, 
                    ////data->customArrangementPoints[pointIndex]->x.get_d(), 
                    ////data->customArrangementPoints[pointIndex]->y.get_d());

            //////printf("PointRational %d (%Zd/%Zd, %Zd/%Zd).\n", 
                    ////data->customArrangementPoints[pointIndex]->index, 
                    ////data->customArrangementPoints[pointIndex]->x.get_num().get_mpz_t(), 
                    ////data->customArrangementPoints[pointIndex]->x.get_den().get_mpz_t(), 
                    ////data->customArrangementPoints[pointIndex]->y.get_num().get_mpz_t(), 
                    ////data->customArrangementPoints[pointIndex]->y.get_den().get_mpz_t());
        //}
        //std::cout << "\n";

        // Add the parts of the segment
        for (int j = 1 ; j < segmentPoints[i].size() ; j++)
        {
            data->customArrangementSegments.push_back(
                    std::make_shared<Segment>(
                        data->customArrangementPoints[segmentPoints[i][j-1].first], 
                        data->customArrangementPoints[segmentPoints[i][j].first], 
                        data->customArrangementSegments.size()));
        }

    }

    //printf("In total there are %ld arrangement vertices and %ld arrangement segments.\n", data->customArrangementPoints.size(), data->customArrangementSegments.size());
}



std::optional<std::tuple<CustomArrangement::PointHandler, CustomArrangement::DataType, CustomArrangement::DataType>> CustomArrangement::computeIntersection(ConstSegmentHandler &s1, ConstSegmentHandler &s2)
{
    const ConstPointHandler A = s1->a;
    const ConstPointHandler B = s1->b;
    const ConstPointHandler C = s2->a;
    const ConstPointHandler D = s2->b;

    const DataType dx1 = B->x - A->x;
    const DataType dy1 = B->y - A->y;
    const DataType dx2 = D->x - C->x;
    const DataType dy2 = D->y - C->y;

    const CustomArrangement::DataType denom = dx1 * dy2 - dy1 * dx2;

    // Lines are parallel or collinear
    if (denom == 0) 
    {
        return std::nullopt;
    }

    const DataType dx = C->x - A->x;
    const DataType dy = C->y - A->y;

    const DataType t_num = dx * dy2 - dy * dx2;
    const DataType u_num = dx * dy1 - dy * dx1;

    const DataType t = t_num / denom;
    const DataType u = u_num / denom;

    // Should not happen, this means the segments share and endpoint
    assert(t != 0 && t != 1 && u != 0 && u != 1);

    if (t > 0 && t < 1  && u > 0 && u < 1) 
    {
        const DataType ix = A->x + t * dx1;
        const DataType iy = A->y + t * dy1;
        return std::make_tuple(std::make_shared<Point>(ix, iy, -1), t, u);
    }

    return std::nullopt;
}
