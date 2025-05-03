#include "./Arrangement.h"
#include "src/CustomArrangementTypes.h"

#include <algorithm>
#include <gmpxx.h>
#include <memory>
#include <vector>

void CustomArrangement::computeArrangementCustom(Data *data)
{
    //
    // Make an array of points,
    //

    // ASSUMPTION - these are already sorted lexicographically
    std::vector<Point> points(data->vertexCoordinatesF.size());

    for (int i = 0 ; i < data->vertexCoordinatesF.size() ; i++)
    {
        const float u = data->vertexCoordinatesF[i];
        const float v = data->vertexCoordinatesG[i];

        //const long uInt = static_cast<long>(std::round(u * 100.0));
        //const long vInt = static_cast<long>(std::round(v * 100.0));

        //points[i].x = mpq_class(uInt, 100);
        //points[i].y = mpq_class(vInt, 100);
        //points[i].index = i;

        points[i].x = u;
        points[i].y = v;
        points[i].index = i;
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

    std::vector<Segment> segments;
    for (const auto& edge : uniqueEdges) 
    {
        // Put in a vector for easy access
        std::vector<int> edgeVector(edge.first.begin(), edge.first.end());
        assert(edgeVector.size() == 2);

        // Swap if the first point is not smaller
        if (points[edgeVector[1]] < points[edgeVector[0]])
        {
            std::swap(edgeVector[0], edgeVector[1]);

        }

        segments.push_back(Segment(points[edgeVector[0]].index, points[edgeVector[1]].index, segments.size()));
    }

    //std::cout << "Printing all points...\n";

    for (int i = 0 ; i < points.size() ; i++)
    {
        //printf("Point %d (%.2f, %.2f).\n", points[i].index, points[i].x.get_d(), points[i].y.get_d());
        //std::cout << "Rational : " << points[i].x.get_num() << "/" << points[i].x.get_den();
        //std::cout << " " << points[i].y.get_num() << "/" << points[i].y.get_den() << "\n";
    }


    //std::cout << "\n\nPrinting all segments...\n";

    for (int i = 0 ; i < segments.size() ; i++)
    {
        //auto pointA = segments[i]->a;
        //auto pointB = segments[i]->b;
        //printf("Segment %d between point %d (%.2f, %.2f) and point %d (%.2f, %.2f).\n", segments[i]->index, pointA->index, pointA->x.get_d(), pointA->y.get_d(), pointB->index, pointB->x.get_d(), pointB->y.get_d());
    }


    // doubling up on memory here, but more efficient 
    std::set<Point> uniquePoints(points.begin(), points.end());
    data->customArrangementPoints = points;


    // All points in the segment, endpoints + intersection points
    // The pair is <pointId, \lambda \in [0,1]>
    std::vector<std::vector<std::pair<int, DataType>>> segmentPoints(segments.size());

    for (int i = 0 ; i < segments.size() ; i++)
    {
        // Push the both endpoints together, they will be sorted later anyway
        segmentPoints[i].push_back({data->customArrangementPoints[segments[i].aIndex].index, 0});
        segmentPoints[i].push_back({data->customArrangementPoints[segments[i].bIndex].index, 1});

        for (int j = i + 1 ; j < segments.size() ; j++)
        {
            // If the two segments have a common endpoint, skip
            if (segments[i].aIndex == segments[j].aIndex || segments[i].aIndex == segments[j].bIndex || segments[i].bIndex == segments[j].bIndex || segments[i].bIndex == segments[j].aIndex)
            {
                continue;
            }

            if (std::optional<std::tuple<Point, DataType, DataType>> result = computeIntersection(segments[i], segments[j], data->customArrangementPoints))
            {
                auto &[intersectionPoint, t, u] = *result;

                // Now we find the index of the point
                int pointIndex = -1;

                auto it = uniquePoints.find(intersectionPoint);

                // If the point has not already been added
                if (it == uniquePoints.end())
                {
                    pointIndex = uniquePoints.size();
                    intersectionPoint.index = pointIndex;
                    uniquePoints.insert(intersectionPoint);
                    data->customArrangementPoints.push_back(intersectionPoint);
                }
                else
                {
                    pointIndex = (*it).index;
                }

                //std::cout << "Intersection at (" << intersectionPoint->x << ", " << intersectionPoint->y << ")\n";

                // Set the index of the point and add it to the rest

                // Add the index of the point to the segments it intersects
                segmentPoints[i].push_back({pointIndex, t});
                segmentPoints[j].push_back({pointIndex, u});
            }
        }
    }

    //printf("\n\n----------------------------\n\n");

    for (int i = 0 ; i < segments.size() ; i++)
    {
        // Sort intersection points along the segment using their lambda value
        std::sort(segmentPoints[i].begin(), segmentPoints[i].end(), [&](const std::pair<int, DataType> &a, std::pair<int, DataType> &b) {
                return a.second < b.second;
                });

        //
        // Print points along the edge
        //
        //std::cout << "There are this many segment points " << segmentPoints[i].size() << " in segment " << i << std::endl;
        //for (int j = 0 ; j < segmentPoints[i].size() ; j++)
        //{
            //const int pointIndex = segmentPoints[i][j].first;

            //printf("Point %d (%.2f, %.2f).\n", 
                    //data->customArrangementPoints[pointIndex].index, 
                    //data->customArrangementPoints[pointIndex].x.get_d(), 
                    //data->customArrangementPoints[pointIndex].y.get_d());


            //std::cout << "Rational : " << data->customArrangementPoints[i].x.get_num() << "/" << data->customArrangementPoints[i].x.get_den();
            //std::cout << " " << data->customArrangementPoints[i].y.get_num() << "/" << data->customArrangementPoints[i].y.get_den() << "\n";
        //}
        //std::cout << "\n";

        // Add the parts of the segment
        for (int j = 1 ; j < segmentPoints[i].size() ; j++)
        {
            data->customArrangementSegments.push_back(
                    Segment(
                        data->customArrangementPoints[segmentPoints[i][j-1].first].index, 
                        data->customArrangementPoints[segmentPoints[i][j].first].index, 
                        data->customArrangementSegments.size()));
        }

    }

    //printf("In total there are %ld arrangement vertices and %ld arrangement segments.\n", data->customArrangementPoints.size(), data->customArrangementSegments.size());
}



std::optional<std::tuple<CustomArrangement::Point, CustomArrangement::DataType, CustomArrangement::DataType>> CustomArrangement::computeIntersection(const Segment &s1, const Segment &s2, const std::vector<Point> &segmentPoints)
{
    const Point &A = segmentPoints[s1.aIndex];
    const Point &B = segmentPoints[s1.bIndex];
    const Point &C = segmentPoints[s2.aIndex];
    const Point &D = segmentPoints[s2.bIndex];

    const DataType dx1 = B.x - A.x;
    const DataType dy1 = B.y - A.y;
    const DataType dx2 = D.x - C.x;
    const DataType dy2 = D.y - C.y;

    const CustomArrangement::DataType denom = dx1 * dy2 - dy1 * dx2;

    // Lines are parallel or collinear
    if (denom == 0) 
    {
        return std::nullopt;
    }

    const DataType dx = C.x - A.x;
    const DataType dy = C.y - A.y;

    const DataType t_num = dx * dy2 - dy * dx2;
    const DataType u_num = dx * dy1 - dy * dx1;

    const DataType t = t_num / denom;
    const DataType u = u_num / denom;

    // Should not happen, this means the segments share and endpoint
    assert(t != 0 && t != 1 && u != 0 && u != 1);

    if (t > 0 && t < 1  && u > 0 && u < 1) 
    {
        const DataType ix = A.x + t * dx1;
        const DataType iy = A.y + t * dy1;
        return std::make_tuple(Point(ix, iy, -1), t, u);
    }

    return std::nullopt;
}
