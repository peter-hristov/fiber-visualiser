#include "./Arrangement.h"
#include "src/CustomArrangementTypes.h"

#include <algorithm>
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
    std::vector<std::vector<int>> segmentPoints(segments.size());

    for (int i = 0 ; i < segments.size() ; i++)
    {
        // Push the both endpoints together, they will be sorted later anyway
        segmentPoints[i].push_back(segments[i]->a->index);
        segmentPoints[i].push_back(segments[i]->b->index);

        for (int j = i + 1 ; j < segments.size() ; j++)
        {
            //printf("Checking segments %d and %d.\n", i, j);
            std::optional<Point> result = computeIntersection(segments[i], segments[j]);

            // Add the new point of intersection to the list of all arrangement points
            if (result) 
            {
                ConstPointHandler intersectionPoint = std::make_shared<Point>(result.value());

                //std::cout << "Intersection at (" << intersectionPoint->x << ", " << intersectionPoint->y << ")\n";

                // Set the index of the point
                intersectionPoint->index = data->customArrangementPoints.size();
                data->customArrangementPoints.push_back(intersectionPoint);

                segmentPoints[i].push_back(intersectionPoint->index);
                segmentPoints[j].push_back(intersectionPoint->index);
            } 
            // Segment was not intersected, keep going
            else 
            {
                //std::cout << "No intersection.\n";
            }
        }
    }

    //printf("\n\n----------------------------\n\n");


    for (int i = 0 ; i < segments.size() ; i++)
    {
        // Sort intersection points along the segment (ChatGPT code)
        const DataType dx = segments[i]->b->x - segments[i]->a->x;
        const DataType dy = segments[i]->b->y - segments[i]->a->y;
        const DataType len2 = dx * dx + dy * dy;

        std::sort(segmentPoints[i].begin(), segmentPoints[i].end(), [&](const int a, const int b) {
                PointHandler p1 = data->customArrangementPoints[a];
                PointHandler p2 = data->customArrangementPoints[b];

                DataType t1 = ((p1->x - segments[i]->a->x) * dx + (p1->y - segments[i]->a->y) * dy) / len2;
                DataType t2 = ((p2->x - segments[i]->a->x) * dx + (p2->y - segments[i]->a->y) * dy) / len2;
                return t1 < t2;
                });

        //std::cout << "There are this many segment points " << segmentPoints[i].size() << " in segment " << i << std::endl;
        for (int j = 0 ; j < segmentPoints[i].size() ; j++)
        {
            const int pointIndex = segmentPoints[i][j];
            //printf("Point %d (%.2f, %.2f).\n", 
                    //data->customArrangementPoints[pointIndex]->index, 
                    //data->customArrangementPoints[pointIndex]->x.get_d(), 
                    //data->customArrangementPoints[pointIndex]->y.get_d());

            ////printf("PointRational %d (%Zd/%Zd, %Zd/%Zd).\n", 
                    //data->customArrangementPoints[pointIndex]->index, 
                    //data->customArrangementPoints[pointIndex]->x.get_num().get_mpz_t(), 
                    //data->customArrangementPoints[pointIndex]->x.get_den().get_mpz_t(), 
                    //data->customArrangementPoints[pointIndex]->y.get_num().get_mpz_t(), 
                    //data->customArrangementPoints[pointIndex]->y.get_den().get_mpz_t());
        }
        //std::cout << "\n";

        // Add the parts of the segment
        for (int j = 1 ; j < segmentPoints[i].size() ; j++)
        {
            data->customArrangementSegments.push_back(
                    std::make_shared<Segment>(
                        data->customArrangementPoints[segmentPoints[i][j-1]], 
                        data->customArrangementPoints[segmentPoints[i][j]], 
                        data->customArrangementSegments.size()));
        }

    }

    //printf("In total there are %ld arrangement vertices and %ld arrangement segments.\n", data->customArrangementPoints.size(), data->customArrangementSegments.size());





}



std::optional<CustomArrangement::Point> CustomArrangement::computeIntersection(ConstSegmentHandler &s1, ConstSegmentHandler &s2)
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

    if (t > 0 && t < 1  && u > 0 && u < 1) 
    {
        const DataType ix = A->x + t * dx1;
        const DataType iy = A->y + t * dy1;
        return Point(ix, iy, -1);  // -1 since it's not from input, we'll set it later
    }

    return std::nullopt;
}
